#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

#include "CommonFunctions/Pandora.h"
#include "CommonFunctions/Scatters.h"
#include "CommonFunctions/Corrections.h"
#include "CommonFunctions/Region.h"
#include "CommonFunctions/Types.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "SignatureTools/SignatureToolBase.h"
#include "SignatureTools/DecayVertexProvider.h"

#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "TDatabasePDG.h"

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <chrono>

class SignatureIntegrityFilter : public art::EDFilter 
{
public:
    explicit SignatureIntegrityFilter(fhicl::ParameterSet const &pset);

    SignatureIntegrityFilter(SignatureIntegrityFilter const &) = delete;
    SignatureIntegrityFilter(SignatureIntegrityFilter &&) = delete;
    SignatureIntegrityFilter &operator=(SignatureIntegrityFilter const &) = delete;
    SignatureIntegrityFilter &operator=(SignatureIntegrityFilter &&) = delete;

    bool filter(art::Event &e) override;

private:
    art::InputTag _HitProducer, _MCPproducer, _MCTproducer, _BacktrackTag;

    int _lower_hit_thresh, _upper_hit_thresh;
    double _displaced_thresh;
    int _channel_proximity;
    double _length_thresh;

    std::string _bad_channel_file;
    bool _veto_bad_channels, _write_output_file;

    calo::CalorimetryAlg* _calo_alg;
    std::vector<std::unique_ptr<::signature::SignatureToolBase>> _signatureToolsVec;

    std::vector<bool> _bad_channel_mask;
    const geo::GeometryCore* _geo;

    double _fid_x_start, _fid_y_start, _fid_z_start;
    double _fid_x_end, _fid_y_end, _fid_z_end;

    double _angular_alignment_limit;

    std::string _output_file;

    void initialiseBadChannelMask();

    bool identifySignatures(art::Event& e, signature::SignatureCollection& sig_coll);
    bool areSignaturesValid(art::Event& e, const signature::SignatureCollection& sig_coll);
    bool validateDecayVertices(art::Event& e, const signature::SignatureCollection& sig_coll);
    bool validateFinalStateParticles(art::Event& e, const signature::SignatureCollection& sig_coll);

    void getFinalStateParticles(art::Event const& e, std::vector<simb::MCParticle>& fsp);
    void getNuVertex(art::Event const& e, TVector3& nu_vtx);
    void getSimulationHits(art::Event const& e, std::vector<art::Ptr<recob::Hit>>& mc_hits);

    bool checkSignatureActivelyBound(const art::Event& e, const signature::SignatureCollection& sig_coll);
    bool checkPositionActivelyBound(const double x[3]);

    bool checkSignatureSensitivity(const art::Event& e, const signature::SignatureCollection& sig_coll); 
    bool checkPositionSensitivity(const TVector3& pos);
    bool checkAngularAlignment(const art::Event& e, const signature::SignatureCollection& sig_coll);

    bool checkGeometricLength(const art::Event& e, const signature::SignatureCollection& sig_coll);
};

SignatureIntegrityFilter::SignatureIntegrityFilter(fhicl::ParameterSet const &pset)
    : EDFilter{pset}
    , _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    , _lower_hit_thresh{pset.get<int>("LowerHitThreshold", 2000)}
    , _upper_hit_thresh{pset.get<int>("UpperHitThreshold", 6000)}
    , _displaced_thresh{pset.get<double>("DisplacedThreshold", 1.0)}
    , _channel_proximity{pset.get<int>("ChannelProximity", 5)}
    , _length_thresh{pset.get<double>("LengthThreshold", 5)}
    , _bad_channel_file{pset.get<std::string>("BadChannelFile", "badchannels.txt")}
    , _veto_bad_channels{pset.get<bool>("VetoBadChannels", true)}
    , _fid_x_start{pset.get<double>("FiducialXStart", 10.0)}
    , _fid_y_start{pset.get<double>("FiducialYStart", 15.0)}
    , _fid_z_start{pset.get<double>("FiducialZStart", 10.0)}
    , _fid_x_end{pset.get<double>("FiducialXEnd", 10.0)}
    , _fid_y_end{pset.get<double>("FiducialYEnd", 15.0)}
    , _fid_z_end{pset.get<double>("FiducialZEnd", 50.0)}
    , _output_file{pset.get<std::string>("OutputFile", "run_subrun_event_log.txt")}
    , _write_output_file{pset.get<bool>("WriteOutputFile", false)}
    , _angular_alignment_limit{pset.get<double>("AngularAlignmentLimit", 0.1)}
{
    _calo_alg = new calo::CalorimetryAlg(pset.get<fhicl::ParameterSet>("CaloAlg"));

    const fhicl::ParameterSet &tool_psets = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (auto const &tool_pset_label : tool_psets.get_pset_names())
    {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_label);
        _signatureToolsVec.push_back(art::make_tool<::signature::SignatureToolBase>(tool_pset));
    };

    _geo = art::ServiceHandle<geo::Geometry>()->provider();
    size_t num_channels = _geo->Nchannels();
    _bad_channel_mask.resize(num_channels, false);

    if (_veto_bad_channels) {
        this->initialiseBadChannelMask();
    }
}

void SignatureIntegrityFilter::initialiseBadChannelMask()
{
    if (!_bad_channel_file.empty()) {
        cet::search_path sp("FW_SEARCH_PATH");
        std::string fullname;
        sp.find_file(_bad_channel_file, fullname);
        if (fullname.empty()) {
            throw cet::exception("SignatureIntegrityFilter") << "Bad channel file not found: " << _bad_channel_file;
        }

        std::ifstream inFile(fullname, std::ios::in);
        std::string line;
        while (std::getline(inFile, line)) {
            if (line.find("#") != std::string::npos) continue;
            std::istringstream ss(line);
            int ch1, ch2;
            ss >> ch1;
            if (!(ss >> ch2)) ch2 = ch1;
            for (int i = ch1; i <= ch2; ++i) {
                _bad_channel_mask[i] = true;
            }
        }
        std::cout << "Loaded bad channels from: " << fullname << std::endl;
    }
}

bool SignatureIntegrityFilter::filter(art::Event &e) 
{
    std::vector<signature::Signature> sig_coll;
    if (!this->identifySignatures(e, sig_coll))
        return false;

    if (!this->areSignaturesValid(e, sig_coll))
        return false;

    std::vector<art::Ptr<recob::Hit>> sim_hits;
    this->getSimulationHits(e, sim_hits);
    if (sim_hits.size() < static_cast<size_t>(_lower_hit_thresh) || sim_hits.size() > static_cast<size_t>(_upper_hit_thresh))
        return false;

    if (!this->validateDecayVertices(e, sig_coll))
        return false;

    if (!this->validateFinalStateParticles(e, sig_coll))
        return false;

    if (_write_output_file) 
    {
        int evt = e.event();
        int sub = e.subRun();
        int run = e.run();

        std::ofstream file(_output_file, std::ios::app);
        if (!file.is_open()) 
            throw cet::exception("SignatureIntegrityFilter") << "Unable to open file: " << _output_file;

        file << run << " " << sub << " " << evt << "\n";

        if (!file.good()) 
            throw cet::exception("SignatureIntegrityFilter") << "Error writing to file: " << _output_file;

        file.close();
    }

    return true; 
}

bool SignatureIntegrityFilter::identifySignatures(art::Event& e, const signature::SignatureCollection& sig_coll)
{
    for (auto &signatureTool : _signatureToolsVec) 
    {
        if (!signatureTool->identifySignalParticles(e, sig_coll))
            return false;
    }

    return true;
}

bool SignatureIntegrityFilter::areSignaturesValid(art::Event &e, const signature::SignatureCollection& sig_coll)
{
    return this->checkGeometricLength(e, sig_coll) &&
           this->checkSignatureActivelyBound(e, sig_coll) &&
           this->checkSignatureSensitivity(e, sig_coll) &&
           this->checkAngularAlignment(e, sig_coll);
}

bool SignatureIntegrityFilter::validateDecayVertices(art::Event &e, const signature::SignatureCollection& sig_coll)
{
    TVector3 nu_vtx;
    this->getNuVertex(e, nu_vtx);

    for (const auto &signatureTool : _signatureToolsVec) 
    {
        auto *decayTool = dynamic_cast<DecayVertexProvider*>(signatureTool.get());
        if (decayTool) 
        {
            std::optional<TVector3> decay_vtx_opt = decayTool->getDecayVertex(e);
            if (decay_vtx_opt) 
            {
                double sep = (nu_vtx - *decay_vtx_opt).Mag();
                if (sep < _displaced_thresh)
                    return false;
                break;
            }
        }
    }

    return true;
}

bool SignatureIntegrityFilter::validateFinalStateParticles(art::Event &e, const signature::SignatureCollection& sig_coll)
{
    std::vector<simb::MCParticle> fsp;
    this->getFinalStateParticles(e, fsp);

    for (const auto &mcp : fsp) 
    {
        if (std::any_of(sig_coll.begin(), sig_coll.end(), [&](const auto &signature) {
                return signature.trckid == mcp.TrackId();
            })) {
            continue;
        }

        TParticlePDG *p_info = TDatabasePDG::Instance()->GetParticle(mcp.PdgCode());
        if (!p_info)
            continue;

        if (std::string(p_info->ParticleClass()) == "meson")
            return false;
    }

    return true;
}

void SignatureIntegrityFilter::getFinalStateParticles(art::Event const& e, std::vector<simb::MCParticle>& fsp)
{
    auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    fsp.clear();
    for (const auto &mcp : *mcp_h) 
    {
        if (mcp.Process() == "primary")
            fsp.push_back(mcp);
    }
}

void SignatureIntegrityFilter::getNuVertex(art::Event const& e, TVector3& nu_vtx)
{
    auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
    const simb::MCTruth& mct = mct_h->at(0);
    if (!mct.NeutrinoSet()) 
        return;

    auto const &neutrino = mct.GetNeutrino();
    auto const &nu = neutrino.Nu();

    nu_vtx.SetXYZ(nu.Vx(), nu.Vy(), nu.Vz());
}

void SignatureIntegrityFilter::getSimulationHits(art::Event const& e, std::vector<art::Ptr<recob::Hit>>& mc_hits)
{
    art::Handle<std::vector<recob::Hit>> hit_handle;
    if (e.getByLabel(_HitProducer, hit_handle))
    {
        std::vector<art::Ptr<recob::Hit>> evt_hits;
        art::fill_ptr_vector(evt_hits, hit_handle);
        auto mcp_bkth_assoc = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hit_handle, e, _BacktrackTag);

        for (const auto& hit : evt_hits)
        {
            if (_veto_bad_channels && _bad_channel_mask[hit->Channel()]) {
                continue; 
            }

            auto assmcp = mcp_bkth_assoc->at(hit.key());
            auto assmdt = mcp_bkth_assoc->data(hit.key());
            for (unsigned int ia = 0; ia < assmcp.size(); ++ia)
            {
                auto amd = assmdt[ia];
                if (amd->isMaxIDE != 1)
                    continue;
                
                mc_hits.push_back(hit);
            }
        }
    }
}

bool SignatureIntegrityFilter::checkGeometricLength(const art::Event& e, const signature::SignatureCollection& sig_coll)
{
    auto const& mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    for (const auto& sig : sig_coll) 
    {
        const simb::MCParticle* mcp = nullptr;
        for (const auto& part : *mcp_h) 
        {
            if (part.TrackId() == sig.trckid) {
                mcp = &part;
                break;
            }
        }
        if (!mcp) 
            continue;

        TVector3 start(mcp->Vx(), mcp->Vy(), mcp->Vz());
        TVector3 end(mcp->EndX(), mcp->EndY(), mcp->EndZ());

        double length = (end - start).Mag();
        if (length < _length_thresh)
            return false;
    }

    return true; 
}

bool SignatureIntegrityFilter::checkSignatureActivelyBound(const art::Event& e, const signature::SignatureCollection& sig_coll)
{
    auto const& mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    for (const auto& sig : sig_coll) 
    {
        const simb::MCParticle* mcp = nullptr;
        for (const auto& part : *mcp_h) 
        {
            if (part.TrackId() == sig.trckid) {
                mcp = &part;
                break;
            }
        }
        if (!mcp) 
            continue;

        TVector3 start(mcp->Vx(), mcp->Vy(), mcp->Vz());
        double pos_start[3] = {start.X(), start.Y(), start.Z()};
        if (!checkPositionActivelyBound(pos_start)) {
            return false;
        }

        if (std::abs(mcp->PdgCode()) != 13) { 
            TVector3 end(mcp->EndX(), mcp->EndY(), mcp->EndZ());
            double pos_end[3] = {end.X(), end.Y(), end.Z()};
            if (!checkPositionActivelyBound(pos_end)) {
                return false;
            }
        }
    }

    return true;
}

bool SignatureIntegrityFilter::checkPositionActivelyBound(const double x[3]) 
{
    geo::TPCGeo const &this_tpc = _geo->TPC();
    geo::BoxBoundedGeo this_tpc_geo = this_tpc.ActiveBoundingBox();

    std::vector<double> bnd = {this_tpc_geo.MinX(), this_tpc_geo.MaxX(), this_tpc_geo.MinY(),
                               this_tpc_geo.MaxY(), this_tpc_geo.MinZ(), this_tpc_geo.MaxZ()};

    bool is_x = x[0] > (bnd[0] + _fid_x_start) && x[0] < (bnd[1] - _fid_x_end);
    bool is_y = x[1] > (bnd[2] + _fid_y_start) && x[1] < (bnd[3] - _fid_y_end);
    bool is_z = x[2] > (bnd[4] + _fid_z_start) && x[2] < (bnd[5] - _fid_z_end);

    return is_x && is_y && is_z;
}

bool SignatureIntegrityFilter::checkSignatureSensitivity(const art::Event& e, const signature::SignatureCollection& sig_coll) 
{
    auto const& mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    for (const auto& sig : sig_coll) 
    {
        const simb::MCParticle* mcp = nullptr;
        for (const auto& part : *mcp_h) 
        {
            if (part.TrackId() == sig.trckid) {
                mcp = &part;
                break;
            }
        }
        if (!mcp) 
            continue;

        TVector3 start(mcp->Vx(), mcp->Vy(), mcp->Vz());
        if (!checkPositionSensitivity(start)) {
            return false;
        }

        if (std::abs(mcp->PdgCode()) != 13) { 
            TVector3 end(mcp->EndX(), mcp->EndY(), mcp->EndZ());
            if (!checkPositionSensitivity(end)) {
                return false;
            }
        }
    }

    return true;
}

bool SignatureIntegrityFilter::checkPositionSensitivity(const TVector3& pos) 
{
    for (geo::PlaneID const& plane : _geo->IteratePlaneIDs()) 
    {
        try {
            geo::WireID wire = _geo->NearestWireID(pos, plane);
            raw::ChannelID_t central_channel = _geo->PlaneWireToChannel(wire);

            for (int offset = -_channel_proximity; offset <= _channel_proximity; ++offset) 
            {
                raw::ChannelID_t neighboring_channel = central_channel + offset;
                if (neighboring_channel < 0 || static_cast<size_t>(neighboring_channel) >= _geo->Nchannels())
                    continue; 

                if (_veto_bad_channels && _bad_channel_mask[neighboring_channel]) 
                    return false; 
            }
        } catch (const cet::exception&) {
            return false;
        }
    }

    return true;
}

bool SignatureIntegrityFilter::checkAngularAlignment(const art::Event& e, const signature::SignatureCollection& sig_coll) 
{
    auto const& mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    for (size_t i = 0; i < sig_coll.size(); ++i) {
        const simb::MCParticle* first_mcp = nullptr;
        for (const auto& part : *mcp_h) {
            if (part.TrackId() == sig_coll[i].trckid) {
                first_mcp = &part;
                break;
            }
        }
        if (!first_mcp) continue;

        TVector3 first_dir = TVector3(first_mcp->EndX(), first_mcp->EndY(), first_mcp->EndZ()) - 
                        TVector3(first_mcp->Vx(), first_mcp->Vy(), first_mcp->Vz());
        first_dir = first_dir.Unit();

        for (size_t j = i + 1; j < sig_coll.size(); ++j) {
            const simb::MCParticle* second_mcp = nullptr;
            for (const auto& part : *mcp_h) {
                if (part.TrackId() == sig_coll[j].trckid) {
                    second_mcp = &part;
                    break;
                }
            }
            if (!second_mcp) continue;

            TVector3 second_dir = TVector3(second_mcp->EndX(), second_mcp->EndY(), second_mcp->EndZ()) - 
                            TVector3(second_mcp->Vx(), second_mcp->Vy(), second_mcp->Vz());
            second_dir = second_dir.Unit();

            double angle = std::acos(first_dir.Dot(second_dir));
            if (angle < _angular_alignment_limit) {
                return false; 
            }
        }
    }

    return true;
}

DEFINE_ART_MODULE(SignatureIntegrityFilter)