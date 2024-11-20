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

    int _hit_threshold;
    double _spatial_tolerance;

    std::string _bad_channel_file;
    bool _veto_bad_channels;

    calo::CalorimetryAlg* _calo_alg;
    std::vector<std::unique_ptr<::signature::SignatureToolBase>> _signatureToolsVec;

    std::vector<bool> _bad_channel_mask;
    const geo::GeometryCore* _geo;

    double _fid_x_start, _fid_y_start, _fid_z_start;
    double _fid_x_end, _fid_y_end, _fid_z_end;

    void initialiseBadChannelMask();

    void getFinalStateParticles(art::Event const& e, std::vector<simb::MCParticle>& fsp);
    void getNuVertex(art::Event const& e, TVector3& nu_vtx);
    void getSimulationHits(art::Event const& e, std::vector<art::Ptr<recob::Hit>>& mc_hits);

    bool checkSignatureActivelyBound(const art::Event& e, const signature::SignatureCollection& sig_coll);
    bool checkPositionActivelyBound(const double x[3]);

    bool checkSignatureSensitivity(const art::Event& e, const signature::SignatureCollection& sig_coll); 
    bool checkPositionSensitivity(const TVector3& pos, int tolerance = 10);
};

SignatureIntegrityFilter::SignatureIntegrityFilter(fhicl::ParameterSet const &pset)
    : EDFilter{pset}
    , _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    , _hit_threshold{pset.get<int>("HitThreshold", 1000)}
    , _spatial_tolerance{pset.get<double>("SpatialTolerance", 0.3)}
    , _bad_channel_file{pset.get<std::string>("BadChannelFile", "badchannels.txt")}
    , _veto_bad_channels{pset.get<bool>("VetoBadChannels", true)}
    , _fid_x_start{pset.get<double>("FiducialXStart", 0.0)}
    , _fid_y_start{pset.get<double>("FiducialYStart", 0.0)}
    , _fid_z_start{pset.get<double>("FiducialZStart", 0.0)}
    , _fid_x_end{pset.get<double>("FiducialXEnd", 0.0)}
    , _fid_y_end{pset.get<double>("FiducialYEnd", 0.0)}
    , _fid_z_end{pset.get<double>("FiducialZEnd", 0.0)}
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
    for (auto& signatureTool : _signatureToolsVec) {
        if (!signatureTool->identifySignalParticles(e, sig_coll))
            return false;
    }

    if (!checkSignatureActivelyBound(e, sig_coll))
        return false;

    if (!checkSignatureSensitivity(e, sig_coll)) 
        return false;

    std::vector<art::Ptr<recob::Hit>> mc_hits;
    this->getSimulationHits(e, mc_hits); 

    if (mc_hits.size() < static_cast<size_t>(_hit_threshold)) 
        return false;

    TVector3 nu_vtx;
    this->getNuVertex(e, nu_vtx);
    for (const auto& signatureTool : _signatureToolsVec) 
    {
        auto* decayTool = dynamic_cast<DecayVertexProvider*>(signatureTool.get());
        if (decayTool) 
        {
            std::optional<TVector3> decay_vtx_opt = decayTool->getDecayVertex(e);
            if (decay_vtx_opt) {
                TVector3 decay_vtx = *decay_vtx_opt;
                double sep = (nu_vtx - decay_vtx).Mag();

                if (sep < _spatial_tolerance)
                    return false;

                break;
            }
        }
    }    

    std::vector<simb::MCParticle> fsp; 
    this->getFinalStateParticles(e, fsp);
    for (const auto& mcp : fsp)
    {
        if (std::any_of(sig_coll.begin(), sig_coll.end(), [&](const auto& sig) {
            return sig.trckid == mcp.TrackId();
        })){
            continue; 
        }

        TParticlePDG* p_info = TDatabasePDG::Instance()->GetParticle(mcp.PdgCode());
        if (!p_info)
            continue;

        auto p_class = p_info->ParticleClass();
        if (std::string(p_class) == "meson")
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

bool SignatureIntegrityFilter::checkSignatureActivelyBound(const art::Event& e, const signature::SignatureCollection& sig_coll)
{
    auto const& mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    for (const auto& sig : sig_coll) 
    {
        const simb::MCParticle* mcp = nullptr;
        for (const auto& particle : *mcp_h) 
        {
            if (particle.TrackId() == sig.trckid) {
                mcp = &particle;
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
    geo::TPCGeo const &thisTPC = _geo->TPC();
    geo::BoxBoundedGeo theTpcGeo = thisTPC.ActiveBoundingBox();

    std::vector<double> bnd = {theTpcGeo.MinX(), theTpcGeo.MaxX(), theTpcGeo.MinY(),
                               theTpcGeo.MaxY(), theTpcGeo.MinZ(), theTpcGeo.MaxZ()};

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
        for (const auto& particle : *mcp_h) 
        {
            if (particle.TrackId() == sig.trckid) {
                mcp = &particle;
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

bool SignatureIntegrityFilter::checkPositionSensitivity(const TVector3& pos, int tolerance) 
{
    for (geo::PlaneID const& plane : _geo->IteratePlaneIDs()) 
    {
        try {
            geo::WireID wire = _geo->NearestWireID(pos, plane);
            raw::ChannelID_t central_channel = _geo->PlaneWireToChannel(wire);

            for (int offset = -tolerance; offset <= tolerance; ++offset) 
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

DEFINE_ART_MODULE(SignatureIntegrityFilter)