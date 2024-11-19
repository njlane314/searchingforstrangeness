#include "art/Framework/Core/EDAnalyzer.h"
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

#include "TDatabasePDG.h"

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cmath>

class TrainingRegionAnalyser : public art::EDAnalyzer 
{
public:
    explicit TrainingRegionAnalyser(fhicl::ParameterSet const &pset);

    void analyze(art::Event const &e) override;
    void beginJob() override;
    void endJob() override;

private:
    art::InputTag _HitProducer, _MCPproducer, _MCTproducer, _BacktrackTag, _PFPproducer, _CLSproducer, _SHRproducer, _SLCproducer, _VTXproducer, _PCAproducer, _TRKproducer;

    int _width, _height;

    float _drift_step;
    float _wire_pitch_u, _wire_pitch_v, _wire_pitch_w;
    std::map<common::PandoraView, float> _wire_pitch;   

    calo::CalorimetryAlg* _calo_alg;
    std::vector<std::unique_ptr<::signature::SignatureToolBase>> _signatureToolsVec;

    std::map<common::PandoraView, std::array<float, 4>> _region_bounds;
    std::vector<art::Ptr<recob::Hit>> _region_hits;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> _mcp_bkth_assoc;

    TTree* _tree;
    int _hit_count_u, _hit_count_v, _hit_count_w;
    float _total_charge_u, _total_charge_v, _total_charge_w;

    void findRegionBounds(art::Event const& evt);
    void getNuVertex(art::Event const& evt, std::array<float, 3>& nu_vtx, bool& found_vertex);
    void calculateChargeCentroid(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hits, std::map<common::PandoraView, std::array<float, 2>>& q_cent_map, std::map<common::PandoraView, float>& tot_q_map);
    std::tuple<float, float, float, float> getBoundsForView(common::PandoraView view) const;

    void fillTree(const std::vector<art::Ptr<recob::Hit>>& region_hits);
};

TrainingRegionAnalyser::TrainingRegionAnalyser(fhicl::ParameterSet const &pset)
    : EDAnalyzer{pset}
    , _width{pset.get<int>("ImageWidth", 512)}
    , _height{pset.get<int>("ImageHeight", 512)}
    , _drift_step{pset.get<float>("DriftStep", 0.5)}
    , _wire_pitch_u{pset.get<float>("WirePitchU", 0.3)}
    , _wire_pitch_v{pset.get<float>("WirePitchU", 0.3)}
    , _wire_pitch_w{pset.get<float>("WirePitchU", 0.3)}
    , _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    , _PFPproducer{pset.get<art::InputTag>("PFPproducer", "pandora")}
    , _CLSproducer{pset.get<art::InputTag>("CLSproducer", "pandora")}
    , _SHRproducer{pset.get<art::InputTag>("SHRproducer", "pandora")}
    , _SLCproducer{pset.get<art::InputTag>("SLCproducer", "pandora")}
    , _VTXproducer{pset.get<art::InputTag>("VTXproducer", "pandora")}
    , _PCAproducer{pset.get<art::InputTag>("PCAproducer", "pandora")}
    , _TRKproducer{pset.get<art::InputTag>("TRKproducer", "pandora")}
{
    _calo_alg = new calo::CalorimetryAlg(pset.get<fhicl::ParameterSet>("CaloAlg"));

    _wire_pitch = {
        {common::TPC_VIEW_U, _wire_pitch_u},
        {common::TPC_VIEW_V, _wire_pitch_v},
        {common::TPC_VIEW_W, _wire_pitch_w}
    };

    const fhicl::ParameterSet &tool_psets = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (auto const &tool_pset_label : tool_psets.get_pset_names())
    {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_label);
        _signatureToolsVec.push_back(art::make_tool<::signature::SignatureToolBase>(tool_pset));
    }

    art::ServiceHandle<art::TFileService> tfs;
    _tree = tfs->make<TTree>("TrainingRegionTree", "Training Region Data");
    _tree->Branch("hit_count_u", &_hit_count_u, "hit_count_u/I");
    _tree->Branch("hit_count_v", &_hit_count_v, "hit_count_v/I");
    _tree->Branch("hit_count_w", &_hit_count_w, "hit_count_w/I");
    _tree->Branch("total_charge_u", &_total_charge_u, "total_charge_u/F");
    _tree->Branch("total_charge_v", &_total_charge_v, "total_charge_v/F");
    _tree->Branch("total_charge_w", &_total_charge_w, "total_charge_w/F");
}

void TrainingRegionAnalyser::beginJob() 
{}

void TrainingRegionAnalyser::analyze(art::Event const &evt) 
{
    _region_bounds.clear();
    _region_hits.clear(); 
    _mcp_bkth_assoc.reset();

    std::vector<signature::Signature> sig_coll;
    for (auto& signatureTool : _signatureToolsVec) {
        if (!signatureTool->identifySignalParticles(evt, sig_coll))
            return;
    }

    art::Handle<std::vector<recob::Hit>> hit_handle;
    if (evt.getByLabel(_HitProducer, hit_handle))
    {
        std::vector<art::Ptr<recob::Hit>> all_hits;
        art::fill_ptr_vector(all_hits, hit_handle);
        _mcp_bkth_assoc = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hit_handle, evt, _BacktrackTag);

        this->findRegionBounds(evt);
        if (_region_bounds.empty())
            return;

        for (const auto& hit : all_hits)
        {
            common::PandoraView view = common::GetPandoraView(hit);
            auto [drift_min, drift_max, wire_min, wire_max] = this->getBoundsForView(view);

            const auto pos = common::GetPandoraHitPosition(evt, hit, view);
            float x = pos.X();
            float z = pos.Z();

            if (x >= drift_min && x <= drift_max && z >= wire_min && z <= wire_max)
                _region_hits.push_back(hit);
        }

        this->fillTree(_region_hits);
    }
}

void TrainingRegionAnalyser::findRegionBounds(art::Event const& evt)
{
    common::ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(evt, _PFPproducer,
                                                        proxy::withAssociated<larpandoraobj::PFParticleMetadata>(_PFPproducer),
                                                        proxy::withAssociated<recob::Cluster>(_CLSproducer),
                                                        proxy::withAssociated<recob::Slice>(_SLCproducer),
                                                        proxy::withAssociated<recob::Track>(_TRKproducer),
                                                        proxy::withAssociated<recob::Vertex>(_VTXproducer),
                                                        proxy::withAssociated<recob::PCAxis>(_PCAproducer),
                                                        proxy::withAssociated<recob::Shower>(_SHRproducer),
                                                        proxy::withAssociated<recob::SpacePoint>(_PFPproducer));

    common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(evt, _CLSproducer,
                                                proxy::withAssociated<recob::Hit>(_CLSproducer));

    auto [nu_slice_hits, nu_slice] = common::getNuSliceHits(pfp_proxy, clus_proxy);
    if (nu_slice_hits.empty())
        return;

    std::map<common::PandoraView, std::array<float, 2>> q_cent_map;
    std::map<common::PandoraView, float> tot_q_map;
    common::initialiseChargeMap(q_cent_map, tot_q_map);
    this->calculateChargeCentroid(evt, nu_slice_hits, q_cent_map, tot_q_map);

    for (const auto& view : {common::TPC_VIEW_U, common::TPC_VIEW_V, common::TPC_VIEW_W}) 
    {
        const auto [x_centroid, z_centroid] = q_cent_map[view];

        float x_min = x_centroid - (_height / 2) * _drift_step;
        float x_max = x_centroid + (_height / 2) * _drift_step;
        float z_min = z_centroid - (_width / 2) * _wire_pitch[view];
        float z_max = z_centroid + (_width / 2) * _wire_pitch[view];

        _region_bounds[view] = {x_min, x_max, z_min, z_max};

        std::cout << "View: " 
            << (view == common::TPC_VIEW_U ? "U" : (view == common::TPC_VIEW_V ? "V" : "W")) 
            << ", X bounds: [" << x_min << ", " << x_max << "]"
            << ", Z bounds: [" << z_min << ", " << z_max << "]" << std::endl;
    }
}

std::tuple<float, float, float, float> TrainingRegionAnalyser::getBoundsForView(common::PandoraView view) const
{
    const auto& bounds = _region_bounds.at(view);  
    float drift_min = bounds[0]; 
    float drift_max = bounds[1]; 
    float wire_min = bounds[2];   
    float wire_max = bounds[3];   
    return std::make_tuple(drift_min, drift_max, wire_min, wire_max);
}

void TrainingRegionAnalyser::calculateChargeCentroid(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hits, std::map<common::PandoraView, std::array<float, 2>>& q_cent_map, std::map<common::PandoraView, float>& tot_q_map)
{
    for (const auto& hit : hits)
    {
        common::PandoraView view = common::GetPandoraView(hit);
        const TVector3 pos = common::GetPandoraHitPosition(evt, hit, view);
        float charge = _calo_alg->ElectronsFromADCArea(hit->Integral(), hit->WireID().Plane);

        q_cent_map[view][0] += pos.X() * charge;  
        q_cent_map[view][1] += pos.Z() * charge;  
        tot_q_map[view] += charge;
    }

    for (auto& [view, charge_center] : q_cent_map)
    {
        if (tot_q_map[view] > 0) 
        {
            charge_center[0] /= tot_q_map[view];
            charge_center[1] /= tot_q_map[view];
        }
    }
}

void TrainingRegionAnalyser::fillTree(const std::vector<art::Ptr<recob::Hit>>& region_hits) 
{
    _hit_count_u = _hit_count_v = _hit_count_w = 0;
    _total_charge_u = _total_charge_v = _total_charge_w = 0.0;

    for (const auto& hit : region_hits) 
    {
        common::PandoraView view = common::GetPandoraView(hit);
        float charge = _calo_alg->ElectronsFromADCArea(hit->Integral(), hit->WireID().Plane);

        if (view == common::TPC_VIEW_U) {
            _hit_count_u++;
            _total_charge_u += charge;
        } else if (view == common::TPC_VIEW_V) {
            _hit_count_v++;
            _total_charge_v += charge;
        } else if (view == common::TPC_VIEW_W) {
            _hit_count_w++;
            _total_charge_w += charge;
        }
    }

    _tree->Fill();
}

void TrainingRegionAnalyser::endJob() 
{}

DEFINE_ART_MODULE(TrainingRegionAnalyser)
