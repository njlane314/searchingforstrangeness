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
#include "CommonFunctions/Visualisation.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "SignatureTools/SignatureToolBase.h"
#include "SignatureTools/VertexToolBase.h"

#include "ClarityTools/ClarityToolBase.h"

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

class PatternClarityFilter : public art::EDFilter 
{
public:
    explicit PatternClarityFilter(fhicl::ParameterSet const &pset);

    PatternClarityFilter(PatternClarityFilter const &) = delete;
    PatternClarityFilter(PatternClarityFilter &&) = delete;
    PatternClarityFilter &operator=(PatternClarityFilter const &) = delete;
    PatternClarityFilter &operator=(PatternClarityFilter &&) = delete;

    bool filter(art::Event &e) override;

private:
    art::InputTag _HitProducer, _MCPproducer, _MCTproducer, _BacktrackTag;

    const geo::GeometryCore* _geo;

    std::string _bad_channel_file;
    std::vector<bool> _bad_channel_mask;

    double _patt_hit_comp_thresh;
    int _patt_hit_thresh;
    double _sig_hit_comp_thresh; 
    int _chan_act_reg;
    double _hit_exclus_thresh;
    double _sig_exclus_thresh;

    calo::CalorimetryAlg* _calo_alg;
    std::vector<std::unique_ptr<::signature::SignatureToolBase>> _signatureToolsVec;
    std::vector<std::unique_ptr<::claritytools::ClarityToolBase>> _clarityToolsVec;
    int _targetDetectorPlane;
    bool _quickVisualise;

    bool filterPatternCompleteness(art::Event &e, signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc);
    bool filterSignatureIntegrity(art::Event &e, signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc);
    bool filterHitExclusivity(art::Event &e, signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc); 
};

PatternClarityFilter::PatternClarityFilter(fhicl::ParameterSet const &pset)
    : EDFilter{pset}
    , _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    //, _bad_channel_file{pset.get<std::string>("BadChannelFile", "badchannels.txt")}
    , _patt_hit_comp_thresh{pset.get<double>("PatternHitCompletenessThreshold", 0.5)}
    , _patt_hit_thresh{pset.get<int>("PatternHitThreshold", 100)}
    , _sig_hit_comp_thresh{pset.get<double>("SignatureHitCompletenessThreshold", 0.1)}
    , _chan_act_reg{pset.get<int>("ChannelActiveRegion", 3)}
    , _hit_exclus_thresh{pset.get<double>("HitExclusivityThreshold", 0.5)}
    , _sig_exclus_thresh{pset.get<double>("SignatureExclusivityThreshold", 0.8)}
    , _targetDetectorPlane{pset.get<int>("TargetDetectorPlane", 2)}
    , _quickVisualise{pset.get<bool>("QuickVisualise", true)}
{
    _calo_alg = new calo::CalorimetryAlg(pset.get<fhicl::ParameterSet>("CaloAlg"));

    const fhicl::ParameterSet &tool_psets = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (auto const &tool_pset_label : tool_psets.get_pset_names())
    {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_label);
        _signatureToolsVec.push_back(art::make_tool<::signature::SignatureToolBase>(tool_pset));
    };

    const fhicl::ParameterSet &claritytool_psets = pset.get<fhicl::ParameterSet>("ClarityTools");
    for (auto const &tool_pset_label : claritytool_psets.get_pset_names())
    {
      auto const tool_pset = claritytool_psets.get<fhicl::ParameterSet>(tool_pset_label);
      _clarityToolsVec.push_back(art::make_tool<::claritytools::ClarityToolBase>(tool_pset));
    };

    /*
    _geo = art::ServiceHandle<geo::Geometry>()->provider();
    size_t num_channels = _geo->Nchannels();
    _bad_channel_mask.resize(num_channels, false);

    if (!_bad_channel_file.empty()) {
        cet::search_path sp("FW_SEARCH_PATH");
        std::string fullname;
        sp.find_file(_bad_channel_file, fullname);
        if (fullname.empty()) 
            throw cet::exception("PatternClarityFilter") << "-- Bad channel file not found: " << _bad_channel_file;

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
    }
   */
}

bool PatternClarityFilter::filter(art::Event &e) 
{
    signature::Pattern patt;
    for (auto &signatureTool : _signatureToolsVec) {
        signature::Signature signature;
        if (!signatureTool->constructSignature(e, signature))
            return false;

        patt.push_back(signature);
    }
/*
    art::Handle<std::vector<recob::Hit>> hit_h;
    if (!e.getByLabel(_HitProducer, hit_h)) 
        return false;

    std::vector<art::Ptr<recob::Hit>> evt_hits;
    art::fill_ptr_vector(evt_hits, hit_h);
    auto mcp_bkth_assoc = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hit_h, e, _BacktrackTag);

    std::vector<art::Ptr<recob::Hit>> mc_hits;
    for (const auto& hit : evt_hits) {
        if (_bad_channel_mask[hit->Channel()]) 
            continue; 

        const geo::WireID& wire_id = hit->WireID(); 
        if (wire_id.Plane != static_cast<unsigned int>(_targetDetectorPlane))
            continue;

        auto assmcp = mcp_bkth_assoc->at(hit.key());
        auto assmdt = mcp_bkth_assoc->data(hit.key());
        for (unsigned int ia = 0; ia < assmcp.size(); ++ia){
            auto amd = assmdt[ia];
            if (amd->isMaxIDEN != 1)
                continue;
            
            mc_hits.push_back(hit);
        }
    }
*/

   for (auto &clarityTool : _clarityToolsVec)
       if(!clarityTool->filter(e, patt)) return false;

    // A clear pattern is defined as requiring that:
    // 1) the interaction topology is dominated by its specific pattern, 
    // 2) that each signature of the pattern retains its integrity within the detector, 
    // 3) and that most of the hits of the signature are exclusive. 
/*
    std::cout << "Testing HitExclusivity" << std::endl;
    if(this->filterHitExclusivity(e, patt, mc_hits, mcp_bkth_assoc)) std::cout << "Pass" << std::endl;
    else std::cout << "Fail" << std::endl;

    std::cout << "Testing PatternCompleteness" << std::endl;
    if(this->filterPatternCompleteness(e, patt, mc_hits, mcp_bkth_assoc)) std::cout << "Pass" << std::endl;
    else std::cout << "Fail" << std::endl;

    std::cout << "Testing SignatureIntegrity" << std::endl;
    if(this->filterSignatureIntegrity(e, patt, mc_hits, mcp_bkth_assoc)) std::cout << "Pass" << std::endl;
    else std::cout << "Fail" << std::endl;
*/

/*     
    if (!this->filterPatternCompleteness(e, patt, mc_hits, mcp_bkth_assoc)){
        return false;
    }

    if (!this->filterSignatureIntegrity(e, patt, mc_hits, mcp_bkth_assoc))
        return false;

    if (!this->filterHitExclusivity(e, patt, mc_hits, mcp_bkth_assoc))
        return false;
*/

    if (_quickVisualise)
    {
        std::string filename = "event_" + std::to_string(e.run()) + "_" + std::to_string(e.subRun()) + "_" + std::to_string(e.event());
        common::visualiseTrueEvent(e, _MCPproducer, _HitProducer, _BacktrackTag, filename);
        common::visualiseSignature(e, _MCPproducer, _HitProducer, _BacktrackTag, patt, filename);
    }

    return true; 
}
/*
bool PatternClarityFilter::filterPatternCompleteness(art::Event &e, signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc)
{
    std::unordered_map<int, int> sig_hit_map;
    double tot_patt_hit = 0; 

    std::vector<art::Ptr<recob::Hit>> patt_hits;
    for (const auto& sig : patt) {
        for (const auto& mcp_s : sig) {
            double sig_hit = 0;

            for (const auto& hit : mc_hits) {
                auto assmcp = mcp_bkth_assoc->at(hit.key());
                auto assmdt = mcp_bkth_assoc->data(hit.key());

                for (unsigned int ia = 0; ia < assmcp.size(); ++ia){
                    auto amd = assmdt[ia];
                    if (assmcp[ia]->TrackId() == mcp_s->TrackId() && amd->isMaxIDEN == 1) {
                        patt_hits.push_back(hit);
                        sig_hit += 1; 
                    }
                }
            }

            sig_hit_map[mcp_s->TrackId()] += sig_hit;
            tot_patt_hit += sig_hit;
        }
    }

    if (mc_hits.empty() || patt_hits.empty()) 
        return false;

    double patt_comp = static_cast<double>(patt_hits.size()) / mc_hits.size();
    //std::cout << "Pattern completeness " << patt_comp << std::endl;
    //std::cout << "Total pattern hits " << tot_patt_hit << std::endl;
    if (patt_comp < _patt_hit_comp_thresh || tot_patt_hit < _patt_hit_thresh)
        return false;

    for (const auto& [_, num_hits] : sig_hit_map) 
    {
        //std::cout << "Signature hit " << num_hits << std::endl;
        if (num_hits / tot_patt_hit < _sig_hit_comp_thresh) 
            return false;       
    }

    return true;
}

bool PatternClarityFilter::filterSignatureIntegrity(art::Event &e, signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc)
{
    auto isChannelRegionActive = [&](const TVector3& point) -> bool {
        for (geo::PlaneID const& plane : _geo->IteratePlaneIDs()) {
            try {
                geo::WireID wire = _geo->NearestWireID(point, plane);
                raw::ChannelID_t central_channel = _geo->PlaneWireToChannel(wire);

                for (int offset = -_chan_act_reg; offset <= _chan_act_reg; ++offset) {
                    raw::ChannelID_t neighboring_channel = central_channel + offset;

                    if (neighboring_channel < 0 || static_cast<size_t>(neighboring_channel) >= _geo->Nchannels())
                        continue; 

                    if (_bad_channel_mask[neighboring_channel])
                        return false; 
                }
            } catch (const cet::exception&) {
                return false; 
            }
        }
        return true;
    };

    for (const auto& sig : patt) {
        for (const auto& mcp_s : sig) {
            TVector3 start(mcp_s->Vx(), mcp_s->Vy(), mcp_s->Vz());
            if (!isChannelRegionActive(start))
                return false;

            if (std::abs(mcp_s->PdgCode()) != 13) { 
                TVector3 end(mcp_s->EndX(), mcp_s->EndY(), mcp_s->EndZ());
                if (!isChannelRegionActive(end))
                    return false;
            }
        }
    }

    return true;
}

bool PatternClarityFilter::filterHitExclusivity(art::Event &e, signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc)
{
    for (const auto& sig : patt) {
        double sig_q_inclusive = 0.0;
        double sig_q_exclusive = 0.0;
        for (const auto& mcp_s : sig) {
            for (const auto& hit : mc_hits) {
                auto assmcp = mcp_bkth_assoc->at(hit.key());
                auto assmdt = mcp_bkth_assoc->data(hit.key());

                for (unsigned int ia = 0; ia < assmcp.size(); ++ia){
                    auto amd = assmdt[ia];
                    if (assmcp[ia]->TrackId() == mcp_s->TrackId()) {
                        sig_q_inclusive += amd->numElectrons * amd->ideNFraction;
                        if (amd->ideNFraction > _hit_exclus_thresh) 
                            sig_q_exclusive += amd->numElectrons * amd->ideNFraction;
                    }
                }
            }
        }

        if (sig_q_exclusive / sig_q_inclusive < _sig_exclus_thresh)
            return false;
    }

    return true;
}
*/
DEFINE_ART_MODULE(PatternClarityFilter)
