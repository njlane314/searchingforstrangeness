#ifndef EVENTWEIGHT_ANALYSIS_CXX
#define EVENTWEIGHT_ANALYSIS_CXX

#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "AnalysisToolBase.h"
#include "Common/ProxyTypes.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace analysis {

class EventWeightAnalysis : public AnalysisToolBase {
public:
    EventWeightAnalysis(const fhicl::ParameterSet &pset);
    ~EventWeightAnalysis(){};
    void configure(fhicl::ParameterSet const & pset);
    void analyseEvent(const art::Event& event, bool is_data) override;
    void analyseSlice(const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec, bool is_data, bool is_selected) override;
    void setBranches(TTree* _tree) override;
    void resetTTree(TTree* _tree) override;

private:
    TTree *_weightstree;
    int _run;
    int _subRun;
    int _evt;

    float _weightSpline;
    float _weightTune;
    float _weightSplineTimesTune;
    float _ppfx_cv;

    double _knobRPAup;
    double _knobRPAdn;
    double _knobCCMECup;
    double _knobCCMECdn;
    double _knobAxFFCCQEup;
    double _knobAxFFCCQEdn;
    double _knobVecFFCCQEup;
    double _knobVecFFCCQEdn;
    double _knobDecayAngMECup;
    double _knobDecayAngMECdn;
    double _knobThetaDelta2Npiup;
    double _knobThetaDelta2Npidn;
    double _knobThetaDelta2NRadup;
    double _knobThetaDelta2NRaddn;
    double _knobNormCCCOHup;
    double _knobNormCCCOHdn;
    double _knobNormNCCOHup;
    double _knobNormNCCOHdn;
    double _knobxsr_scc_Fv3up;
    double _knobxsr_scc_Fv3dn;
    double _knobxsr_scc_Fa3up;
    double _knobxsr_scc_Fa3dn;
    double _RootinoFix;

    std::vector<unsigned short> _vecWeightsGenie;
    std::vector<double> _vecWeightsGenieD;
    std::vector<unsigned short> _vecWeightFlux;
    std::vector<double> _vecWeightFluxD;
    std::vector<unsigned short> _vecWeightsReint;
    std::vector<double> _vecWeightsReintD;
    std::vector<unsigned short> _vecWeightsPPFX;
    std::vector<double> _vecWeightsPPFXD;

    std::map<std::string, std::vector<double>> _mapWeight;
    std::vector<double> _vecWeightsGenie_vec;
    std::vector<int> _vecWeightsGenie_nam;
    std::vector<unsigned short> _vecWeightsGenieUp;
    std::vector<unsigned short> _vecWeightsGenieDn;

    bool _createDedicatedTree;
    bool _createMapBranch;
    bool _createFluxBranch;
    bool _createGenieBranch;
    bool _createReintBranch;
    bool _createSplineBranch;
    bool _createTuneBranch;
    bool _createSplineTimesTuneBranch;
    bool _createPPFXBranch;
    bool _SaveAllFlux;
    bool _createGenieUpDnVecs;
    int _genieAllUniverses;
    bool _makeNuMItuple;
    bool _useReweightedFlux;

    std::string _event_weight_process_name_00;
    std::string _event_weight_process_name_01;
    std::string _event_weight_process_name_02;
    std::string _event_weight_process_name_03;
    std::string _event_weight_process_name_04;
    std::string _event_weight_process_name_05;
};

EventWeightAnalysis::EventWeightAnalysis(const fhicl::ParameterSet &p) {
    _createDedicatedTree = p.get<bool>("createDedicatedTree");
    _createMapBranch = p.get<bool>("createMapBranch");
    _createFluxBranch = p.get<bool>("createFluxBranch");
    _createGenieBranch = p.get<bool>("createGenieBranch");
    _createReintBranch = p.get<bool>("createReintBranch");
    _createSplineBranch = p.get<bool>("createSplineBranch");
    _createTuneBranch = p.get<bool>("createTuneBranch");
    _createSplineTimesTuneBranch = p.get<bool>("createSplineTimesTuneBranch");
    _createPPFXBranch = p.get<bool>("createPPFXBranch", true);
    _SaveAllFlux = p.get<bool>("SaveAllFlux", true);
    _createGenieUpDnVecs = p.get<bool>("createGenieUpDnVecs", true);
    _genieAllUniverses = p.get<int>("genieAllUniverses", 500);
    _makeNuMItuple = p.get<bool>("makeNuMINtuple", true);
    _useReweightedFlux = p.get<bool>("useReweightedFlux", true);

    _event_weight_process_name_00 = p.get<std::string>("eventWeightProcessName00", "EventWeightSep24");
    _event_weight_process_name_01 = p.get<std::string>("eventWeightProcessName01", "EventWeightSep24ExtraGENIE1");
    _event_weight_process_name_02 = p.get<std::string>("eventWeightProcessName02", "EventWeightSep24ExtraGENIE2");
    _event_weight_process_name_03 = p.get<std::string>("eventWeightProcessName03", "EventWeightSep24ExtraGENIE3");
    _event_weight_process_name_04 = p.get<std::string>("eventWeightProcessName04", "EventWeightSep24ExtraGENIE4");
    _event_weight_process_name_05 = p.get<std::string>("eventWeightProcessName05", "EventWeightSep24ExtraGENIE5");

    if(_createDedicatedTree) {
        art::ServiceHandle<art::TFileService> tfs;
        _weightstree = tfs->make<TTree>("EventWeights", "EventWeights TTree");
    } else {
        _weightstree = nullptr;
    }
}

void EventWeightAnalysis::configure(fhicl::ParameterSet const & p) {}

void EventWeightAnalysis::analyseEvent(const art::Event& event, bool is_data) {
    _run = event.run();
    _subRun = event.subRun();
    _evt = event.event();

    _vecWeightsGenie  = std::vector<unsigned short>(_genieAllUniverses,1);
    _vecWeightsGenieD = std::vector<double>(_genieAllUniverses,1.0);
    _vecWeightFluxD.clear();
    _vecWeightsReintD.clear();
    _vecWeightsGenieUp.clear();
    _vecWeightsGenieDn.clear();

    if (_makeNuMItuple) {
        _vecWeightsPPFX  = std::vector<unsigned short>(600,1);
        _vecWeightsPPFXD = std::vector<double>(600,1.0);
    }

    std::vector<art::InputTag> vecTag;
    art::InputTag eventweight_tag_00("eventweightSep24","",_event_weight_process_name_00);
    art::InputTag eventweight_tag_01("eventweightSep24","",_event_weight_process_name_01);
    art::InputTag eventweight_tag_02("eventweightSep24","",_event_weight_process_name_02);
    art::InputTag eventweight_tag_03("eventweightSep24","",_event_weight_process_name_03);
    art::InputTag eventweight_tag_04("eventweightSep24","",_event_weight_process_name_04);
    art::InputTag eventweight_tag_knobs("eventweightGenieKnobs");
    
    vecTag.push_back(eventweight_tag_00);
    vecTag.push_back(eventweight_tag_01);
    vecTag.push_back(eventweight_tag_02);
    vecTag.push_back(eventweight_tag_03);
    vecTag.push_back(eventweight_tag_04);
    if (_makeNuMItuple) {
        art::InputTag eventweight_tag_05("eventweightSep24","",_event_weight_process_name_05);
        vecTag.push_back(eventweight_tag_05);
    }
    vecTag.push_back(eventweight_tag_knobs);
    
    int GenieCounter = 0;
    int PPFXCounter = 0;

    for(auto& thisTag : vecTag){
        art::Handle<std::vector<evwgh::MCEventWeight>> eventweights_handle;
        event.getByLabel(thisTag, eventweights_handle);

        if(eventweights_handle.isValid()) {
            std::vector<art::Ptr<evwgh::MCEventWeight>> eventweights;
            art::fill_ptr_vector(eventweights, eventweights_handle);

            if (eventweights.empty()) continue;

            if (_createGenieUpDnVecs && thisTag.label()=="eventweightGenieKnobs") {
                std::vector<std::string> knobList = {"AGKYpT1pi_UBGenie","AGKYxF1pi_UBGenie","AhtBY_UBGenie","AxFFCCQEshape_UBGenie","BhtBY_UBGenie",
                            "CV1uBY_UBGenie","CV2uBY_UBGenie","DecayAngMEC_UBGenie","EtaNCEL_UBGenie","FrAbs_N_UBGenie",
                            "FrAbs_pi_UBGenie","FrCEx_N_UBGenie","FrCEx_pi_UBGenie","FrInel_N_UBGenie","FrInel_pi_UBGenie",
                            "FrPiProd_N_UBGenie","FrPiProd_pi_UBGenie","FracDelta_CCMEC_UBGenie","FracPN_CCMEC_UBGenie","MFP_N_UBGenie",
                            "MFP_pi_UBGenie","MaCCQE_UBGenie","MaCCRES_UBGenie","MaNCEL_UBGenie","MaNCRES_UBGenie",
                            "MvCCRES_UBGenie","MvNCRES_UBGenie","NonRESBGvbarnCC1pi_UBGenie","NonRESBGvbarnCC2pi_UBGenie","NonRESBGvbarnNC1pi_UBGenie",
                            "NonRESBGvbarnNC2pi_UBGenie","NonRESBGvbarpCC1pi_UBGenie","NonRESBGvbarpCC2pi_UBGenie","NonRESBGvbarpNC1pi_UBGenie","NonRESBGvbarpNC2pi_UBGenie",
                            "NonRESBGvnCC1pi_UBGenie","NonRESBGvnCC2pi_UBGenie","NonRESBGvnNC1pi_UBGenie","NonRESBGvnNC2pi_UBGenie","NonRESBGvpCC1pi_UBGenie",
                            "NonRESBGvpCC2pi_UBGenie","NonRESBGvpNC1pi_UBGenie","NonRESBGvpNC2pi_UBGenie","NormCCMEC_UBGenie","NormNCMEC_UBGenie",
                            "RDecBR1eta_UBGenie","RDecBR1gamma_UBGenie","RPA_CCQE_UBGenie","Theta_Delta2Npi_UBGenie","TunedCentralValue_UBGenie",
                            "VecFFCCQEshape_UBGenie","XSecShape_CCMEC_UBGenie","splines_general_Spline"};

                std::map<std::string, std::vector<double>> evtwgt_map_knobs = eventweights.at(0)->fWeights;
                for (size_t count=0; count<knobList.size(); count++) {
                    bool knobFound = false;
                    for(std::map<std::string, std::vector<double>>::iterator it=evtwgt_map_knobs.begin(); it!=evtwgt_map_knobs.end(); ++it){
                        if (it->first != knobList[count]) continue;
                        knobFound = true;
                        std::vector<double>& vals = it->second;
                        if (vals.size()==1) {
                            float w0 = vals[0];
                            unsigned short w0short = (unsigned short)(w0*1000.);
                            if (w0 * 1000. > static_cast<double>(std::numeric_limits<unsigned short>::max()))  { w0short = std::numeric_limits<unsigned short>::max(); }
                            if (w0 < 0) { w0short = 0; }
                            _vecWeightsGenieUp.push_back(w0short);
                            _vecWeightsGenieDn.push_back((unsigned short)(1000)); 
                        }
                        else if (vals.size()==2) {
                            float w0 = vals[0];
                            unsigned short w0short = (unsigned short)(w0*1000.);
                            if (w0 * 1000. > static_cast<double>(std::numeric_limits<unsigned short>::max()))  { w0short = std::numeric_limits<unsigned short>::max(); }
                            if (w0 < 0) { w0short = 0; }
                            _vecWeightsGenieUp.push_back(w0short);
                            float w1 = vals[1];
                            unsigned short w1short = (unsigned short)(w1*1000.);
                            if (w1 * 1000. > static_cast<double>(std::numeric_limits<unsigned short>::max()))  { w1short = std::numeric_limits<unsigned short>::max(); }
                            if (w1 < 0) { w1short = 0; }
                            _vecWeightsGenieDn.push_back(w1short);
                        }
                        else {}
                    }
                    if (knobFound==false) {}
                }
            }

            std::map<std::string, std::vector<double>> evtwgt_map = eventweights.at(0)->fWeights;

            std::cout << "--- Weights from art::Event label: " << thisTag.label() << " ---" << std::endl;
            for (const auto& pair : evtwgt_map) {
                std::cout << "    Found weight name: " << pair.first << std::endl;
            }

            if (evtwgt_map.count("RPA_CCQE_UBGenie")) { 
                _knobRPAup = evtwgt_map.at("RPA_CCQE_UBGenie")[0]; 
                _knobRPAdn = evtwgt_map.at("RPA_CCQE_UBGenie")[1]; 
            }
            if (evtwgt_map.count("XSecShape_CCMEC_UBGenie")) { 
                _knobCCMECup = evtwgt_map.at("XSecShape_CCMEC_UBGenie")[0]; 
                _knobCCMECdn = evtwgt_map.at("XSecShape_CCMEC_UBGenie")[1]; 
            }
            if (evtwgt_map.count("AxFFCCQEshape_UBGenie")) { 
                _knobAxFFCCQEup = evtwgt_map.at("AxFFCCQEshape_UBGenie")[0]; 
                _knobAxFFCCQEdn = evtwgt_map.at("AxFFCCQEshape_UBGenie")[1]; 
            }
            if (evtwgt_map.count("VecFFCCQEshape_UBGenie")) { 
                _knobVecFFCCQEup = evtwgt_map.at("VecFFCCQEshape_UBGenie")[0]; 
                _knobVecFFCCQEdn = evtwgt_map.at("VecFFCCQEshape_UBGenie")[1]; 
            }
            if (evtwgt_map.count("DecayAngMEC_UBGenie")) { 
                _knobDecayAngMECup = evtwgt_map.at("DecayAngMEC_UBGenie")[0]; 
                _knobDecayAngMECdn = evtwgt_map.at("DecayAngMEC_UBGenie")[1]; 
            }
            if (evtwgt_map.count("Theta_Delta2Npi_UBGenie")) { 
                _knobThetaDelta2Npiup = evtwgt_map.at("Theta_Delta2Npi_UBGenie")[0]; 
                _knobThetaDelta2Npidn = evtwgt_map.at("Theta_Delta2Npi_UBGenie")[1]; 
            }
            if (evtwgt_map.count("ThetaDelta2NRad_UBGenie")) { 
                _knobThetaDelta2NRadup = evtwgt_map.at("ThetaDelta2NRad_UBGenie")[0]; 
                _knobThetaDelta2NRaddn = evtwgt_map.at("ThetaDelta2NRad_UBGenie")[1]; 
            }
            if (evtwgt_map.count("NormCCCOH_UBGenie")) { 
                _knobNormCCCOHup = evtwgt_map.at("NormCCCOH_UBGenie")[0]; 
                _knobNormCCCOHdn = evtwgt_map.at("NormCCCOH_UBGenie")[1]; 
            }
            if (evtwgt_map.count("NormNCCOH_UBGenie")) { 
                _knobNormNCCOHup = evtwgt_map.at("NormNCCOH_UBGenie")[0]; 
                _knobNormNCCOHdn = evtwgt_map.at("NormNCCOH_UBGenie")[1]; 
            }
            if (evtwgt_map.count("xsr_scc_Fv3_SCC")) { 
                _knobxsr_scc_Fv3up = evtwgt_map.at("xsr_scc_Fv3_SCC")[0]; 
                _knobxsr_scc_Fv3dn = evtwgt_map.at("xsr_scc_Fv3_SCC")[1]; 
            }
            if (evtwgt_map.count("xsr_scc_Fa3_SCC")) { 
                _knobxsr_scc_Fa3up = evtwgt_map.at("xsr_scc_Fa3_SCC")[0]; 
                _knobxsr_scc_Fa3dn = evtwgt_map.at("xsr_scc_Fa3_SCC")[1]; 
            }
            if (evtwgt_map.count("RootinoFix_UBGenie")) { 
                _RootinoFix = evtwgt_map.at("RootinoFix_UBGenie")[0]; 
            }

            if(evtwgt_map.count("splines_general_Spline")) _weightSpline = evtwgt_map.at("splines_general_Spline")[0];
            if(evtwgt_map.count("TunedCentralValue_UBGenie")) _weightTune = evtwgt_map.at("TunedCentralValue_UBGenie")[0];
            
            if(_weightSpline != -1.f && _weightTune != -1.f) _weightSplineTimesTune = _weightSpline * _weightTune;
            else _weightSplineTimesTune = -1.f;


            if (_makeNuMItuple){
                if (_useReweightedFlux) {
                    if(evtwgt_map.count("ppfx_oldrw_cv_UBOLDPPFXCV")) _ppfx_cv = evtwgt_map.at("ppfx_oldrw_cv_UBOLDPPFXCV")[0];
                } 
                else {
                    if(evtwgt_map.count("ppfx_cv_UBPPFXCV")) _ppfx_cv = evtwgt_map.at("ppfx_cv_UBPPFXCV")[0];
                }
            } 

            bool isFirstVectorFlux   = _vecWeightFluxD.empty();
            bool isFirstVectorReint  = _vecWeightsReintD.empty();

            for(std::map<std::string, std::vector<double>>::iterator it=evtwgt_map.begin(); it!=evtwgt_map.end(); ++it){
                std::string keyname = it->first;

                if (keyname.find("All_UBGenie") != std::string::npos || (keyname.find("All") != std::string::npos && keyname.find("UBGenie") != std::string::npos) ) {
                    for(unsigned int i = 0; i < it->second.size(); ++i) {
                        if ( (i + (100 * GenieCounter) ) < _vecWeightsGenieD.size()) 
                            _vecWeightsGenieD[i + (100 * GenieCounter) ] *= it->second[i];
                    }
                    GenieCounter += 1;
                }
                else if ( (keyname.find("ppfx_ms_UBPPFX") != std::string::npos) && _makeNuMItuple && !_useReweightedFlux ) {
                    for(unsigned int i = 0; i < it->second.size(); ++i) {
                        if ( (i + (100 * PPFXCounter) ) < _vecWeightsPPFXD.size()) 
                            _vecWeightsPPFXD[i + (100 * PPFXCounter) ] *= it->second[i];
                    }
                    PPFXCounter += 1;
                }
                else if ( (keyname.find("ppfx_oldrw_ms_UBOLDPPFX") != std::string::npos) && _makeNuMItuple && _useReweightedFlux) {
                    for(unsigned int i = 0; i < it->second.size(); ++i) {
                        if ( (i + (100 * PPFXCounter) ) < _vecWeightsPPFXD.size())
                            _vecWeightsPPFXD[i + (100 * PPFXCounter) ] *= it->second[i];
                    }
                    PPFXCounter += 1;
                } 
                else if(!_makeNuMItuple &&
                    (keyname.find("horncurrent") != std::string::npos ||
                    keyname.find("expskin") != std::string::npos ||
                    keyname.find("piplus") != std::string::npos ||
                    keyname.find("piminus") != std::string::npos ||
                    keyname.find("kplus") != std::string::npos ||
                    keyname.find("kzero") != std::string::npos ||
                    keyname.find("kminus") != std::string::npos ||
                    keyname.find("pioninexsec") != std::string::npos ||
                    keyname.find("pionqexsec") != std::string::npos ||
                    keyname.find("piontotxsec") != std::string::npos ||
                    keyname.find("nucleontotxsec") != std::string::npos ||
                    keyname.find("nucleonqexsec") != std::string::npos ||
                    keyname.find("nucleoninexsec") != std::string::npos))
                    {
                        if (_SaveAllFlux) 
                            _mapWeight.insert(*it);
                        else {
                            if(isFirstVectorFlux){
                                _vecWeightFluxD = it->second;
                                isFirstVectorFlux = false;
                            }
                            else{
                                if ( (it->second).size() == _vecWeightFluxD.size() ) {
                                    for(unsigned int i = 0; i < it->second.size(); ++i) 
                                        _vecWeightFluxD[i] *= it->second[i];
                                }
                            }
                        }
                    }
                else if ( keyname == "reinteractions_piplus_Geant4" || keyname == "reinteractions_piminus_Geant4" || keyname == "reinteractions_proton_Geant4" ) {
                    if(isFirstVectorReint){
                        _vecWeightsReintD = it->second;
                        isFirstVectorReint = false;
                    }
                    else{
                        if ( (it->second).size() == _vecWeightsReintD.size() ) {
                            for(unsigned int i = 0; i < it->second.size(); ++i)
                                _vecWeightsReintD[i] *= it->second[i];
                        }
                    }
                }
                else{ 
                    if (keyname != "splines_general_Spline" && keyname != "TunedCentralValue_UBGenie" &&
                        !(keyname.find("RPA_CCQE_UBGenie") != std::string::npos) && 
                        !(keyname.find("XSecShape_CCMEC_UBGenie") != std::string::npos) &&
                        !(keyname.find("AxFFCCQEshape_UBGenie") != std::string::npos) &&
                        !(keyname.find("VecFFCCQEshape_UBGenie") != std::string::npos) &&
                        !(keyname.find("DecayAngMEC_UBGenie") != std::string::npos) &&
                        !(keyname.find("Theta_Delta2Npi_UBGenie") != std::string::npos) &&
                        !(keyname.find("ThetaDelta2NRad_UBGenie") != std::string::npos) &&
                        !(keyname.find("NormCCCOH_UBGenie") != std::string::npos) &&
                        !(keyname.find("NormNCCOH_UBGenie") != std::string::npos) &&
                        !(keyname.find("xsr_scc_Fv3_SCC") != std::string::npos) &&
                        !(keyname.find("xsr_scc_Fa3_SCC") != std::string::npos) &&
                        !(keyname.find("RootinoFix_UBGenie") != std::string::npos) &&
                        !(keyname.find("ppfx_oldrw_cv_UBOLDPPFXCV") != std::string::npos && _makeNuMItuple && _useReweightedFlux) &&
                        !(keyname.find("ppfx_cv_UBPPFXCV") != std::string::npos && _makeNuMItuple && !_useReweightedFlux)
                        ) {
                            _mapWeight.insert(*it);
                        }
                }
            }
        }
    }
    if (!_SaveAllFlux && !_makeNuMItuple) {
        _mapWeight.insert( std::pair<std::string,std::vector<double> >("flux_all",_vecWeightFluxD) );
    }
    _mapWeight.insert( std::pair<std::string,std::vector<double> >("All_UBGenie",_vecWeightsGenieD) );
    _mapWeight.insert( std::pair<std::string,std::vector<double> >("reint_all",_vecWeightsReintD) );
    if(_makeNuMItuple){
            _mapWeight.insert( std::pair<std::string,std::vector<double> >("ppfx_all",_vecWeightsPPFXD) );
    }


    _vecWeightFlux.assign(_vecWeightFluxD.size(), 1);
    for (size_t i=0; i < _vecWeightFluxD.size(); i++) {
        auto w = _vecWeightFluxD[i];
        unsigned short wshort = (unsigned short)(w*1000.);
        if (w * 1000. > static_cast<double>(std::numeric_limits<unsigned short>::max()))  { wshort = std::numeric_limits<unsigned short>::max(); }
        if (w < 0) { wshort = 0; }
        _vecWeightFlux[i] = wshort;
    }

    _vecWeightsGenie.assign(_vecWeightsGenieD.size(), 1);
    for (size_t i=0; i < _vecWeightsGenieD.size(); i++) {
        auto w = _vecWeightsGenieD[i];
        unsigned short wshort = (unsigned short)(w*1000.);
        if (w * 1000. > static_cast<double>(std::numeric_limits<unsigned short>::max()))  { wshort = std::numeric_limits<unsigned short>::max(); }
        if (w < 0) { wshort = 0; }
        _vecWeightsGenie[i] = wshort;
    }

    _vecWeightsReint.assign(_vecWeightsReintD.size(), 1);
    for (size_t i=0; i < _vecWeightsReintD.size(); i++) {
        auto w = _vecWeightsReintD[i];
        unsigned short wshort = (unsigned short)(w*1000.);
        if (w * 1000. > static_cast<double>(std::numeric_limits<unsigned short>::max()))  { wshort = std::numeric_limits<unsigned short>::max(); }
        if (w < 0) { wshort = 0; }
        _vecWeightsReint[i] = wshort;
    }

    if (_makeNuMItuple) {
        _vecWeightsPPFX.assign(_vecWeightsPPFXD.size(), 1);
        for (size_t i=0; i < _vecWeightsPPFXD.size(); i++) {
            auto w = _vecWeightsPPFXD[i];
            unsigned short wshort = (unsigned short)(w*1000.);
            if (w * 1000. > static_cast<double>(std::numeric_limits<unsigned short>::max()))  { wshort = std::numeric_limits<unsigned short>::max(); }
            if (w < 0) { wshort = 0; }
            _vecWeightsPPFX[i] = wshort;
        }
    }

    if(_createDedicatedTree && _weightstree)
        _weightstree->Fill();
}

void EventWeightAnalysis::analyseSlice(const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec, bool is_data, bool is_selected) {}

void EventWeightAnalysis::setBranches(TTree *_tree){
    if(!_tree) return;

    if(_createDedicatedTree && _weightstree){
        _weightstree->Branch("weights", "std::map<std::string, std::vector<double>>", &_mapWeight);
        _weightstree->Branch("run",&_run,"run/I");
        _weightstree->Branch("subRun",&_subRun,"subRun/I");
        _weightstree->Branch("evt",&_evt,"evt/I");
    }
    if(_createMapBranch) _tree->Branch("weights", "std::map<std::string, std::vector<double>>", &_mapWeight);

    if(_createSplineBranch) _tree->Branch("weightSpline",&_weightSpline,"weightSpline/F");
    if(_createTuneBranch) _tree->Branch("weightTune",&_weightTune,"weightTune/F");
    if(_createSplineTimesTuneBranch) _tree->Branch("weightSplineTimesTune",&_weightSplineTimesTune,"weightSplineTimesTune/F");
    if(_createPPFXBranch) _tree->Branch("ppfx_cv",&_ppfx_cv,"ppfx_cv/F");

    if (_createGenieBranch) {
        _tree->Branch("knobRPAup",&_knobRPAup,"knobRPAup/D");
        _tree->Branch("knobRPAdn",&_knobRPAdn,"knobRPAdn/D");
        _tree->Branch("knobCCMECup",&_knobCCMECup,"knobCCMECup/D");
        _tree->Branch("knobCCMECdn",&_knobCCMECdn,"knobCCMECdn/D");
        _tree->Branch("knobAxFFCCQEup",&_knobAxFFCCQEup,"knobAxFFCCQEup/D");
        _tree->Branch("knobAxFFCCQEdn",&_knobAxFFCCQEdn,"knobAxFFCCQEdn/D");
        _tree->Branch("knobVecFFCCQEup",&_knobVecFFCCQEup,"knobVecFFCCQEup/D");
        _tree->Branch("knobVecFFCCQEdn",&_knobVecFFCCQEdn,"knobVecFFCCQEdn/D");
        _tree->Branch("knobDecayAngMECup",&_knobDecayAngMECup,"knobDecayAngMECup/D");
        _tree->Branch("knobDecayAngMECdn",&_knobDecayAngMECdn,"knobDecayAngMECdn/D");
        _tree->Branch("knobThetaDelta2Npiup",&_knobThetaDelta2Npiup,"knobThetaDelta2Npiup/D");
        _tree->Branch("knobThetaDelta2Npidn",&_knobThetaDelta2Npidn,"knobThetaDelta2Npidn/D");
        _tree->Branch("knobThetaDelta2NRadup",&_knobThetaDelta2NRadup,"knobThetaDelta2NRadup/D");
        _tree->Branch("knobThetaDelta2NRaddn",&_knobThetaDelta2NRaddn,"knobThetaDelta2NRaddn/D");
        _tree->Branch("knobNormCCCOHup",&_knobNormCCCOHup,"knobNormCCCOHup/D");
        _tree->Branch("knobNormCCCOHdn",&_knobNormCCCOHdn,"knobNormCCCOHdn/D");
        _tree->Branch("knobNormNCCOHup",&_knobNormNCCOHup,"knobNormNCCOHup/D");
        _tree->Branch("knobNormNCCOHdn",&_knobNormNCCOHdn,"knobNormNCCOHdn/D");
        _tree->Branch("knobxsr_scc_Fv3up",&_knobxsr_scc_Fv3up,"knobxsr_scc_Fv3up/D");
        _tree->Branch("knobxsr_scc_Fv3dn",&_knobxsr_scc_Fv3dn,"knobxsr_scc_Fv3dn/D");
        _tree->Branch("knobxsr_scc_Fa3up",&_knobxsr_scc_Fa3up,"knobxsr_scc_Fa3up/D");
        _tree->Branch("knobxsr_scc_Fa3dn",&_knobxsr_scc_Fa3dn,"knobxsr_scc_Fa3dn/D");
        _tree->Branch("RootinoFix",&_RootinoFix,"RootinoFix/D");
    }

    if(_createFluxBranch) _tree->Branch("weightsFlux", "std::vector<unsigned short>", &_vecWeightFlux);
    if(_createGenieBranch) _tree->Branch("weightsGenie", "std::vector<unsigned short>", &_vecWeightsGenie);
    if(_createReintBranch) _tree->Branch("weightsReint", "std::vector<unsigned short>", &_vecWeightsReint);
    if(_createPPFXBranch && _makeNuMItuple) _tree->Branch("weightsPPFX", "std::vector<unsigned short>", &_vecWeightsPPFX);

    if(_createGenieUpDnVecs) _tree->Branch("weightsGenieUp", "std::vector<unsigned short>", &_vecWeightsGenieUp);
    if(_createGenieUpDnVecs) _tree->Branch("weightsGenieDn", "std::vector<unsigned short>", &_vecWeightsGenieDn);
}

void EventWeightAnalysis::resetTTree(TTree *_tree){
    _mapWeight.clear();

    _weightSpline = -1.f;
    _weightTune = -1.f;
    _weightSplineTimesTune = -1.f;
    _ppfx_cv = -1.f;

    _knobRPAup = 1.0; _knobRPAdn = 1.0;
    _knobCCMECup = 1.0; _knobCCMECdn = 1.0;
    _knobAxFFCCQEup = 1.0; _knobAxFFCCQEdn = 1.0;
    _knobVecFFCCQEup = 1.0; _knobVecFFCCQEdn = 1.0;
    _knobDecayAngMECup = 1.0; _knobDecayAngMECdn = 1.0;
    _knobThetaDelta2Npiup = 1.0; _knobThetaDelta2Npidn = 1.0;
    _knobThetaDelta2NRadup = 1.0; _knobThetaDelta2NRaddn = 1.0;
    _knobNormCCCOHup = 1.0; _knobNormCCCOHdn = 1.0;
    _knobNormNCCOHup = 1.0; _knobNormNCCOHdn = 1.0;
    _knobxsr_scc_Fv3up = 1.0; _knobxsr_scc_Fv3dn = 1.0;
    _knobxsr_scc_Fa3up = 1.0; _knobxsr_scc_Fa3dn = 1.0;
    _RootinoFix = 1.0;

    _vecWeightFlux.clear(); _vecWeightFluxD.clear();
    _vecWeightsGenie.clear(); _vecWeightsGenieD.clear();
    _vecWeightsReint.clear(); _vecWeightsReintD.clear();
    _vecWeightsPPFX.clear(); _vecWeightsPPFXD.clear();

    _vecWeightsGenie_vec.clear();
    _vecWeightsGenie_nam.clear();
    _vecWeightsGenieUp.clear();
    _vecWeightsGenieDn.clear();

    _run = -1;
    _subRun = -1;
    _evt = -1;
}

DEFINE_ART_CLASS_TOOL(EventWeightAnalysis)

}

#endif
