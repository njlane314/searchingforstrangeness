#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include "TreeUtilities.h"
#include <TTree.h>

namespace analysis {
    enum class InteractionType {
        Unknown = -1,
        QE = 0,
        Res = 1,
        DIS = 2,
        Coh = 3,
        CohElastic = 4,
        ElectronScattering = 5,
        IMDAnnihilation = 6,
        InverseBetaDecay = 7,
        GlashowResonance = 8,
        AMNuGamma = 9,
        MEC = 10,
        Diffractive = 11,
        EM = 12,
        WeakMix = 13,
        CCQE = 1001,
        NCQE = 1002,
        ResCCNuProtonPiPlus = 1003,
        ResCCNuNeutronPi0 = 1004,
        ResCCNuNeutronPiPlus = 1005,
        ResNCNuProtonPi0 = 1006,
        ResNCNuProtonPiPlus = 1007,
        ResNCNuNeutronPi0 = 1008,
        ResNCNuNeutronPiMinus = 1009,
        ResCCNuBarNeutronPiMinus = 1010,
        ResCCNuBarProtonPi0 = 1011,
        ResCCNuBarProtonPiMinus = 1012,
        ResNCNuBarProtonPi0 = 1013,
        ResNCNuBarProtonPiPlus = 1014,
        ResNCNuBarNeutronPi0 = 1015,
        ResNCNuBarNeutronPiMinus = 1016,
        ResCCNuDeltaPlusPiPlus = 1017,
        ResCCNuDelta2PlusPiMinus = 1021,
        ResCCNuBarDelta0PiMinus = 1028,
        ResCCNuBarDeltaMinusPiPlus = 1032,
        ResCCNuProtonRhoPlus = 1039,
        ResCCNuNeutronRhoPlus = 1041,
        ResCCNuBarNeutronRhoMinus = 1046,
        ResCCNuBarNeutronRho0 = 1048,
        ResCCNuSigmaPlusKaonPlus = 1053,
        ResCCNuSigmaPlusKaon0 = 1055,
        ResCCNuBarSigmaMinusKaon0 = 1060,
        ResCCNuBarSigma0Kaon0 = 1062,
        ResCCNuProtonEta = 1067,
        ResCCNuBarNeutronEta = 1070,
        ResCCNuKaonPlusLambda0 = 1073,
        ResCCNuBarKaon0Lambda0 = 1076,
        ResCCNuProtonPiPlusPiMinus = 1079,
        ResCCNuProtonPi0Pi0 = 1080,
        ResCCNuBarNeutronPiPlusPiMinus = 1085,
        ResCCNuBarNeutronPi0Pi0 = 1086,
        ResCCNuBarProtonPi0Pi0 = 1090,
        CCDIS = 1091,
        NCDIS = 1092,
        NCCOH = 1096,
        CCCOH = 1097,
        NuElectronElastic = 1098,
        InverseMuDecay = 1099,
        MEC2p2h = 1100
    };

    static const std::vector<InteractionType> allInteractionTypes = {
        InteractionType::Unknown,
        InteractionType::QE,
        InteractionType::Res,
        InteractionType::DIS,
        InteractionType::Coh,
        InteractionType::CohElastic,
        InteractionType::ElectronScattering,
        InteractionType::IMDAnnihilation,
        InteractionType::InverseBetaDecay,
        InteractionType::GlashowResonance,
        InteractionType::AMNuGamma,
        InteractionType::MEC,
        InteractionType::Diffractive,
        InteractionType::EM,
        InteractionType::WeakMix,
        InteractionType::CCQE,
        InteractionType::NCQE,
        InteractionType::ResCCNuProtonPiPlus,
        InteractionType::ResCCNuNeutronPi0,
        InteractionType::ResCCNuNeutronPiPlus,
        InteractionType::ResNCNuProtonPi0,
        InteractionType::ResNCNuProtonPiPlus,
        InteractionType::ResNCNuNeutronPi0,
        InteractionType::ResNCNuNeutronPiMinus,
        InteractionType::ResCCNuBarNeutronPiMinus,
        InteractionType::ResCCNuBarProtonPi0,
        InteractionType::ResCCNuBarProtonPiMinus,
        InteractionType::ResNCNuBarProtonPi0,
        InteractionType::ResNCNuBarProtonPiPlus,
        InteractionType::ResNCNuBarNeutronPi0,
        InteractionType::ResNCNuBarNeutronPiMinus,
        InteractionType::ResCCNuDeltaPlusPiPlus,
        InteractionType::ResCCNuDelta2PlusPiMinus,
        InteractionType::ResCCNuBarDelta0PiMinus,
        InteractionType::ResCCNuBarDeltaMinusPiPlus,
        InteractionType::ResCCNuProtonRhoPlus,
        InteractionType::ResCCNuNeutronRhoPlus,
        InteractionType::ResCCNuBarNeutronRhoMinus,
        InteractionType::ResCCNuBarNeutronRho0,
        InteractionType::ResCCNuSigmaPlusKaonPlus,
        InteractionType::ResCCNuSigmaPlusKaon0,
        InteractionType::ResCCNuBarSigmaMinusKaon0,
        InteractionType::ResCCNuBarSigma0Kaon0,
        InteractionType::ResCCNuProtonEta,
        InteractionType::ResCCNuBarNeutronEta,
        InteractionType::ResCCNuKaonPlusLambda0,
        InteractionType::ResCCNuBarKaon0Lambda0,
        InteractionType::ResCCNuProtonPiPlusPiMinus,
        InteractionType::ResCCNuProtonPi0Pi0,
        InteractionType::ResCCNuBarNeutronPiPlusPiMinus,
        InteractionType::ResCCNuBarNeutronPi0Pi0,
        InteractionType::ResCCNuBarProtonPi0Pi0,
        InteractionType::CCDIS,
        InteractionType::NCDIS,
        InteractionType::NCCOH,
        InteractionType::CCCOH,
        InteractionType::NuElectronElastic,
        InteractionType::InverseMuDecay,
        InteractionType::MEC2p2h
    };

    inline std::string GetInteractionType(InteractionType type) {
        static const std::unordered_map<InteractionType, std::string> typeNames = {
            {InteractionType::Unknown, "Unknown"},
            {InteractionType::QE, "QE"},
            {InteractionType::Res, "Res"},
            {InteractionType::DIS, "DIS"},
            {InteractionType::Coh, "Coh"},
            {InteractionType::CohElastic, "CohElastic"},
            {InteractionType::ElectronScattering, "ElectronScattering"},
            {InteractionType::IMDAnnihilation, "IMDAnnihilation"},
            {InteractionType::InverseBetaDecay, "InverseBetaDecay"},
            {InteractionType::GlashowResonance, "GlashowResonance"},
            {InteractionType::AMNuGamma, "AMNuGamma"},
            {InteractionType::MEC, "MEC"},
            {InteractionType::Diffractive, "Diffractive"},
            {InteractionType::EM, "EM"},
            {InteractionType::WeakMix, "WeakMix"},
            {InteractionType::CCQE, "CCQE"},
            {InteractionType::NCQE, "NCQE"},
            {InteractionType::ResCCNuProtonPiPlus, "ResCCNuProtonPiPlus"},
            {InteractionType::ResCCNuNeutronPi0, "ResCCNuNeutronPi0"},
            {InteractionType::ResCCNuNeutronPiPlus, "ResCCNuNeutronPiPlus"},
            {InteractionType::ResNCNuProtonPi0, "ResNCNuProtonPi0"},
            {InteractionType::ResNCNuProtonPiPlus, "ResNCNuProtonPiPlus"},
            {InteractionType::ResNCNuNeutronPi0, "ResNCNuNeutronPi0"},
            {InteractionType::ResNCNuNeutronPiMinus, "ResNCNuNeutronPiMinus"},
            {InteractionType::ResCCNuBarNeutronPiMinus, "ResCCNuBarNeutronPiMinus"},
            {InteractionType::ResCCNuBarProtonPi0, "ResCCNuBarProtonPi0"},
            {InteractionType::ResCCNuBarProtonPiMinus, "ResCCNuBarProtonPiMinus"},
            {InteractionType::ResNCNuBarProtonPi0, "ResNCNuBarProtonPi0"},
            {InteractionType::ResNCNuBarProtonPiPlus, "ResNCNuBarProtonPiPlus"},
            {InteractionType::ResNCNuBarNeutronPi0, "ResNCNuBarNeutronPi0"},
            {InteractionType::ResNCNuBarNeutronPiMinus, "ResNCNuBarNeutronPiMinus"},
            {InteractionType::ResCCNuDeltaPlusPiPlus, "ResCCNuDeltaPlusPiPlus"},
            {InteractionType::ResCCNuDelta2PlusPiMinus, "ResCCNuDelta2PlusPiMinus"},
            {InteractionType::ResCCNuBarDelta0PiMinus, "ResCCNuBarDelta0PiMinus"},
            {InteractionType::ResCCNuBarDeltaMinusPiPlus, "ResCCNuBarDeltaMinusPiPlus"},
            {InteractionType::ResCCNuProtonRhoPlus, "ResCCNuProtonRhoPlus"},
            {InteractionType::ResCCNuNeutronRhoPlus, "ResCCNuNeutronRhoPlus"},
            {InteractionType::ResCCNuBarNeutronRhoMinus, "ResCCNuBarNeutronRhoMinus"},
            {InteractionType::ResCCNuBarNeutronRho0, "ResCCNuBarNeutronRho0"},
            {InteractionType::ResCCNuSigmaPlusKaonPlus, "ResCCNuSigmaPlusKaonPlus"},
            {InteractionType::ResCCNuSigmaPlusKaon0, "ResCCNuSigmaPlusKaon0"},
            {InteractionType::ResCCNuBarSigmaMinusKaon0, "ResCCNuBarSigmaMinusKaon0"},
            {InteractionType::ResCCNuBarSigma0Kaon0, "ResCCNuBarSigma0Kaon0"},
            {InteractionType::ResCCNuProtonEta, "ResCCNuProtonEta"},
            {InteractionType::ResCCNuBarNeutronEta, "ResCCNuBarNeutronEta"},
            {InteractionType::ResCCNuKaonPlusLambda0, "ResCCNuKaonPlusLambda0"},
            {InteractionType::ResCCNuBarKaon0Lambda0, "ResCCNuBarKaon0Lambda0"},
            {InteractionType::ResCCNuProtonPiPlusPiMinus, "ResCCNuProtonPiPlusPiMinus"},
            {InteractionType::ResCCNuProtonPi0Pi0, "ResCCNuProtonPi0Pi0"},
            {InteractionType::ResCCNuBarNeutronPiPlusPiMinus, "ResCCNuBarNeutronPiPlusPiMinus"},
            {InteractionType::ResCCNuBarNeutronPi0Pi0, "ResCCNuBarNeutronPi0Pi0"},
            {InteractionType::ResCCNuBarProtonPi0Pi0, "ResCCNuBarProtonPi0Pi0"},
            {InteractionType::CCDIS, "CCDIS"},
            {InteractionType::NCDIS, "NCDIS"},
            {InteractionType::NCCOH, "NCCOH"},
            {InteractionType::CCCOH, "CCCOH"},
            {InteractionType::NuElectronElastic, "NuElectronElastic"},
            {InteractionType::InverseMuDecay, "InverseMuDecay"},
            {InteractionType::MEC2p2h, "MEC2p2h"}
        };

        auto it = typeNames.find(type);
        if (it != typeNames.end()) {
            return it->second;
        }
        return "Unknown"; 
    }

    enum class InteractionGroup {
        QE,
        Res,
        DIS,
        Coh,
        MEC,
        Other
    };

    inline std::string GetInteractionGroupName(InteractionGroup group) {
        switch (group) {
            case InteractionGroup::QE: return "QE";
            case InteractionGroup::Res: return "Res";
            case InteractionGroup::DIS: return "DIS";
            case InteractionGroup::Coh: return "Coh";
            case InteractionGroup::MEC: return "MEC";
            case InteractionGroup::Other: return "Other";
            default: return "Unknown";
        }
    }

    inline InteractionGroup GetInteractionGroup(InteractionType type) {
    switch (type) {
        case InteractionType::QE:
        case InteractionType::CCQE:
        case InteractionType::NCQE:
            return InteractionGroup::QE;
        case InteractionType::Res:
        case InteractionType::ResCCNuProtonPiPlus:
        case InteractionType::ResCCNuNeutronPi0:
        case InteractionType::ResCCNuNeutronPiPlus:
        case InteractionType::ResNCNuProtonPi0:
        case InteractionType::ResNCNuProtonPiPlus:
        case InteractionType::ResNCNuNeutronPi0:
        case InteractionType::ResNCNuNeutronPiMinus:
        case InteractionType::ResCCNuBarNeutronPiMinus:
        case InteractionType::ResCCNuBarProtonPi0:
        case InteractionType::ResCCNuBarProtonPiMinus:
        case InteractionType::ResNCNuBarProtonPi0:
        case InteractionType::ResNCNuBarProtonPiPlus:
        case InteractionType::ResNCNuBarNeutronPi0:
        case InteractionType::ResNCNuBarNeutronPiMinus:
        case InteractionType::ResCCNuDeltaPlusPiPlus:
        case InteractionType::ResCCNuDelta2PlusPiMinus:
        case InteractionType::ResCCNuBarDelta0PiMinus:
        case InteractionType::ResCCNuBarDeltaMinusPiPlus:
        case InteractionType::ResCCNuProtonRhoPlus:
        case InteractionType::ResCCNuNeutronRhoPlus:
        case InteractionType::ResCCNuBarNeutronRhoMinus:
        case InteractionType::ResCCNuBarNeutronRho0:
        case InteractionType::ResCCNuSigmaPlusKaonPlus:
        case InteractionType::ResCCNuSigmaPlusKaon0:
        case InteractionType::ResCCNuBarSigmaMinusKaon0:
        case InteractionType::ResCCNuBarSigma0Kaon0:
        case InteractionType::ResCCNuProtonEta:
        case InteractionType::ResCCNuBarNeutronEta:
        case InteractionType::ResCCNuKaonPlusLambda0:
        case InteractionType::ResCCNuBarKaon0Lambda0:
        case InteractionType::ResCCNuProtonPiPlusPiMinus:
        case InteractionType::ResCCNuProtonPi0Pi0:
        case InteractionType::ResCCNuBarNeutronPiPlusPiMinus:
        case InteractionType::ResCCNuBarNeutronPi0Pi0:
        case InteractionType::ResCCNuBarProtonPi0Pi0:
            return InteractionGroup::Res;
        case InteractionType::DIS:
        case InteractionType::CCDIS:
        case InteractionType::NCDIS:
            return InteractionGroup::DIS;
        case InteractionType::Coh:
        case InteractionType::CohElastic:
        case InteractionType::NCCOH:
        case InteractionType::CCCOH:
            return InteractionGroup::Coh;
        case InteractionType::MEC:
        case InteractionType::MEC2p2h:
            return InteractionGroup::MEC;
        default:
            return InteractionGroup::Other;
    }
}

    struct Event {
        int run, sub, evt;

        int event_type;
        tree::ManagedPointer<std::vector<int>> signature_types; 
        int is_signal;

        int nu_pdg;
        int ccnc;
        float nu_e;
        float nu_theta;
        float nu_pt;
        InteractionType interaction; 
        int target_nucleus;
        int hit_nucleon;
        float W;
        float X;
        float Y;
        float QSqr;

        tree::ManagedPointer<std::vector<std::vector<float>>> wire_input_plane_images;
        tree::ManagedPointer<std::vector<std::vector<float>>> wire_label_plane_images;
        tree::ManagedPointer<std::vector<std::vector<float>>> wire_truth_plane_images;

        float slice_total_charge;

       std::array<tree::ManagedPointer<std::vector<double>>, 3> completeness_hit;
        std::array<tree::ManagedPointer<std::vector<double>>, 3> completeness_total_hits;
        std::array<tree::ManagedPointer<std::vector<double>>, 3> exclusivity_ratio;
        std::array<tree::ManagedPointer<std::vector<double>>, 3> hit_exclusivity_fraction;
        std::array<tree::ManagedPointer<std::vector<int>>, 3> particle_start_active;
        std::array<tree::ManagedPointer<std::vector<int>>, 3> particle_end_active;

        tree::ManagedPointer<std::vector<double>> vertex_distance;

        void SetBranches(TTree* tree) {
            tree->SetBranchAddress("run", &run);
            tree->SetBranchAddress("sub", &sub);
            tree->SetBranchAddress("evt", &evt);
            tree->SetBranchAddress("event_type", &event_type);
            tree::set_object_input_branch_address(*tree, "signature_types", signature_types);
            tree->SetBranchAddress("is_signal", &is_signal);
            tree->SetBranchAddress("nu_pdg", &nu_pdg);
            tree->SetBranchAddress("ccnc", &ccnc);
            tree->SetBranchAddress("nu_e", &nu_e);
            tree->SetBranchAddress("nu_theta", &nu_theta);
            tree->SetBranchAddress("nu_pt", &nu_pt);
            tree->SetBranchAddress("interaction", &interaction); 
            tree->SetBranchAddress("target_nucleus", &target_nucleus);
            tree->SetBranchAddress("hit_nucleon", &hit_nucleon);
            tree->SetBranchAddress("W", &W);
            tree->SetBranchAddress("X", &X);
            tree->SetBranchAddress("Y", &Y);
            tree->SetBranchAddress("QSqr", &QSqr);

            tree::set_object_input_branch_address(*tree, "slice_wire_images", wire_input_plane_images);
            tree::set_object_input_branch_address(*tree, "slice_label_wire_images", wire_label_plane_images);
            tree::set_object_input_branch_address(*tree, "slice_truth_wire_images", wire_truth_plane_images);

            tree->SetBranchAddress("slice_total_charge", &slice_total_charge);

            tree::set_object_input_branch_address(*tree, "completeness_hit_U", completeness_hit[0]);
            tree::set_object_input_branch_address(*tree, "completeness_hit_V", completeness_hit[1]);
            tree::set_object_input_branch_address(*tree, "completeness_hit_W", completeness_hit[2]);
            tree::set_object_input_branch_address(*tree, "completeness_total_hits_U", completeness_total_hits[0]);
            tree::set_object_input_branch_address(*tree, "completeness_total_hits_V", completeness_total_hits[1]);
            tree::set_object_input_branch_address(*tree, "completeness_total_hits_W", completeness_total_hits[2]);
            tree::set_object_input_branch_address(*tree, "exclusivity_ratio_U", exclusivity_ratio[0]);
            tree::set_object_input_branch_address(*tree, "exclusivity_ratio_V", exclusivity_ratio[1]);
            tree::set_object_input_branch_address(*tree, "exclusivity_ratio_W", exclusivity_ratio[2]);
            tree::set_object_input_branch_address(*tree, "hit_exclusivity_fraction_U", hit_exclusivity_fraction[0]);
            tree::set_object_input_branch_address(*tree, "hit_exclusivity_fraction_V", hit_exclusivity_fraction[1]);
            tree::set_object_input_branch_address(*tree, "hit_exclusivity_fraction_W", hit_exclusivity_fraction[2]);
            tree::set_object_input_branch_address(*tree, "particle_start_active_U", particle_start_active[0]);
            tree::set_object_input_branch_address(*tree, "particle_start_active_V", particle_start_active[1]);
            tree::set_object_input_branch_address(*tree, "particle_start_active_W", particle_start_active[2]);
            tree::set_object_input_branch_address(*tree, "particle_end_active_U", particle_end_active[0]);
            tree::set_object_input_branch_address(*tree, "particle_end_active_V", particle_end_active[1]);
            tree::set_object_input_branch_address(*tree, "particle_end_active_W", particle_end_active[2]);
            tree::set_object_input_branch_address(*tree, "vertex_distance", vertex_distance);
        }
    };
}