#pragma once

#include <string>

namespace AnalysisFramework {
namespace Selections {

// NuMu Preselection (NUMUPRESEL from Python)
// Based on the provided Python snippet
const std::string NuMuPreselection = 
    "nslice == 1 && "
    "( (_opfilter_pe_beam > 0 && _opfilter_pe_veto < 20) || bnbdata == 1 || extdata == 1) && "
    "reco_nu_vtx_sce_x > 5 && reco_nu_vtx_sce_x < 251. && "
    "reco_nu_vtx_sce_y > -110 && reco_nu_vtx_sce_y < 110. && "
    "reco_nu_vtx_sce_z > 20 && reco_nu_vtx_sce_z < 986. && "
    "(reco_nu_vtx_sce_z < 675 || reco_nu_vtx_sce_z > 775) && " // Beam-gate dead region
    "topological_score > 0.06";

// NuMu Preselection with CRT Veto (NUMUPRESELCRT from Python)
const std::string NuMuPreselectionCRT = 
    NuMuPreselection + " && (crtveto != 1 || crthitpe < 100) && _closestNuCosmicDist > 5.";

// You can add other selection strings here as needed, for example:
// const std::string NuMuSelected = NuMuPreselectionCRT + " && n_muon_candidates > 0";

// Example for a specific strangeness signal region (conceptual)
// This would typically be applied *after* ProcessNuMuVariables and AddEventCategories
// const std::string NuMuCC_SingleKPlus_Selected = 
//     NuMuPreselectionCRT + " && "
//     "selected_muon_idx != -1 && " // A valid muon was selected
//     "event_category == 10";       // Corresponds to NuMu CC Single K+ in our new scheme

} // namespace Selections
} // namespace AnalysisFramework
