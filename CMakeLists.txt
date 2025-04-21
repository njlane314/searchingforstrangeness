cet_build_plugin(
  EventSelectionFilter art::EDFilter
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  ROOT::Tree
  larevt::SpaceChargeService
  larpandora::LArPandoraInterface
  ROOT::Physics
  larsim::Simulation
  ROOT::Core
  lardataobj::MCBase
)

add_subdirectory(SelectionTools)
add_subdirectory(AnalysisTools)

install_headers()
install_fhicl()
install_source()