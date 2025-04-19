cet_build_plugin(
  StrangenessSelectionFilter art::EDFilter
  LIBRARIES
  PRIVATE
  lardata::Utilities
  lardataobj::RecoBase
  nusimdata::SimulationBase
  art_root_io::TFileService_service
  ROOT::Tree
  larevt::SpaceChargeService
  larpandora::LArPandoraInterface
)

add_subdirectory(SelectionTools)
add_subdirectory(AnalysisTools)
add_subdirectory(SignatureTools)
add_subdirectory(ClarityTools)
add_subdirectory(job)

install_headers()
install_fhicl()
install_source()