cet_enable_asserts()

cet_make_library(
  LIBRARY_NAME StrangenessSelectionToolBase INTERFACE
  SOURCE SelectionToolBase.h
  LIBRARIES
  INTERFACE
    art::Framework_Core
    art::Framework_Services_Registry
    art_root_io::TFileService_service
    fhiclcpp::fhiclcpp
    ROOT::Tree
    lardataobj::RecoBase  
    lardata::RecoBaseProxy
    nusimdata::SimulationBase
)

cet_write_plugin_builder(
  ubana_tool::SelectionToolBase art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
    ubana::SelectionToolBase
)

include(ubana_tool::SelectionToolBase)

cet_build_plugin(
  EmptySelection ubana_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    lardata::RecoBaseProxy
    larcore::Geometry_Geometry_service
    nusimdata::SimulationBase
    larpandora::LArPandoraInterface
    larevt::SpaceChargeService
)

cet_build_plugin(
  TruthBackgroundSelection ubana_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    lardata::RecoBaseProxy
    larcore::Geometry_Geometry_service
    nusimdata::SimulationBase
    larpandora::LArPandoraInterface
    larevt::SpaceChargeService
)

cet_build_plugin(
  TruthSignalSelection ubana_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    lardata::RecoBaseProxy
    larcore::Geometry_Geometry_service
    nusimdata::SimulationBase
    larpandora::LArPandoraInterface
    larevt::SpaceChargeService
)

install_headers()
install_source()
install_fhicl()