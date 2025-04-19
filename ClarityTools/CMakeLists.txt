cet_enable_asserts()

cet_make_library(
  LIBRARY_NAME ClarityToolBase INTERFACE
  SOURCE ClarityToolBase.h
  LIBRARIES
  INTERFACE
    art::Framework_Core
    art::Framework_Principal
    art::Persistency_Common
    art::Utilities
    canvas::canvas
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    larcore::Geometry_Geometry_service
    lardataobj::RecoBase
    nusimdata::SimulationBase
    ROOT::Core
    ROOT::Physics
    art_root_io::TFileService_service
)

cet_write_plugin_builder(
  ubana_tool::ClarityToolBase art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
    ubana::ClarityToolBase
)

include(ubana_tool::ClarityToolBase)

cet_build_plugin(
  HitExclusivity ubana_tool::ClarityToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    nusimdata::SimulationBase
    larcore::Geometry_Geometry_service
    larpandora::LArPandoraInterface
    art_root_io::TFileService_service
    larevt::SpaceChargeService
    lardataobj::MCBase
)

cet_build_plugin(
  PatternCompleteness ubana_tool::ClarityToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    nusimdata::SimulationBase
    larcore::Geometry_Geometry_service
    larpandora::LArPandoraInterface
    art_root_io::TFileService_service
    larevt::SpaceChargeService
    lardataobj::MCBase
)

cet_build_plugin(
  SignatureIntegrity ubana_tool::ClarityToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    nusimdata::SimulationBase
    larcore::Geometry_Geometry_service
    larpandora::LArPandoraInterface
    art_root_io::TFileService_service
    larevt::SpaceChargeService
    larevt::SpaceChargeService
    lardataobj::MCBase
)

install_headers()
install_source()
install_fhicl()