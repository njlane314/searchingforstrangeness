cet_enable_asserts()

cet_make_library(
  LIBRARY_NAME SignatureToolBase INTERFACE
  SOURCE SignatureToolBase.h
  LIBRARIES
  INTERFACE
    art::Framework_Core
    art::Framework_Services_Registry
    art_root_io::TFileService_service
    fhiclcpp::fhiclcpp
    ROOT::Tree
    nusimdata::SimulationBase
    larcore::Geometry_Geometry_service
    larevt::SpaceChargeService
)

cet_write_plugin_builder(
  ubana_tool::SignatureToolBase art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
    ubana::SignatureToolBase
)

include(ubana_tool::SignatureToolBase)

cet_build_plugin(
  KaonShortSignature ubana_tool::SignatureToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    nusimdata::SimulationBase
    larcore::Geometry_Geometry_service
    lardata::RecoBaseProxy
    larpandora::LArPandoraInterface
)

cet_build_plugin(
  LambdaSignature ubana_tool::SignatureToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    nusimdata::SimulationBase
    larcore::Geometry_Geometry_service
    lardata::RecoBaseProxy
    larpandora::LArPandoraInterface
)

cet_build_plugin(
  ChargedKaonSignature ubana_tool::SignatureToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    nusimdata::SimulationBase
    larcore::Geometry_Geometry_service
    lardata::RecoBaseProxy
    larpandora::LArPandoraInterface
)

cet_build_plugin(
  ChargedSigmaSignature ubana_tool::SignatureToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    nusimdata::SimulationBase
    larcore::Geometry_Geometry_service
    lardata::RecoBaseProxy
    larpandora::LArPandoraInterface
)

cet_build_plugin(
  MuonSignature ubana_tool::SignatureToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    nusimdata::SimulationBase
    larcore::Geometry_Geometry_service
    lardata::RecoBaseProxy
    larpandora::LArPandoraInterface
)

install_headers()
install_source()
install_fhicl()