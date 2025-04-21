cet_enable_asserts()

cet_make_library(
  LIBRARY_NAME _AnalysisToolBase INTERFACE
  SOURCE AnalysisToolBase.h
  LIBRARIES
  INTERFACE
  art::Framework_Services_Registry
  art_root_io::TFileService_service
  fhiclcpp::fhiclcpp
  ROOT::Tree
)

cet_write_plugin_builder(
  ubana_tool::AnalysisToolBase art::tool Modules
  INSTALL_BUILDER
  LIBRARIES CONDITIONAL
  ubana::AnalysisToolBase
  lardata::RecoBaseProxy
  lardataobj::RecoBase
  lardataobj::MCBase
)

include(ubana_tool::AnalysisToolBase)

cet_build_plugin(
  EventWeightAnalysis ubana_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
    larsim::EventWeight_Base
)

cet_build_plugin(
  NeutrinoAnalysis ubana_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    nusimdata::SimulationBase
    ROOT::Physics
    larevt::SpaceChargeService
    lardataalg::DetectorInfo
    lardata::Utilities
    larpandora::LArPandoraInterface
    larsim::Simulation           
    ROOT::Core      
    lardataobj::MCBase             
)

cet_build_plugin(
  NuGraphAnalysis ubana_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
    lardata::Utilities
    lardataobj::AnalysisBase
)

cet_build_plugin(
  TruthAnalysis ubana_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
    larcore::Geometry_Geometry_service
    nusimdata::SimulationBase
)

install_headers()
install_source()
install_fhicl()