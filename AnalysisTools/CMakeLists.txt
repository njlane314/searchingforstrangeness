cet_enable_asserts()

cet_make_library(
  LIBRARY_NAME StrangenessAnalysisToolBase INTERFACE
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
  ClassificationAnalysis ubana_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core              
    nusimdata::SimulationBase        
    larcore::Geometry_Geometry_service  
    larpandora::LArPandoraInterface
    larevt::SpaceChargeService
    larsim::Simulation           
    ROOT::Core                   
    ROOT::Physics
    lardataobj::MCBase
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
  PatternAnalysis ubana_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    nusimdata::SimulationBase
    lardataobj::RecoBase
    larpandora::LArPandoraInterface
    ROOT::Physics
    larevt::SpaceChargeService
    larsim::Simulation           
    ROOT::Core    
    lardataobj::MCBase               
)

cet_build_plugin(
  SliceAnalysis ubana_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    larcore::Geometry_Geometry_service
    nusimdata::SimulationBase
    lardataobj::RecoBase
    Eigen3::Eigen
    larevt::SpaceChargeService
    lardata::DetectorClocksService
    larpandora::LArPandoraInterface
    larsim::Simulation           
    ROOT::Core                   
    ROOT::Physics
    lardataobj::MCBase
)

cet_build_plugin(
  WireImageAnalysis ubana_tool::AnalysisToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    nusimdata::SimulationBase
    larcore::Geometry_Geometry_service
    lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
    ROOT::Tree
    larpandora::LArPandoraInterface
    larevt::SpaceChargeService
    larsim::Simulation           
    ROOT::Core                   
    ROOT::Physics
    lardataobj::MCBase
)

install_headers()
install_source()
install_fhicl()