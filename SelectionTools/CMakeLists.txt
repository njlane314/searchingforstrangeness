cet_enable_asserts()

find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
include_directories(${TORCH_INCLUDE_DIRS})

cet_make_library(
  LIBRARY_NAME _SelectionToolBase INTERFACE
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
    lardataobj::MCBase
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
    larsim::Simulation
    lardataobj::MCBase
)

cet_build_plugin(
  FiducialVolumeSelection ubana_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    lardata::RecoBaseProxy
    larcore::Geometry_Geometry_service
    nusimdata::SimulationBase
    larpandora::LArPandoraInterface
    larevt::SpaceChargeService
    larsim::Simulation
    lardataobj::MCBase
)

cet_build_plugin(
  TopologySelection ubana_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    lardata::RecoBaseProxy
    larcore::Geometry_Geometry_service
    nusimdata::SimulationBase
    larpandora::LArPandoraInterface
    larevt::SpaceChargeService
    larsim::Simulation
    lardataobj::MCBase
    ubana::XGBoost
)

cet_build_plugin(
  ImageSelection ubana_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    lardata::RecoBaseProxy
    larcore::Geometry_Geometry_service
    nusimdata::SimulationBase
    larpandora::LArPandoraInterface
    larevt::SpaceChargeService
    larsim::Simulation
    lardataobj::MCBase
)

cet_build_plugin(
  ContrastiveSelection ubana_tool::SelectionToolBase
  LIBRARIES
  PRIVATE
    art::Framework_Core
    lardataobj::RecoBase
    lardata::RecoBaseProxy
    larcore::Geometry_Geometry_service
    nusimdata::SimulationBase
    larpandora::LArPandoraInterface
    larevt::SpaceChargeService
    larsim::Simulation
    lardataobj::MCBase
    ${TORCH_LIBRARIES}
)

set_source_files_properties(
    ContrastiveSelection_tool.cc
    PROPERTIES COMPILE_FLAGS "-Wno-error=variadic-macros"
)

install_headers()
install_source()
install_fhicl()