include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )
include_directories( $ENV{SEARCH_TOP} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(/cvmfs/larsoft.opensciencegrid.org/products/python/v2_7_14b/Linux64bit+3.10-2.17/include/python2.7)

set(PYTHON_LIB_DIR /cvmfs/larsoft.opensciencegrid.org/products/python/v2_7_14b/Linux64bit+3.10-2.17/lib)
set(PYTHON_LIBRARY ${PYTHON_LIB_DIR}/libpython2.7.so)

find_ups_product(lardata)
include_directories(${LARDATA_INC})

add_subdirectory(Common)
add_subdirectory(Products)
add_subdirectory(SelectionTools)
add_subdirectory(AnalysisTools)
add_subdirectory(Imaging)
add_subdirectory(job)

message("Checking system platform: ${CMAKE_SYSTEM_NAME}")
message("Compiler: ${CMAKE_CXX_COMPILER}")
message("Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

art_make( MODULE_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation nutools_ParticleNavigation
                           lardataobj_Simulation
                           lardataobj_RecoBase
                           lardataobj_RawData
                           lardataobj_AnalysisBase
                           lardataobj_MCBase
                           lardataalg_DetectorInfo
                           lardata_Utilities
                           lardata_RecoObjects
                           lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
                           lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                           nusimdata_SimulationBase
                           larevt_Filters
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           larevt_CalibrationDBI_Providers
                           nutools_RandomUtils_NuRandomService_service
                           ubreco_SelectionTool_OpT0FinderBase
                           ubreco_SelectionTool_OpT0FinderAlgorithms
                           ubreco_BlipRecoAlg
                           larreco_RecoAlg
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_BASIC}
                           ${ART_FRAMEWORK_SERVICES_SYSTEM}
                           art_Persistency_Common
                           art_Persistency_Provenance
                           art_Utilities
                           canvas
                           ${MF_MESSAGELOGGER}
                           ${MF_UTILITIES}
                           ${FHICLCPP}
                           ${CETLIB}
                           cetlib_except
                           ${CLHEP}
                           ${ROOT_GEOM}
                           ${ROOT_XMLIO}
                           ${ROOT_GDML}
                           ${ROOT_BASIC_LIB_LIST}
                           ${ROOT_EG}
                           ${ROOT_MINUIT}
                           ${PANDORASDK}
                           larpandora_LArPandoraInterface
                           ${PYTHON_LIBRARY}
                           larreco_Calorimetry
                           pthread
        )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")

install_headers()
install_fhicl()
install_source()

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/assets/ DESTINATION assets)
