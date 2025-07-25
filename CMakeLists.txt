include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )
include_directories( $ENV{SEARCH_TOP} )

find_package(Torch REQUIRED)
include_directories($ENV{LIBTORCH_FQ_DIR}/lib/python2.7/site-packages/torch/lib/include/torch/csrc/api/include)
include_directories($ENV{LIBTORCH_FQ_DIR}/lib/python2.7/site-packages/torch/lib/include)
message("TORCH_LIBRARIES: ${TORCH_LIBRARIES}")

set(TORCH_LIBRARIES $ENV{LIBTORCH_FQ_DIR}/lib/python2.7/site-packages/torch/lib/libtorch.so)
set(LIBTORCH_LIB_DIR /cvmfs/uboone.opensciencegrid.org/products/libtorch/v1_0_1/Linux64bit+3.10-2.17-e17-prof/lib/python2.7/site-packages/torch/lib)
file(GLOB LIBTORCH_LIBRARIES "${LIBTORCH_LIB_DIR}/*.so")

include_directories(/cvmfs/larsoft.opensciencegrid.org/products/python/v2_7_14b/Linux64bit+3.10-2.17/include/python2.7)

set(PYTHON_LIB_DIR /cvmfs/larsoft.opensciencegrid.org/products/python/v2_7_14b/Linux64bit+3.10-2.17/lib)
set(PYTHON_LIBRARY ${PYTHON_LIB_DIR}/libpython2.7.so)

add_subdirectory(CommonDefs)
add_subdirectory(SelectionTools)
add_subdirectory(AnalysisTools)

message("Checking system platform: ${CMAKE_SYSTEM_NAME}")
message("Compiler: ${CMAKE_CXX_COMPILER}")
message("Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

art_make( MODULE_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           larsim_Simulation nutools_ParticleNavigation
                           lardataobj_Simulation
                           lardataobj_RecoBase
                           lardataobj_MCBase	
                           lardata_Utilities
			   nusimdata_SimulationBase
                           larevt_Filters
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           lardataobj_RawData
                           larevt_CalibrationDBI_Providers
                           nutools_RandomUtils_NuRandomService_service
                           ${ART_FRAMEWORK_CORE}
                           ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_SERVICES_REGISTRY}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
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
                           ${TORCH_LIBRARIES}
                           ${LIBTORCH_LIBRARIES}
                           ${PYTHON_LIBRARY}
                           larreco_Calorimetry
                           pthread
        )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")

install_headers()
install_fhicl()
install_source()