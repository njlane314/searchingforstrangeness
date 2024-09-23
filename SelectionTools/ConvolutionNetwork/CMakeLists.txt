cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(main)

set(CMAKE_CXX_STANDARD 17)

find_package(Torch REQUIRED)
find_package(TBB REQUIRED)
find_package(Protobuf REQUIRED)
find_package(ROOT REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

add_executable(main main.cpp)
target_link_libraries(main ${TORCH_LIBRARIES} ${TBB_LIBRARIES} ${ROOT_LIBRARIES} protobuf::libprotobuf pthread)

set_target_properties(main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin
    OUTPUT_NAME "main"
    SUFFIX "" 
)

set_property(TARGET main PROPERTY CXX_STANDARD 17)