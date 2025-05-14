cmake_minimum_required(VERSION 3.10)
project(AnalysisFramework CXX)

set(CMAKE_CXX_STANDARD 17) # Or 14 if needed
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Net Tree Graf Hist RDataFrame) # Add components you use
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIRS}) # For older CMake/ROOT, new uses target_link_libraries

# Define your source files
set(FRAMEWORK_SOURCES
    DataLoader.cpp
    EventCategories.cpp
    RDataFrameExporter.cpp
    # Add other .cpp files here
    AnalysisFramework_Dict.cxx # The dictionary file that will be generated
)

# Define your header files (for IDEs and dependency tracking)
set(FRAMEWORK_HEADERS
    DataLoader.h
    EventCategories.h
    Selections.h
    RDataFrameExporter.h
    LinkDef.h
    # Add other .h files here
)

# Generate the ROOT Dictionary
ROOT_GENERATE_DICTIONARY(AnalysisFramework_Dict ${FRAMEWORK_HEADERS} LINKDEF LinkDef.h)

# Add the library
add_library(MyAnalysisFramework SHARED ${FRAMEWORK_SOURCES})
target_link_libraries(MyAnalysisFramework PUBLIC ${ROOT_LIBRARIES})

# Optional: Install step (if you want to install the library and headers)
# install(TARGETS MyAnalysisFramework DESTINATION lib)
# install(FILES ${FRAMEWORK_HEADERS} DESTINATION include/AnalysisFramework)