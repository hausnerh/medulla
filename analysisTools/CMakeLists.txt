cmake_minimum_required(VERSION 3.14)
project(analysisTools LANGUAGES CXX)

# Set optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
#set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Tree)
include(${ROOT_USE_FILE})

# Add library for analysis tools
add_library(anatools SHARED
    src/cut_sequence.cc
    src/plot_structs.cc
    src/analysis_tree.cc
)

# Set include directories
target_include_directories(anatools PRIVATE
    ./
)

# Set link directories
target_link_directories(anatools PRIVATE
)

# Link required libraries
target_link_libraries(anatools INTERFACE
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Tree
    ROOT::Graf
    ROOT::Gui
)

# Set compile features
target_compile_features(anatools PRIVATE cxx_std_17)

# Add analysis executable(s)
add_executable(run_analysis_gOre src/run_analysis_gOre.cc)
target_link_libraries(run_analysis_gOre PRIVATE anatools)
target_include_directories(run_analysis_gOre PRIVATE . include/)

# Optionally build documentation
option(BUILD_DOCS "Build documentation" OFF)

if(BUILD_DOCS)
    find_package(Doxygen QUIET)
    find_package(Sphinx QUIET)
    if(DOXYGEN_FOUND)
        message(STATUS "Doxygen found: enabling C++ API documentation.")
    else()
        message(WARNING "Doxygen not found: skipping C++ API documentation.")
    endif()

    if(SPHINX_FOUND)
        message(STATUS "Sphinx found: enabling Python + unified documentation.")
    else()
        message(WARNING "Sphinx not found: skipping Sphinx-based documentation.")
    endif()

    if(DOXYGEN_FOUND AND SPHINX_FOUND)
        add_subdirectory(docs)
    endif()
endif()