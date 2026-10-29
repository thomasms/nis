cmake_minimum_required (VERSION 3.20)

project (NIS)

# cmake options
option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_PY_API "Build python bindings" OFF)
option(BUILD_TOOLS "Build the nis programs" ON)
option(USE_SYSTEM_CATCH2 "Use a pre-installed version of Catch2" OFF)
# Set release build by default
set (CMAKE_BUILD_TYPE Release CACHE STRING "Select which configuration to build.")

set (NIS_LIB_NAME nislib)

set (SOURCE_DIR_BASE ${CMAKE_CURRENT_LIST_DIR}/src)
set (THIRD_PARTY_DIR ${CMAKE_CURRENT_LIST_DIR}/dependencies)
set (NIS_BUILD_DIR ${CMAKE_BINARY_DIR})

# Get the NIS version, set in .VERSION file
file ( STRINGS "${CMAKE_SOURCE_DIR}/.VERSION" VERSION )
string( REPLACE "." ";" VERSION_LIST ${VERSION} )
list(GET VERSION_LIST 0 VERSION_MAJOR)
list(GET VERSION_LIST 1 VERSION_MINOR)
list(GET VERSION_LIST 2 VERSION_PATCH)
set(NIS_PROJECT_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
message ( STATUS "CMake build configuration for NIS ${NIS_PROJECT_VERSION}" )

# C++ flags
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED YES )
set( CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -Wno-unused-parameter" )# -Werror" )
set( CMAKE_CXX_FLAGS_DEBUG "-O0 -g" )
set( CMAKE_CXX_FLAGS_RELEASE "-O3" )
set( CMAKE_EXE_LINKER_FLAGS " ${CMAKE_EXE_LINKER_FLAGS}")# -static-libgcc -static-libstdc++ -static")

set(LINK_SEARCH_START_STATIC TRUE)
set(LINK_SEARCH_END_STATIC TRUE)

# set the output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

# CLI is header only
set (CLI_INCLUDES ${THIRD_PARTY_DIR}/CLI11/include)
set (ENDFTK_INCLUDES ${THIRD_PARTY_DIR}/ENDFtk/src)
set (FASTFLOAT_INCLUDES ${THIRD_PARTY_DIR}/fast_float/include)
set (RANGE_INCLUDES ${THIRD_PARTY_DIR}/range-v3/include)
set (TOOLS_INCLUDES ${THIRD_PARTY_DIR}/tools/src)
set (SPDLOG_INCLUDES ${THIRD_PARTY_DIR}/spdlog/include)
set (EIGEN_INCLUDES ${THIRD_PARTY_DIR}/eigen/)

if(BUILD_TESTS)
    # find_package(Catch2 3 REQUIRED)
    add_subdirectory(${THIRD_PARTY_DIR}/Catch2)
endif(BUILD_TESTS)

# find python
if(BUILD_PY_API)
    set(Python_VIRTUALENV FIRST)

    # give -DPython_EXECUTABLE=/path/to/python to build for that python
    # (pip and uv do this when they build the package from pyproject.toml)
    find_package(Python 3 COMPONENTS Interpreter Development.Module REQUIRED)
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 CONFIG REQUIRED)

    message(STATUS "Using PyBind")
    set(PYBIND_INC_DIR ${PYBIND_INCLUDES})

    # python output directory
    set(PY_OUTPUT_DIR ${NIS_BUILD_DIR}/bin/python3)
endif()

# set includes
include_directories(${CLI_INCLUDES}
                    ${FASTFLOAT_INCLUDES}
                    ${TOOLS_INCLUDES}
                    ${SPDLOG_INCLUDES}
                    ${RANGE_INCLUDES}
                    ${EIGEN_INCLUDES}
                    ${ENDFTK_INCLUDES}
                    ${SOURCE_DIR_BASE})

# sources
add_subdirectory(src)
