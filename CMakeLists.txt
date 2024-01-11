cmake_minimum_required (VERSION 3.20)

project (NIS)

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
set( common_flags "-Wall" "-Wextra" "-Wpedantic" )
set( strict_flags "-Werror" )
set( release_flags "-O3" )
set( debug_flags "-O0" "-g" )

# set the output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

# CLI is header only
set (CLI_INCLUDES ${THIRD_PARTY_DIR}/CLI11/include)
set (ENDFTK_INCLUDES ${THIRD_PARTY_DIR}/ENDFtk/src)
set (LOG_INCLUDES ${THIRD_PARTY_DIR}/Log/src)
set (DISCO_INCLUDES ${THIRD_PARTY_DIR}/disco/src)
set (RANGE_INCLUDES ${THIRD_PARTY_DIR}/range-v3/include)
set (HU_INCLUDES ${THIRD_PARTY_DIR}/header-utilities/src)
set (SPDLOG_INCLUDES ${THIRD_PARTY_DIR}/spdlog/include)
set (EIGEN_INCLUDES ${THIRD_PARTY_DIR}/eigen/)

# set includes
include_directories(${CLI_INCLUDES}
                    ${LOG_INCLUDES}
                    ${DISCO_INCLUDES}
                    ${HU_INCLUDES}
                    ${SPDLOG_INCLUDES}
                    ${RANGE_INCLUDES}
                    ${EIGEN_INCLUDES}
                    ${ENDFTK_INCLUDES}
                    ${SOURCE_DIR_BASE})

# sources
add_subdirectory(src)
