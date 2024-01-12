project(NIS)

# program name and sources
set (NIS_EXE_NAME nis)
set (NIS_EXE_SOURCES nis.cpp)

# Make the executable
add_executable (${NIS_EXE_NAME} ${NIS_EXE_SOURCES})

# Set the exe properties
set_target_properties(${NIS_EXE_NAME}
                      PROPERTIES OUTPUT_NAME ${NIS_EXE_NAME}
                      VERSION ${NIS_PROJECT_VERSION}
                      RUNTIME_OUTPUT_DIRECTORY ${NIS_BUILD_DIR}/bin)
