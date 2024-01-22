project(NIS)

# program name and sources
set (NIS_EXE_NAME nis)
set (NIS_EXE_SOURCES nis.cpp)

# Make the executable
link_directories (${NIS_BUILD_DIR}/lib)
add_executable (${NIS_EXE_NAME} ${NIS_EXE_SOURCES})

target_link_libraries(${NIS_EXE_NAME} PRIVATE ${NIS_LIB_NAME})

# Set the exe properties
set_target_properties(${NIS_EXE_NAME}
                      PROPERTIES OUTPUT_NAME ${NIS_EXE_NAME}
                      VERSION ${NIS_PROJECT_VERSION}
                      RUNTIME_OUTPUT_DIRECTORY ${NIS_BUILD_DIR}/bin)
