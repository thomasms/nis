project(NIS_INPUT_WRITER)

# program name and sources
set (WRITER_EXE_NAME input_writer)

# Make the executable
link_directories (${NIS_BUILD_DIR}/lib)
add_executable (${WRITER_EXE_NAME} main.cpp)

target_link_libraries(${WRITER_EXE_NAME} PRIVATE ${NIS_LIB_NAME})

# Set the exe properties
set_target_properties(${WRITER_EXE_NAME}
                      PROPERTIES OUTPUT_NAME ${WRITER_EXE_NAME}
                      VERSION ${NIS_PROJECT_VERSION}
                      RUNTIME_OUTPUT_DIRECTORY ${NIS_BUILD_DIR}/bin)
