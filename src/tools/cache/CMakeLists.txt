project(NIS_REACTION_CACHE)

# program name and sources
set (CACHE_EXE_NAME nis_cache)

# Make the executable
link_directories (${NIS_BUILD_DIR}/lib)
add_executable (${CACHE_EXE_NAME} main.cpp)

target_link_libraries(${CACHE_EXE_NAME} PRIVATE ${NIS_LIB_NAME})

# Set the exe properties
set_target_properties(${CACHE_EXE_NAME}
                      PROPERTIES OUTPUT_NAME ${CACHE_EXE_NAME}
                      VERSION ${NIS_PROJECT_VERSION}
                      RUNTIME_OUTPUT_DIRECTORY ${NIS_BUILD_DIR}/bin)
