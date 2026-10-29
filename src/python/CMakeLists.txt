set (NIS_PY_NAME pynis)

set(PY_API_SOURCES
    # ${NIS_CPP_SOURCES}
    pynis.cpp
)

pybind11_add_module(${NIS_PY_NAME} MODULE ${PY_API_SOURCES})

if(APPLE)
    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

link_directories (${NIS_BUILD_DIR}/lib)
target_link_libraries(${NIS_PY_NAME} PRIVATE pybind11::module ${NIS_LIB_NAME})

set_target_properties(${NIS_PY_NAME}
                      PROPERTIES OUTPUT_NAME ${NIS_PY_NAME}
                      POSITION_INDEPENDENT_CODE TRUE
                      LIBRARY_OUTPUT_DIRECTORY ${PY_OUTPUT_DIR}
                      ARCHIVE_OUTPUT_DIRECTORY ${PY_OUTPUT_DIR}
                      RUNTIME_OUTPUT_DIRECTORY ${PY_OUTPUT_DIR}
                      INTERPROCEDURAL_OPTIMIZATION ON
                      CXX_VISIBILITY_PRESET "hidden"
                      CUDA_VISIBILITY_PRESET "hidden")

# a pip or uv build (scikit-build-core) puts the module at the top of the wheel
if(SKBUILD)
    install(TARGETS ${NIS_PY_NAME} LIBRARY DESTINATION .)
endif()
