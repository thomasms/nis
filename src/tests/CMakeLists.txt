set(CPP_UNIT_TESTS_NAME nistests)

set(CPP_UNIT_TESTS_SOURCES
  test_decayheat.cpp
  test_maths.cpp
  test_solver.cpp
  test_spectrum.cpp
  test_util.cpp
  test_zai.cpp
  endf/test_decay.cpp
  endf/test_reaction.cpp
  io/test_helpers.cpp
)

add_executable(${CPP_UNIT_TESTS_NAME} ${CPP_UNIT_TESTS_SOURCES})

target_link_libraries(${CPP_UNIT_TESTS_NAME}
  PRIVATE
    Catch2::Catch2WithMain
    ${NIS_LIB_NAME}
)

# Set the exe properties
set_target_properties(${CPP_UNIT_TESTS_NAME}
                      PROPERTIES OUTPUT_NAME ${CPP_UNIT_TESTS_NAME}
                      VERSION ${NIS_PROJECT_VERSION}
                      RUNTIME_OUTPUT_DIRECTORY ${NIS_BUILD_DIR}/bin)
