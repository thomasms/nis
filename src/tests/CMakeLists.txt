set(CPP_UNIT_TESTS_NAME nistests)

set(CPP_UNIT_TESTS_SOURCES
  test_util.cpp
  test_zai.cpp
)

add_executable(${CPP_UNIT_TESTS_NAME} ${CPP_UNIT_TESTS_SOURCES})

target_link_libraries(${CPP_UNIT_TESTS_NAME}
  PRIVATE
    Catch2::Catch2WithMain
)

# Set the exe properties
set_target_properties(${CPP_UNIT_TESTS_NAME}
                      PROPERTIES OUTPUT_NAME ${CPP_UNIT_TESTS_NAME}
                      VERSION ${NIS_PROJECT_VERSION}
                      RUNTIME_OUTPUT_DIRECTORY ${NIS_BUILD_DIR}/bin)
