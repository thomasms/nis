# main library
add_subdirectory(nis)

# executables
add_subdirectory(tools)

# unit tests
if(BUILD_TESTS)
  add_subdirectory(tests)
endif(BUILD_TESTS)

# python API
if(BUILD_PY_API)
  add_subdirectory(python)
endif(BUILD_PY_API)