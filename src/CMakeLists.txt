# main library
add_subdirectory(nis)

# executables
if(BUILD_TOOLS)
  add_subdirectory(tools)
endif(BUILD_TOOLS)

# unit tests
if(BUILD_TESTS)
  add_subdirectory(tests)
endif(BUILD_TESTS)

# python API
if(BUILD_PY_API)
  add_subdirectory(python)
endif(BUILD_PY_API)