project (${NIS_LIB_NAME})

set(NIS_CPP_SOURCES
    endf/core.cpp
    endf/decay.cpp
)

add_library(${NIS_LIB_NAME} ${NIS_CPP_SOURCES} )

set_target_properties(${NIS_LIB_NAME}
    PROPERTIES OUTPUT_NAME ${NIS_LIB_NAME}
    VERSION ${NIS_PROJECT_VERSION}
    LIBRARY_OUTPUT_DIRECTORY ${NIS_BUILD_DIR}/lib
    ARCHIVE_OUTPUT_DIRECTORY ${NIS_BUILD_DIR}/lib
)
