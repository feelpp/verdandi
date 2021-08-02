
file(GLOB CPP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)
add_library(verdandi SHARED ${CPP_FILES})
target_include_directories(verdandi PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/share>
      $<INSTALL_INTERFACE:include/ops> )
target_link_libraries(verdandi PUBLIC seldon ops lua)