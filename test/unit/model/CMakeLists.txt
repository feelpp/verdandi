add_executable(model_test run.cpp)
configure_file(configuration.lua configuration.lua COPYONLY)
ADD_CUSTOM_TARGET(model_test_sync ALL
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/configuration ${CMAKE_CURRENT_BINARY_DIR}/configuration
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/result ${CMAKE_CURRENT_BINARY_DIR}/result)
add_dependencies(model_test model_test_sync)
target_link_libraries(model_test verdandi gtest_main)

include(GoogleTest)
gtest_discover_tests(model_test)