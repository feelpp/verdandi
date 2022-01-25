add_executable(performance_test run.cpp)
configure_file(configuration.lua configuration.lua COPYONLY)
ADD_CUSTOM_TARGET(performance_test_sync ALL
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/configuration ${CMAKE_CURRENT_BINARY_DIR}/configuration
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/result ${CMAKE_CURRENT_BINARY_DIR}/result
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/module_performance ${CMAKE_CURRENT_BINARY_DIR}/module_performance)
add_dependencies(performance_test performance_test_sync)
target_link_libraries(performance_test verdandi gtest_main)

include(GoogleTest)
gtest_discover_tests(performance_test)