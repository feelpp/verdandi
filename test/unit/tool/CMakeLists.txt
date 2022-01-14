enable_testing()

add_executable(
  tool_test
  run.cpp)

target_link_libraries(
    tool_test 
    verdandi gtest_main
 )

include(GoogleTest)
gtest_discover_tests(tool_test)