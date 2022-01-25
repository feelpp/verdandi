# add_executable(TestInstantiator TestInstantiator.cxx)
# target_link_libraries(TestInstantiator vtkCommon)
# add_test(NAME TestInstantiator
#          COMMAND TestInstantiator)

# set(VERDANDI_GTEST_MODEL "Verdandi::QuadraticModel<real>")
# set(VERDANDI_GTEST_CONFIG_PATH "\'\"configuration.lua\"\'")

add_executable(method_test run.cpp)
configure_file(configuration.lua configuration.lua COPYONLY)
ADD_CUSTOM_TARGET(method_test_sync ALL
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/configuration ${CMAKE_CURRENT_BINARY_DIR}/configuration
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/result ${CMAKE_CURRENT_BINARY_DIR}/result)
add_dependencies(method_test method_test_sync)

target_link_libraries( method_test seldon verdandi gtest_main )

include(GoogleTest)
gtest_discover_tests(method_test)