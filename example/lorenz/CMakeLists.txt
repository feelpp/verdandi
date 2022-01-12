add_custom_target(lorenz)
ADD_CUSTOM_TARGET(lorenz_sync ALL
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/configuration ${CMAKE_CURRENT_BINARY_DIR}/configuration
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/result ${CMAKE_CURRENT_BINARY_DIR}/result)
add_dependencies(lorenz lorenz_sync)

add_executable(lorenz_forward forward.cpp)
target_link_libraries(lorenz_forward verdandi)
add_dependencies(lorenz lorenz_forward)