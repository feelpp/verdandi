add_custom_target(clamped_bar)
ADD_CUSTOM_TARGET(clamped_bar_sync ALL
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/configuration ${CMAKE_CURRENT_BINARY_DIR}/configuration
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/result ${CMAKE_CURRENT_BINARY_DIR}/result)
add_dependencies(clamped_bar clamped_bar_sync)

add_executable(clamped_bar_forward forward.cpp)
target_link_libraries(clamped_bar_forward verdandi)
add_dependencies(clamped_bar clamped_bar_forward)