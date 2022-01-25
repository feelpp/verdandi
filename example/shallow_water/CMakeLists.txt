add_custom_target(shallow_water)
ADD_CUSTOM_TARGET(shallow_water_sync ALL
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/configuration ${CMAKE_CURRENT_BINARY_DIR}/configuration
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/result ${CMAKE_CURRENT_BINARY_DIR}/result)
add_dependencies(shallow_water shallow_water_sync)

add_executable(shallow_water_forward forward.cpp)
target_link_libraries(shallow_water_forward verdandi)
add_dependencies(shallow_water shallow_water_forward)

add_executable(shallow_water_optimal_interpolation optimal_interpolation.cpp)
target_link_libraries(shallow_water_optimal_interpolation verdandi)
add_dependencies(shallow_water shallow_water_optimal_interpolation)


add_executable(shallow_water_enkf ensemble_kalman_filter.cpp)
target_link_libraries(shallow_water_enkf verdandi)
add_dependencies(shallow_water shallow_water_enkf)