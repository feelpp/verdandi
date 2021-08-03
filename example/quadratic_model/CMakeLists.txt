add_custom_target(quadratic_model)
ADD_CUSTOM_TARGET(quadratic_model_sync ALL
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/configuration ${CMAKE_CURRENT_BINARY_DIR}/configuration
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/result ${CMAKE_CURRENT_BINARY_DIR}/result)
add_dependencies(quadratic_model quadratic_model_sync)

add_executable(quadratic_model_forward forward.cpp)
target_link_libraries(quadratic_model_forward verdandi)
add_dependencies(quadratic_model quadratic_model_forward)

add_executable(quadratic_model_optimal_interpolation optimal_interpolation.cpp)
target_link_libraries(quadratic_model_optimal_interpolation verdandi)
add_dependencies(quadratic_model quadratic_model_optimal_interpolation)

add_executable(quadratic_model_nudging nudging.cpp)
target_link_libraries(quadratic_model_nudging verdandi)
add_dependencies(quadratic_model quadratic_model_nudging)
