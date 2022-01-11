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

add_executable(quadratic_model_enkf ensemble_kalman_filter.cpp)
target_link_libraries(quadratic_model_enkf verdandi)
add_dependencies(quadratic_model quadratic_model_enkf)

add_executable(quadratic_model_extkf extended_kalman_filter.cpp)
target_link_libraries(quadratic_model_extkf verdandi)
add_dependencies(quadratic_model quadratic_model_extkf)

add_executable(quadratic_model_extminmax extended_minimax_filter.cpp)
target_link_libraries(quadratic_model_extminmax verdandi)
add_dependencies(quadratic_model quadratic_model_extminmax)

add_executable(quadratic_model_mc monte_carlo.cpp)
target_link_libraries(quadratic_model_mc verdandi)
add_dependencies(quadratic_model quadratic_model_mc)

add_executable(quadratic_model_ukf unscented_kalman_filter.cpp)
target_link_libraries(quadratic_model_ukf verdandi)
add_dependencies(quadratic_model quadratic_model_ukf)

# add_executable(quadratic_model_hjb hjb.cpp)
# target_link_libraries(quadratic_model_hjb verdandi)
# add_dependencies(quadratic_model quadratic_model_hjb)