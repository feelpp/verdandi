add_custom_target(petsc_clamped_bar)
ADD_CUSTOM_TARGET(petsc_clamped_bar_sync ALL
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/configuration ${CMAKE_CURRENT_BINARY_DIR}/configuration
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/result ${CMAKE_CURRENT_BINARY_DIR}/result)
add_dependencies(petsc_clamped_bar petsc_clamped_bar_sync)

find_package(MPI REQUIRED)
if (MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
else (MPI_FOUND)
    message(SEND_ERROR "This application cannot compile without MPI")
endif (MPI_FOUND)

# PkgConfig 
find_package(PkgConfig)

# PETSc
if (PKG_CONFIG_FOUND)
    pkg_check_modules(PETSC PETSc)
endif()

if (PETSC_FOUND)
    list(APPEND COMPILE_OPTIONS ${PETSC_CFLAGS})

    include_directories(${PETSC_INCLUDE_DIRS})
    
    set(LINK_FLAGS "${LINK_FLAGS} ${PETSC_LDFLAGS}")   
    
    list(APPEND LIBRARIES ${PETSC_LINK_LIBRARIES})   
    
    set(CMAKE_REQUIRED_FLAGS ${PETSC_CFLAGS})   
    
    set(CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDE_DIRS}")
else (PETSC_FOUND)
    message(SEND_ERROR "This application cannot compile without PETSC")
endif (PETSC_FOUND)

add_executable(petsc_clamped_bar_forward forward.cpp)
target_link_libraries(petsc_clamped_bar_forward verdandi ${MPI_LIBRARIES} ${PETSC_LIBRARIES})
add_dependencies(petsc_clamped_bar petsc_clamped_bar_forward)

# add_executable(petsc_clamped_bar_optimal_interpolation optimal_interpolation.cpp)
# target_link_libraries(petsc_clamped_bar_optimal_interpolation verdandi)
# add_dependencies(petsc_clamped_bar petsc_clamped_bar_optimal_interpolation)


# add_executable(petsc_clamped_bar_enkf ensemble_kalman_filter.cpp)
# target_link_libraries(petsc_clamped_bar_enkf verdandi)
# add_dependencies(petsc_clamped_bar petsc_clamped_bar_enkf)