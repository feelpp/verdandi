cmake_minimum_required(VERSION 3.16)
project(verdandi CXX)

if(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  option(ENABLE_BUILD_WITH_TIME_TRACE "Enable -ftime-trace to generate time tracing .json files on clang" OFF)
  if(ENABLE_BUILD_WITH_TIME_TRACE)
    target_compile_options(project_options INTERFACE -ftime-trace)
  endif()
endif()

option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" OFF)

#set(BLA_VENDOR Intel10_64lp)
find_package(BLAS)
find_package(LAPACK)

include(FetchContent)

# get seldon from github
FetchContent_Declare(seldon GIT_REPOSITORY https://github.com/libseldon/seldon.git) 
FetchContent_GetProperties(seldon)
if(NOT seldon_POPULATED)
    FetchContent_Populate(seldon)
endif()
add_library(seldon SHARED ${seldon_SOURCE_DIR}/Seldon.cpp)  
target_compile_features(seldon PUBLIC cxx_std_17)
target_link_libraries(seldon ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

ADD_CUSTOM_TARGET(seldon_link_target ALL
                  COMMAND ${CMAKE_COMMAND} -E create_symlink ${seldon_SOURCE_DIR} ${seldon_SOURCE_DIR}/seldon)
add_dependencies(seldon seldon_link_target)

#set_target_properties(seldon PROPERTIES OUTPUT_NAME "seldon")
target_include_directories(seldon PUBLIC
      $<BUILD_INTERFACE:${seldon_SOURCE_DIR}/>
      $<INSTALL_INTERFACE:include/seldon> )
 
 find_library(MATH_LIBRARY_FOR_LUA m)
message(STATUS "${MATH_LIBRARY_FOR_LUA}")
find_package(Lua 5.2) 
message(STATUS "${LUA_LIBRARY_lua}")
message(STATUS "${LUA_INCLUDE_DIR}")
add_library(lua INTERFACE)  
target_include_directories(lua INTERFACE
      $<BUILD_INTERFACE:${LUA_INCLUDE_DIR}/>
      $<INSTALL_INTERFACE:${LUA_INCLUDE_DIR}> )
target_link_libraries(lua INTERFACE ${LUA_LIBRARIES})

#
# Ops
#
FetchContent_Declare(libops GIT_REPOSITORY https://github.com/feelpp/libops.git GIT_TAG be65e19ab337f52305d4a98b7361acc8e7d6db1f) 
FetchContent_GetProperties(libops)
if(NOT libops_POPULATED)
    FetchContent_Populate(libops)
endif()
file(GLOB OPS_CPP_FILES ${libops_SOURCE_DIR}/*.cxx)
add_library(ops SHARED ${OPS_CPP_FILES})  
#target_compile_definitions(ops PUBLIC -DOPS_WITH_EXCEPTION)
ADD_CUSTOM_TARGET(ops_link_target ALL
                  COMMAND ${CMAKE_COMMAND} -E create_symlink ${libops_SOURCE_DIR} ${libops_SOURCE_DIR}/ops)
add_dependencies(ops ops_link_target)
target_link_libraries(ops PUBLIC lua)
#set_target_properties(seldon PROPERTIES OUTPUT_NAME "seldon")
target_include_directories(ops PUBLIC
      $<BUILD_INTERFACE:${libops_SOURCE_DIR}/>
      $<INSTALL_INTERFACE:include/ops> )


add_subdirectory(share)

add_subdirectory(python)
add_subdirectory(example)