cmake_minimum_required(VERSION 3.14)
project(test)

# GoogleTest requires at least C++11
set(CMAKE_CXX_STANDARD 11)


include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG e2239ee6043f73722e7aa812a459f54a28552929 )

FetchContent_MakeAvailable(googletest)
# FetchContent_GetProperties(googletest)
# if(NOT googletest_POPULATED)
#     FetchContent_Populate(googletest)
# endif()
enable_testing()

add_custom_target(test)
#include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/..)

add_subdirectory(unit)
#add_subdirectory(performance)