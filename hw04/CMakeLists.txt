cmake_minimum_required(VERSION 3.16.3)

project(hw04)

set(CMAKE_CXX_FLAGS "-Wall -Werror -Wpedantic")

option(PARALLEL "PARALLEL" OFF)

add_subdirectory(analyser)
add_subdirectory(tests)