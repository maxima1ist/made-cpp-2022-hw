cmake_minimum_required(VERSION 3.16.3)

project(hw02)

add_subdirectory(tests)

add_library(matrix INTERFACE)
target_include_directories(matrix INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/matrix/")