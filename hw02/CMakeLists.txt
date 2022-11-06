cmake_minimum_required(VERSION 3.16.3)

project(hw01)

add_subdirectory(tests)

add_library(matrix INTERFACE)
target_include_directories(matrix SYSTEM INTERFACE matrix)