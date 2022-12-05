cmake_minimum_required(VERSION 3.16.3)

project(analyser)

if(PARALLEL)
  set(ANALYZER_PATH ${CMAKE_CURRENT_BINARY_DIR})
  configure_file(config.h.in config.h @ONLY)

  add_executable(counter counter.cpp utils.hpp utils.cpp)

  add_library(parallel_analyser SHARED analyser.hpp parallel_analyser.cpp utils.hpp utils.cpp)
  target_include_directories(
    parallel_analyser
    PRIVATE 
    ${CMAKE_CURRENT_BINARY_DIR}
  )

  find_package(Boost 1.64.0 REQUIRED COMPONENTS filesystem)
  target_link_libraries(parallel_analyser Boost::filesystem)
else()
  add_library(sequential_analyser STATIC analyser.hpp sequential_analyser.cpp utils.hpp utils.cpp)
  target_include_directories(
    sequential_analyser
    PRIVATE 
    ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()
