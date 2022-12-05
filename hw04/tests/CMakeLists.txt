cmake_minimum_required(VERSION 3.16.3)

project(tests)

include(CTest)
enable_testing()
find_package(GTest REQUIRED)

add_executable(tests tests.cpp)
gtest_add_tests(TARGET tests)

target_include_directories(
  tests
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../analyser
)

if(PARALLEL)
  target_link_libraries(tests parallel_analyser GTest::Main)
else()
  target_link_libraries(tests sequential_analyser GTest::Main)
endif()
