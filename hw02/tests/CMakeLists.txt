cmake_minimum_required(VERSION 3.16.3)

project(tests)

set(SOURCES unit_test_matrix.cc integration_test_matrix.cc)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} PRIVATE matrix GTest::gtest_main)

find_package(GTest REQUIRED)

include(CTest)
enable_testing()
gtest_add_tests(TARGET ${PROJECT_NAME})