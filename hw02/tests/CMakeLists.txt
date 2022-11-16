cmake_minimum_required(VERSION 3.16.3)

project(tests)

include(CTest)
enable_testing()
find_package(GTest REQUIRED)

add_executable(unit_test_matrix unit_test_matrix.cc)
target_link_libraries(unit_test_matrix PRIVATE matrix GTest::gtest_main)
gtest_add_tests(TARGET unit_test_matrix)

add_executable(unit_test_vector unit_test_vector.cc)
target_link_libraries(unit_test_vector PRIVATE matrix GTest::gtest_main)
gtest_add_tests(TARGET unit_test_vector)

add_executable(integration_test_matrix integration_test_matrix.cc)
target_link_libraries(integration_test_matrix PRIVATE matrix GTest::gtest_main)
gtest_add_tests(TARGET integration_test_matrix)