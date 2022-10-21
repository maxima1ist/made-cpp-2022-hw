cmake_minimum_required(VERSION 3.16.3)

project(tests)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-fprofile-arcs -ftest-coverage -O0")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -coverage -lgcov" )

file(GLOB SOURCES *.cpp)

find_package(GTest REQUIRED)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC ${TAG_LIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${TAG_LIB_LIBRARIES} GTest::gtest_main)

include (CTest)
enable_testing()
gtest_add_tests(TARGET ${PROJECT_NAME}
                TEST_SUFFIX .tag
                TEST_LIST   tagTests)

# add_compile_options(-fsanitize=address,undefined)
# add_link_options(-fsanitize=address,undefined)