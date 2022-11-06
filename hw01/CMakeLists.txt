cmake_minimum_required(VERSION 3.16.3)

project(hw01)

add_subdirectory(tag)

message("TAG_LIB_INCLUDE_DIRS = ${TAG_LIB_INCLUDE_DIRS}")
message("TAG_LIB_LIBRARIES = ${TAG_LIB_LIBRARIES}")

add_executable(${PROJECT_NAME} main.c)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -coverage -lgcov" )

target_include_directories(${PROJECT_NAME} PUBLIC ${TAG_LIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${TAG_LIB_LIBRARIES})

add_subdirectory(tests)