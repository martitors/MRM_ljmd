cmake_minimum_required(VERSION 3.16)

project(LJMD VERSION 1.0 LANGUAGES C CXX)

add_library(mdlib src/utilities.c src/input.c src/output.c src/force_compute.c src/cleanup.c src/verlet.c src/constants.c)

target_link_libraries(mdlib m)
target_include_directories(mdlib PUBLIC ${CMAKE_SOURCE_DIR}/include)


add_executable(main src/main.cpp)
target_link_libraries(main PUBLIC mdlib m)

include(GTest.cmake)

add_gtest_test(test_utils test_utils.cpp)
add_gtest_test(test_verlet test_verlet.cpp)

