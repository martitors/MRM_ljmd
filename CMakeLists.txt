cmake_minimum_required(VERSION 3.16)

project(LJMD VERSION 1.0 LANGUAGES C)
add_library(ljmd src/utilities.c src/input.c src/output.c src/force_compute.c src/cleanup.c src/verlet.c src/constants.c)
target_link_libraries(ljmd m)
target_include_directories(ljmd PUBLIC ${CMAKE_SOURCE_DIR}/include )


add_executable(main src/main.c)
target_link_libraries(main PUBLIC ljmd m)