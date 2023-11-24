cmake_minimum_required(VERSION 3.16)

project(LJMD VERSION 1.0 LANGUAGES C CXX)

add_library(mdlib src/utilities.c src/input.c src/output.c src/force_compute.c src/cleanup.c src/verlet.c src/constants.c)

target_link_libraries(mdlib m)
target_include_directories(mdlib PUBLIC ${CMAKE_SOURCE_DIR}/include)


add_executable(main src/main.cpp)
target_link_libraries(main PUBLIC mdlib m)


option(ENABLE_TESTING "Enable building unit tests" ON)
if(ENABLE_TESTING)
  enable_testing()
  enable_language(CXX)

  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  include(GTest)
  add_executable(test_verlet tests/test_verlet.cpp)
  target_include_directories(test_verlet PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_verlet PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Verlet COMMAND test_verlet)

  add_executable(test_utils tests/test_utils.cpp)
  target_include_directories(test_utils PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_utils PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Utils COMMAND test_utils)
endif()

