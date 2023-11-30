cmake_minimum_required(VERSION 3.16)

project(LJMD VERSION 1.0 LANGUAGES C CXX)

option(USE_OPENMP "Build with OPENMP support" ON)

option( USE_MPI "Build with MPI support" ON )

if (USE_OPENMP)
	find_package(OpenMP REQUIRED)
	add_compile_definitions(LJMD_OMP=yes)
endif()


if (USE_MPI) 

    find_package(MPI REQUIRED)
    message( "building with MPI" )
    add_compile_definitions( _MPI=yes )
    #include_directories( ${MPI_C_INCLUDE_PATH})
    #include_directories( ${MPI_CXX_INCLUDE_PATH} )

endif()

######################## library ##########################################
add_library(mdlib SHARED 
            src/utilities.c 
            src/input.c 
            src/output.c 
            src/force_compute.c
            src/cleanup.c 
            src/verlet.c 
            src/constants.c 
            src/allocate.c 
            src/mpi_funcs.c)
target_compile_options(mdlib PRIVATE -O3 -ffast-math -fexpensive-optimizations -msse3)
if (USE_OPENMP)
	target_link_libraries(mdlib OpenMP::OpenMP_C)
endif()
if (USE_MPI) 
    target_link_libraries(mdlib MPI::MPI_C)
endif()
target_link_libraries(mdlib m)
target_include_directories(mdlib PUBLIC ${CMAKE_SOURCE_DIR}/include )

##################### main #################################

add_executable(main src/main.cpp)
target_compile_options(main PUBLIC -O3 -ffast-math -fexpensive-optimizations -msse3)

if (USE_OPENMP)
        target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)
        target_compile_options(main PUBLIC -fopenmp)
endif()

if (USE_MPI) 
    target_link_libraries(main PUBLIC MPI::MPI_CXX)
endif()

target_link_libraries(main PRIVATE mdlib)
target_include_directories(main PUBLIC ${CMAKE_SOURCE_DIR}/include )


########################## test ###############################
option(ENABLE_TESTING "Enable building unit tests" ON)
if(ENABLE_TESTING)

  enable_testing()
  enable_language(CXX)

  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  include(GTest)

  add_executable(test_force tests/test_force.cpp)
  target_include_directories(test_force PRIVATE ${CMAKE_SOURCE_DIR}/src)
  
  if (USE_OPENMP)
	  target_link_libraries(test_force PUBLIC OpenMP::OpenMP_CXX)
	  target_compile_options(test_force PUBLIC -fopenmp)
  endif()

    if (USE_MPI) 
    target_link_libraries(test_force PUBLIC MPI::MPI_CXX)
    endif()

  target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Force COMMAND test_force)



  add_executable(test_verlet tests/test_verlet.cpp)
  target_include_directories(test_verlet PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_verlet PRIVATE GTest::GTestMain GTest::GTest mdlib m)
  add_test(NAME Verlet COMMAND test_verlet)

  add_executable(test_utils tests/test_utils.cpp)
  target_include_directories(test_utils PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_utils PRIVATE GTest::GTestMain GTest::GTest mdlib)
 add_test(NAME Utils COMMAND test_utils)

  add_executable(test_ekin tests/test_ekin.cpp)
  target_include_directories(test_ekin PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_ekin PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Ekin COMMAND test_ekin)

  add_executable(test_read tests/test_read.cpp)
  target_include_directories(test_read PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(test_read PRIVATE GTest::GTestMain GTest::GTest mdlib)
  add_test(NAME Read COMMAND test_read)

endif()

