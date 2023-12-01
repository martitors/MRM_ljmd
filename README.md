This package contains refactored MD code for simulating atoms with a Lennard-Jones potential. Its development involved three major steps , which were first  undertaken individualy then merged to the final version :
1. Serial optimisation using compiler flags
2. Distributed memory multiprocessing with MPI
3. Shared memory multiprocessing with 
OpenMP

## Collaborators

```
  
- Regina Mumbi Gachomba   (MumbiGachomba04)
- Martina Torsello        (martitors)
- Mohammad Enayati        (MEnayati)

```

## Code structure

All the optimisation features were implemented on a simplified version of the program that was provided in the beginning.The file ljmd.c was split into multiple files ,i.e. force compute, verlet function, input, output, utilities(pbc, timing), cleanup, contants and header for data structures and prototypes(types.h); These files made up the library mdlib which was linked with the main executable. Simple unit tests for the functions were also set up using the Googletest framework . 
The files and tests were then integrated to the CMakeLists.txt file.
For the purpose of flexibility and portability , the features (MPI,OMP and testing ) are defined as build options in the cmake file .These options can be turned on or off .

```
option(USE_OPENMP "Build with OPENMP support" ON)

option( USE_MPI "Build with MPI support" OFF )

if (USE_OPENMP)
	find_package(OpenMP REQUIRED)
	add_compile_definitions(LJMD_OMP=yes)
endif()


if (USE_MPI) 

    find_package(MPI REQUIRED)
    message( "building with MPI" )
    add_compile_definitions( _MPI=yes )
 
endif()



option(ENABLE_TESTING "Enable building unit tests" ON)

```

## Serial optimisation
 
The program was refactored in the following steps which were performed **sequentially**:

```
NB: Reference: The simplified code without any modifications  took 24.3 seconds for 108 atoms ,10000 steps

1. Adding -O3 compiler optimisation flag inorder to inline the pbc function therefore reducing number of calls
   Time : 17.13 seconds --> 1.4 X faster

2. Adding more compiler flags (-ffast-math -fexpensive-optimizations -msse3)
   Time : 5.5518 seconds --> 4.4 X faster

3. Applying Newton's third law (Fij = -Fji)
   Time  : 2.92 seconds  --> 8.2 times faster

4. Avoiding expensive math (pow,sqrt,division)
   Time : 2.44 seconds  --> 9.8 times faster

```
Below is a graph representation of the speedup and time of the serial optimisation features , run with 108,2916 and 78732 atoms.

![Speedup vs optimisation features](/plots/Serial_bm_1.png)
![Time vs optimisation features](/plots/Serial_bm_2.png)


## Serial Optimisation + MPI

-----TODO : Martina, add a brief explanation on what you did in mpi------

![Speedup MPI](/plots/MPI_bm_1.png)
![Time MPI](/plots/MPI_bm_2.png)


## Serial Optimisation + OpenMP


-----TODO : Mohammad, add a brief explanation on what you did in openmp------

![Speedup OpenMP](/plots/OpenMP_bm_1.png)
![Time OpenMP](/plots/OpenMP_bm_2.png)


## Hybrid (Serial Optimisation + MPI + OpenMP)

![Speedup Hybrid,1 Node](/plots/Hybrid_bm_1.png)
![Speedup Hybrid,More Nodes](/plots/Hybrid_bm_2.png)
![Time Hybrid](/plots/Hybrid_bm_3.png)

## Total Speedup Comparison

![Speedup total](/plots/Total.png)



