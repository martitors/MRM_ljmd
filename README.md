This package contains refactored MD code for simulating atoms with a Lennard-Jones potential. Its development involved three major steps , which were first  undertaken individualy then merged to the final version :
1. Serial optimisation using compiler flags
2. Distributed memory multiprocessing with MPI
3. Shared memory multiprocessing with 
OpenMP

Benchmarks have been performed on Leonardo cluster.
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
The profiling results are as follows: 

|                                          |       force       |      pbc     |
|------------------------------------------|-------------------|--------------|
|    Serial (without any modifications)    |    67.71%         |    24.45%    |  
|    Add O3                                |    65.97%         |    28.47%    | 
|    More comiler flags                    |    77.93%         |    22.07%    | 
|    Newton's third law                    |    70.80%         |    29.11%    | 
|    Avoid expensive math                  |    70.91%         |    29.09%    |

Below is a graph representation of the speedup and time of the serial optimisation features , run with 108,2916 and 78732 atoms.

![Speedup vs optimisation features](/plots/Serial_bm_1.png)
![Time vs optimisation features](/plots/Serial_bm_2.png)


## Serial Optimisation + MPI

When MPI is active, the workflow to calculate the force is divided among the different tasks. The master processor will read the input file and brodcast the necessary information to the other processors, including the number of atoms and the initial position of the particles. Each processor will calculate the force for a set particles and will collect the potential energy. At the end all the informations will be redirected to the master, which will also provide the print the output.

Below we report the performance of the code following the steps of the previous case. The code has been tested on 1 up to 8 nodes.

  - **108 atoms**: With 8 processors the code reaches the best performance. A further increasing of the number of processor will only make the parallel overhead more significant. 

  - **2916 atoms**: The speedup of the code in this case saturates with 2 nodes.
  
  - **78732 atoms**: Up the 8 nodes, the code is keep scaling almost perfectly.


![Speedup MPI](/plots/MPI_bm_1.png)
![Time MPI](/plots/MPI_bm_2.png)


## Serial Optimisation + OpenMP

Within our computational framework, we use parallel computation to calculate force interactions within a particle set and collect the corresponding potential energy (see the force function inside force_compute file). To do so, using OpenMP's reduction clause can be really useful, especially in scenarios where you need to perform computations across multiple threads and aggregate the results into a single (scalar) value, like calculating potential energy here. But, for force array, we cannot use reduction clause because reduction on arrays is just supported for system with OpenMP version 4.5 and higher. So, to ensure compatibility not only with various OpenMP versions but also with MPI, we divide atoms among different threads and implement a piece of code to use threads to parallelize the reductions. Moreover, see the following plots for scalability. 

![Speedup OpenMP](/plots/OpenMP_bm_1.png)
![Time OpenMP](/plots/OpenMP_bm_2.png)


## Hybrid (Serial Optimisation + MPI + OpenMP)

The code has been structured to allow the user to use both MPI and OpenMP simultaneously. In this way each slave processor will initialize a parallel region within which to subdivide the computation of the forces and potential energy of its own task.

First we test this hybrid parallelization scheme on a single node: we measured the timings using different combination of MPI tasks and Threads, keeping their product equal to 32 (maximum number for a single node). As we can see in the first graph below, for systems with a larger number of particles there is no specific combination of MPI TASK and threads that leads to a better result. On the other hand, for the system with the smallest size this happens with 8 MPI tasks and 4 threads. With a further increasing of the number of threads, performance deteriorates considerably, suggesting an increasing of the parallelization overhead. Smaller size problems indeed can exacerbate issues like workload imbalance due to frequent synchronization points (barriers, reductions), and also sensitivity to cache contention with a system shared simultaneously by an higher number of threads.

![Speedup Hybrid,1 Node](/plots/Hybrid_bm_1.png)

In the next plots we show the speedup and timing results performed on multiple nodes, keeping the same hybrid configuration (8 MPI TASK x 4 Threads) on the single nodes. The system with the lowest number of atoms does not seem to benefit from the increase in nodes as we expected. For the middle system, after 2 nodes the speedup does not increase significantly, instead for the bigger size system performance are linearly improving up to 8 nodes used.

![Speedup Hybrid,More Nodes](/plots/Hybrid_bm_2.png)
![Time Hybrid](/plots/Hybrid_bm_3.png)

## Total Speedup Comparison

Here we report the total speedup comparison between the different version of the program, with different level of optimization and parallelization. For the MPI and OpenMP (on single node) and Hybrid version (on multiple nodes), we report the results of the configuration that leads to the best performance for each system size:
 
|                 |       MPI       |      OPENMP      |                HYBRID                |
|-----------------|-----------------|------------------|--------------------------------------|
|     108 ATOMS   |    16 procs     |    8 threads     |  1 NODES, 8 MPI x 4 Threads per Node |
|    2916 ATOMS   |    16 procs     |    32 threads    |  8 NODES, 8 MPI x 4 Threads per Node |
|   78732 ATOMS   |    32 procs     |    32 threads    |  8 NODES, 8 MPI x 4 Threads per Node |

The advantages of using the hybrid approce increase with the size of the computed system and so reduce significantly the performance for smaller size systems. 
On a single node, the performance obtained paralleling the code with MPI or with OPENMP are very similar between each other for the system we have studied.


![Speedup total](/plots/Total.png)


## Extras

# Python Wrapper

A python wrapper (ljmd.py) for the main function has been created . All indidvidual functions are redefined in funcs_and_types.py which is then imported into ljmd.py. One more C source file was added ;mpi_funcs.c , which encapsulates functions in main that would be implemented only if MPI is on.
NB: **The wrapper works in serial and MPI only. It has not been programmed to handle OpenMP**

The wrapper implementation is in the *Regina_extra* branch

How to run the wrapper:
- Build with cmake as you would do in C
- Copy ljmd.py and funcs_and_types.py into the build file
- Move to the build folder and run the command: ```mpirun -np <num_of_processes> python ljmd.py < ../examples/argon_108.inp```

# Cell list

If the tag **USE CELL** is active, the code will performe the force computation using the cell list optimization. In the *cell.c* file all the functions needed are collected, including the building, the updating and the freeing of the total grid.
The size of the cell has been chosen as 2 time the cutoff, but changing the ratio between this two is an easy and user-friendly operation (just changing the first value of the building function):
```
double cell_rcut_ratio = 2.0;
```
Each system has a *clist* and *plist* which respectivelly contains the list of the cells and the list of pairs of close cells. In the building function, we store the number of atoms inside each cells, the center coordinates of the cell and the global indexes of the atoms included.
In the updating function, we restore the clist and insert each atoms in the cell which contains it.
Finally, we compute the force splitting the global loop in two cases: first we compute the force between atoms which belong to the same cell, then using the pair list we compute the force between particles in different cells.
Since for the input data available the time for updating the cells are considerably smaller than the total time, this operation has been performed at every step.

As we can see in the plot below, the optimization due to the cell list confguration is evident when we are dealing with large systems of particles(the blue line represents the classic approch while the green line is the one obtained using the cell list).


![Cell_list](/plots/Cell.png)

#Morse potential

Note that: this part is not merge in the "main" branch, you can find it in the "Morse" branch.

The code implements two different types of interatomic potentials: Lennard-Jones potential and Morse potential. The distinction between which potential to use is determined by a flag provided in the command-line argument. Indeed, when the flag is set to '0' in the command-line argument, the code executes calculations using the Lennard-Jones potential. On the other hand, when the flag is set to '1' in the command-line argument, the code switches to utilize the Morse potential. The implementation also involves OpenMP and MPI considerations. This approach provides flexibility and versatility to the code, enabling users to switch between different potential models (Lennard-Jones or Morse) based on the requirements of their simulations or calculations through a simple command-line flag selection. Moreover, the following is the scalability based on the OpenMP for Morse potential:

for 108 atoms:

threads   |          time
---------------------------
 1        |          1.450
 2        |          0.611
 4        |          0.366
 8        |          0.221
 16       |          0.286
 32       |          0.536

for 2916 atoms:

threads    |          time
----------------------------
 1         |          47.399
 2         |          24.220
 4         |          12.466
 8         |          6.618
 16        |          3.506
 32        |          2.070

Similar to the behavior observed with the Lennard-Jones potential, the code achieves optimal performance with 8 threads when handling 108 atoms. Additionally, scalability remains linear for 2916 atoms.

