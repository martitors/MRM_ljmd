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

All the optimisation features were implemented on a simplified version of the program that was provided in the beginning.The file ljmd.c was split into multiple files ,i.e. force compute, verlet function, input, output, utilities(pbc, timing), cleanup, contants and header for data structures and prototypes(types.h); These files made up the library mdlib which was linked with the main executable. Simple unit tests for the functions were also set up using the Googletest framework . The files and tests were then integrated to the CMakeLists.txt file.






