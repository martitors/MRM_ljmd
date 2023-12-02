import ctypes
#from ctypes import byref
import funcs_and_types as ft
import sys
import timeit


if __name__ == "__main__":
    # Extract command-line arguments passed to the Python script
    arguments = [ctypes.c_char_p(arg.encode('utf-8')) for arg in sys.argv[1:]]  # Encode string arguments as bytes
    
    # Create an instance of mdsys_t
    pysys = ft.MDSys()

    # İnit mpi
    ft.initialise_mpi_omp(ctypes.c_int(len(arguments)), (ctypes.c_char_p * len(arguments))(*arguments), pysys) 
    start = 0.0
    ergfile = "ergf"
    trajfile = "trajf"
    restfile = "restf"
    nprint=0

    if pysys.rank == 0:
        version = 1.0
        print("LJMD version " + str(version) )
        start = timeit.default_timer()
        with  sys.stdin as input_file:
                restfile, trajfile, ergfile, nprint = ft.read_input(input_file, pysys)


    ft.bcasts(pysys)
    ft.allocate(pysys)
  
    if pysys.rank == 0:
        # Read restart
        restPath = f"../examples/{restfile}"

        try:
            with open(restPath, "r") as fp:
                for i in range(pysys.natoms):
                    pysys.rx[i], pysys.ry[i], pysys.rz[i] = map(float, fp.readline().split())
                for i in range(pysys.natoms):
                    pysys.vx[i], pysys.vy[i], pysys.vz[i] = map(float, fp.readline().split())
        except IOError as e:
            print("cannot read restart file:", e)
            sys.exit(3)
    
    pysys.nfi=0
    ft.force_compute(pysys)
 

    if pysys.rank == 0:
        ft.ekinetic(pysys)
 
        ergPath = f"../examples/{ergfile}"
        trajPath = f"../examples/{trajfile}"

        with open(ergPath, "w") as erg, open(trajPath, "w") as traj:
            print(f"Startup time: {timeit.default_timer() - start:10.3f}s") ## ------time diff taken
            print(f"Starting simulation with {pysys.natoms} atoms for {pysys.nsteps} steps.")
            print("     NFI            TEMP            EKIN                 EPOT              ETOT")
            ft.output(pysys, erg, traj)
            start = timeit.default_timer()

    for pysys.nfi in range(1, pysys.nsteps + 1): 
        if  pysys.rank == 0 and pysys.nfi % nprint == 0:
            ergPath = f"../examples/{ergfile}"
            trajPath = f"../examples/{trajfile}"
            with open(ergPath, "w") as erg, open(trajPath, "w") as traj:
                 ft.output(pysys, erg, traj)

        ft.verlet_1(pysys)
        ft.force_compute(pysys)
        ft.verlet_2(pysys)
        ft.ekinetic(pysys)

    
    if pysys.rank == 0:
        print(f"Simulation Done. Run time:{timeit.default_timer() - start:10.3f}s")
    
    ergPath = f"../examples/{ergfile}"
    trajPath = f"../examples/{trajfile}"
    if pysys.rank == 0:
        with open(ergPath, "w") as erg, open(trajPath, "w") as traj:
             ft.close_files(pysys, erg, traj)
    ft.mpi_finalise(pysys)