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
    print("ok1")
    # Pass the arguments and mdsys_t instance by reference to the function
    ft.initialise_mpi_omp(ctypes.c_int(len(arguments)), (ctypes.c_char_p * len(arguments))(*arguments), pysys)
    print("ok2")
    line=None 
    restfile = None
    trajfile = None
    ergfile = None
    nprint = None
    start=None

    if pysys.rank == 0:
        version = 1.0
        print("LJMD version " + str(version) )
        start = timeit.default_timer()
        ft.read_file(line,restfile, trajfile, ergfile, pysys, nprint)
        print("ok3")
    
    ft.bcasts(pysys)
    print("ok4")
    ft.allocate(pysys)
    print("ok5")

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
        if pysys.nfi % nprint == 0 and pysys.rank == 0:
             ft.output(pysys, erg, traj)

        ft.verlet_1(pysys)
        ft.force_compute(pysys)
        ft.verlet_2(pysys)
        ft.ekinetic(pysys)
    
    if pysys.rank == 0:
        print(f"Simulation Done. Run time:{timeit.default_timer() - start:10.3f}s")
    
    ft.clean(pysys, erg, traj)
    ft.mpi_finalise(pysys)