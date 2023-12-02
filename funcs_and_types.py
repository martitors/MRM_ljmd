import ctypes
import sys

ljmd_lib = ctypes.CDLL('./libmdlib.so') #,ctypes.RTLD_GLOBAL

class MDSys(ctypes.Structure):
    _fields_ = [
        ("natoms", ctypes.c_int),
        ("nfi", ctypes.c_int),
        ("nsteps", ctypes.c_int),
        ("dt", ctypes.c_double),
        ("mass", ctypes.c_double),
        ("epsilon", ctypes.c_double),
        ("sigma", ctypes.c_double),
        ("box", ctypes.c_double),
        ("rcut", ctypes.c_double),
        ("ekin", ctypes.c_double),
        ("epot", ctypes.c_double),
        ("temp", ctypes.c_double),
        ("rx", ctypes.POINTER(ctypes.c_double)),
        ("ry", ctypes.POINTER(ctypes.c_double)),
        ("rz", ctypes.POINTER(ctypes.c_double)),
        ("vx", ctypes.POINTER(ctypes.c_double)),
        ("vy", ctypes.POINTER(ctypes.c_double)),
        ("vz", ctypes.POINTER(ctypes.c_double)),
        ("fx", ctypes.POINTER(ctypes.c_double)),
        ("fy", ctypes.POINTER(ctypes.c_double)),
        ("fz", ctypes.POINTER(ctypes.c_double)),
        ("nthreads", ctypes.c_int),
        ("rank", ctypes.c_int),
        ("npes", ctypes.c_int),
        ("cx", ctypes.POINTER(ctypes.c_double)),
        ("cy", ctypes.POINTER(ctypes.c_double)),
        ("cz", ctypes.POINTER(ctypes.c_double)),
    ]

#prototypes for c functions
ljmd_lib.init_mpi_omp.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(MDSys)]
ljmd_lib.init_mpi_omp.restype = None
 

ljmd_lib.force.argtypes = [ctypes.POINTER(MDSys)]
ljmd_lib.force.restype = None

ljmd_lib.velverlet_1.argtypes = [ctypes.POINTER(MDSys)]
ljmd_lib.velverlet_1.restype = None

ljmd_lib.velverlet_2.argtypes = [ctypes.POINTER(MDSys)]
ljmd_lib.velverlet_2.restype = None

ljmd_lib.ekin.argtypes = [ctypes.POINTER(MDSys)]
ljmd_lib.ekin.restype = None

ljmd_lib.mpi_bcasts.argtypes = [ctypes.POINTER(MDSys)]
ljmd_lib.mpi_bcasts.restype = None

ljmd_lib.mpi_fin.argtypes = [ctypes.POINTER(MDSys)]
ljmd_lib.mpi_fin.restype = None

# ljmd_lib.cleanup.argtypes = [MDSys, ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_void_p)]
# ljmd_lib.cleanup.restype = None


# ljmd_lib.read_input.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,ctypes.c_char_p, ctypes.POINTER(MDSys), ctypes.POINTER(ctypes.c_int)]
# ljmd_lib.read_input.restype = None
def ekinetic(pysys):
    ljmd_lib.ekin(ctypes.byref(pysys))

def verlet_1(pysys):
    ljmd_lib.velverlet_1(ctypes.byref(pysys))

def verlet_2(pysys):
    ljmd_lib.velverlet_2(ctypes.byref(pysys))

def force_compute(pysys):
    ljmd_lib.force(ctypes.byref(pysys))

def initialise_mpi_omp(argc, argv, pysys):
    ljmd_lib.init_mpi_omp(argc, argv, ctypes.byref(pysys))

def bcasts(pysys):
    ljmd_lib.mpi_bcasts(ctypes.byref(pysys))

def close_files(sys, erg, traj):
    #ljmd_lib.cleanup(sys, ctypes.byref(erg), ctypes.byref(traj))
    erg.close()
    traj.close()

def mpi_finalise(pysys):
    ljmd_lib.mpi_fin(ctypes.byref(pysys))

# def read_file(line,restfile, trajfile, ergfile, pysys, nprint):
#     ljmd_lib.read_input(line,restfile, trajfile, ergfile,ctypes.byref(pysys),ctypes.byref(nprint))

def get_a_line(ifile):
    line = ifile.readline()
    line = line.split('#')[0]
    line = line.strip()
    return line

def read_input(input_file, pysys):
    pysys.natoms = int(get_a_line(input_file))
    pysys.mass = float(get_a_line(input_file))
    pysys.epsilon = float(get_a_line(input_file))
    pysys.sigma = float(get_a_line(input_file))
    pysys.rcut = float(get_a_line(input_file))
    pysys.box = float(get_a_line(input_file))
    restfile = get_a_line(input_file)
    trajfile = get_a_line(input_file)
    ergfile = get_a_line(input_file)
    pysys.nsteps = int(get_a_line(input_file))
    pysys.dt = float(get_a_line(input_file))
    nprint = int(get_a_line(input_file))

    return restfile, trajfile, ergfile, nprint


def allocate(pysys):
    pysys.rx = (ctypes.c_double * pysys.natoms)()
    pysys.ry = (ctypes.c_double * pysys.natoms)()
    pysys.rz = (ctypes.c_double * pysys.natoms)()
    pysys.vx = (ctypes.c_double * pysys.natoms)()
    pysys.vy = (ctypes.c_double * pysys.natoms)()
    pysys.vz = (ctypes.c_double * pysys.natoms)()
    pysys.fx = (ctypes.c_double * pysys.natoms)()
    pysys.fy = (ctypes.c_double * pysys.natoms)()
    pysys.fz = (ctypes.c_double * pysys.natoms)()

def output(pysys, erg, traj):
    print("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n" % (pysys.nfi, pysys.temp, pysys.ekin, pysys.epot, pysys.ekin+pysys.epot))
    erg.write("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n" % (pysys.nfi, pysys.temp, pysys.ekin, pysys.epot, pysys.ekin+pysys.epot))
    traj.write("%d\n nfi=%d etot=%20.8f\n"% (pysys.natoms, pysys.nfi, pysys.ekin+pysys.epot))
    for i in range(pysys.natoms):
       traj.write("Ar  %20.8f %20.8f %20.8f\n" % (pysys.rx[i], pysys.ry[i], pysys.rz[i])); 

     

