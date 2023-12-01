import ctypes
import sys

ljmd_lib = ctypes.CDLL('./libmdlib.so')


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
ljmd_lib.init_mpi_omp.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(mdsys_t)]
ljmd_lib.init_mpi_omp.restype = None

ljmd_lib.force.argtypes = [ctypes.POINTER(mdsys_t)]
ljmd_lib.force.restype = None

ljmd_lib.velverlet_1.argtypes = [ctypes.POINTER(mdsys_t)]
ljmd_lib.velverlet_1.restype = None

ljmd_lib.velverlet_2.argtypes = [ctypes.POINTER(mdsys_t)]
ljmd_lib.velverlet_2.restype = None

ljmd_lib.ekin.argtypes = [ctypes.POINTER(mdsys_t)]
ljmd_lib.ekin.restype = None

ljmd_lib.mpi_bcasts.argtypes = [ctypes.POINTER(mdsys_t)]
ljmd_lib.mpi_bcasts.restype = None

ljmd_lib.mpi_fin.argtypes = [ctypes.POINTER(mdsys_t)]
ljmd_lib.mpi_fin.restype = None

ljmd_lib.cleanup.argtypes = [mdsys_t, ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_void_p)]
ljmd_lib.cleanup.restype = None

def ekin(pysys):
    ljmd_lib.ekin(ctypes.byref(pysys))

def velverlet_1(pysys):
    ljmd_lib.velverlet_1(ctypes.byref(pysys))

def velverlet_2(pysys):
    ljmd_lib.velverlet_2(ctypes.byref(pysys))

def force(pysys):
    ljmd_lib.force(ctypes.byref(pysys))

def init_mpi_omp(argc, argv, pysys):
    ljmd_lib.init_mpi_omp(argc, argv, ctypes.byref(pysys))

def mpi_bcasts(pysys):
    ljmd_lib.mpi_bcasts(ctypes.byref(pysys))

def cleanup(sys, erg, traj):
    ljmd_lib.cleanup(sys, ctypes.byref(erg), ctypes.byref(traj))

def mpi_fin(pysys):
    ljmd_lib.mpi_fin(ctypes.byref(pysys))

def get_a_line(ifile) :
  line = ifile.readline()
  line = line.partition('#')[0]
  line = line.rstrip()
  return line

def read_input(line, restfile, trajfile, ergfile, pysys, nprint):
    if get_a_line(sys.stdin, line):
        exit(1)
    pysys.natoms = int(line)

    if get_a_line(sys.stdin, line):
        exit(1)
    pysys.mass = float(line)

    if get_a_line(sys.stdin, line):
        exit(1)
    pysys.epsilon = float(line)

    if get_a_line(sys.stdin, line):
        exit(1)
    pysys.sigma = float(line)

    if get_a_line(sys.stdin, line):
        exit(1)
    pysys.rcut = float(line)

    if get_a_line(sys.stdin, line):
        exit(1)
    pysys.box = float(line)

    if get_a_line(sys.stdin, restfile):
        exit(1)
    if get_a_line(sys.stdin, trajfile):
        exit(1)
    if get_a_line(sys.stdin, ergfile):
        exit(1)

    if get_a_line(sys.stdin, line):
        exit(1)
    pysys.nsteps = int(line)

    if get_a_line(sys.stdin, line):
        exit(1)
    pysys.dt = float(line)

    if get_a_line(sys.stdin, line):
        exit(1)
    nprint[0] = int(line)

def allocate(pysys):
    pysys.rx = (ct.c_double * pysys.natoms)()
    pysys.ry = (ct.c_double * pysys.natoms)()
    pysys.rz = (ct.c_double * pysys.natoms)()
    pysys.vx = (ct.c_double * pysys.natoms)()
    pysys.vy = (ct.c_double * pysys.natoms)()
    pysys.vz = (ct.c_double * pysys.natoms)()
    pysys.fx = (ct.c_double * pysys.natoms)()
    pysys.fy = (ct.c_double * pysys.natoms)()
    pysys.fz = (ct.c_double * pysys.natoms)()

def output(pysys, erg, traj):
    print("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n" % (pysys.nfi, pysys.temp, pysys.ekin, pysys.epot, pysys.ekin+pysys.epot))
    erg.write("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n" % (pysys.nfi, pysys.temp, pysys.ekin, pysys.epot, pysys.ekin+pysys.epot))
    traj.write("%d\n nfi=%d etot=%20.8f\n"% (pysys.natoms, pysys.nfi, pysys.ekin+pysys.epot))
    for i in range(pysys.natoms):
       traj.write("Ar  %20.8f %20.8f %20.8f\n" % (pysys.rx[i], pysys.ry[i], pysys.rz[i])); 

     

