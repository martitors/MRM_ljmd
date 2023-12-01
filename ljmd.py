import ctypes

ljmd = ctypes.CDLL('./libmdlib.so')

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

ljmd.init_mpi_omp.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p), ctypes.POINTER(mdsys_t)]
ljmd.init_mpi_omp.restype = None
