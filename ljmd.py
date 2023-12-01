import ctypes
import sys


if __name__ == "__main__":
    # Extract command-line arguments passed to the Python script
    arguments = [ctypes.c_char_p(arg.encode('utf-8')) for arg in sys.argv[1:]]  # Encode string arguments as bytes
    
    # Create an instance of mdsys_t
    pysys = mdsys_t()
    
    # Pass the arguments and mdsys_t instance by reference to the function
    initialize_mpi_omp(ctypes.c_int(len(arguments)), (ctypes.c_char_p * len(arguments))(*arguments), pysys)

    if pysys.
