# wdmapp_coupling

Adios2-based coupler for XGC and GENE

## Dependencies

- CMake 3.13+
- MPI
- Adios2 2.5
- Kokkos 3.0+
- FFTW 3.3.8+

## Build Instructions

Details instructions for a few systems are available on the wiki.

Another Executing Approach: One would also comment out the BUILDING_TESTING in CMakeFiles.txt included in test folder; 
Assign the full path of testdatas to test_dir in the test_init.cc file; Use "mpirun -np 4 bin/test_init" for the execution.  

## building on summit for gem-xgc coupling

module load cuda/11.1.1  gcc/9.1.0 
