readme.txt: This file explains the code and the structure of the algorithm.

I. Motivation

This simulator is able to simulate a family of tissue-like P systems with cell division solving SAT problem in linear time. The input of the algorithm is given by a DIMACS CNF file, which codifies an instance of SAT problem (a CNF formula) to be automatically converted to an input multiset. It is based on the stages detected during the computation of the P systems:

1. Generation.
2. Exchange.
3. Synchronization.
4. Check-in.
5. Output.


II. Installation: 

1. Install the CUDA SDK version 4.X.

2. Install the counterslib: extract the file counterslib.tar.gz inside the common folder of the CUDA SDK.

3. Copy the folder into the source folder of the CUDA SDK, and type "make".


III. Usage:

Type ./tsp -h to list the different options.
* A sequential simulation: ./tsp -m 2 -f file.cnf
* A parallel simulation on the GPU: ./tsp -m 4 -f file.cnf


IV. Source:

The objective of each file is the following:

main.cpp: contains the main function, calling the different algorithms.

SATcnf.cpp, SATcnf.h: parse the input files.

sequential_solver_x.cpp/.h: the sequential solvers. Sequential_solver_c.cpp/.h is the most updated version. Rest of files were used to test several optimization approaches.

sequential_hybrid_solver.cpp/.h: deprecated.

object.cpp, object.cuH, object.h: implementation of objects for both CPU and GPU.

GPU_solver.cu/.h, gpu_division_evolution_kernel.cu, gpu_checkout_kernel.cu, gpu_synchronization_kernel_x.cu: host and kernels of the GPU implementation. gpu_synchronization_kernel_d.cu is the most updated version. Rest of files were used to test several optimization approaches.

V. Warning:

Many function names and comments are in Spanish! This simulator was made as a final project of the computer engineering degree (similar to master thesis).

/*$Id: readme.txt 2012-12-11 16:02:12 mdelamor $*/
