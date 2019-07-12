#
#    tspSAT-GPU: Simulating an efficient solution to SAT with Tissue P systems on CUDA 
#
#    This simulator is part of the final degree project entitled:
#    "Aceleración de simulaciones de sistemas celulares en soluciones del problema SAT
#     usando GPUs" ("Acceleration of cellular systems simulations on solutions to SAT 
#     problem using GPUs")
#    Jesús Pérez-Carrasco, June 2012, Dpt. Comput. Sci. & A.I. (University of Seville)
#
#    tspSAT-GPU is a subproject of PMCGPU (Parallel simulators for Membrane 
#                                        Computing on the GPU)   
# 
#    Copyright (c) 2012 Jesús Pérez-Carrasco (University of Seville)
#                       Miguel Á. Martínez-del-Amor (RGNC, University of Seville)
#    
#    This file is part of tspSAT-GPU.
#  
#    tspSAT-GPU is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    tspSAT-GPU is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with tspSAT-GPU.  If not, see <http://www.gnu.org/licenses/>.


################################################################################
#
# Build script for project
#
################################################################################

# Add source files here
EXECUTABLE	:= tsp

# Cuda source files (compiled with cudacc)
#CUFILES		:= gpuhybridsolver.cu gpusolver.cu

# CUDA dependency files
#CU_DEPS		:= 

# Cuda source files (compiled with cudacc)
CUFILES_sm_13   :=  GPU_solver.cu gpu_division_evolution_kernel.cu gpu_synchronization_kernel_a.cu gpu_synchronization_kernel_c.cu

# C/C++ source files (compiled with gcc / c++)
CCFILES		:= \
	main.cpp \
	Object.cpp \
	SATcnf.cpp \
	Sequential_solver_a.cpp \
	Sequential_solver_b.cpp \
	Sequential_solver_c.cpp \
	Sequential_solver_d.cpp \
	Sequential_hybrid_solver.cpp \


################################################################################
# Rules and targets

include ../../common/common.mk

INCLUDES += -I../../common/counterslib/inc
LIB	+= -L../../common/counterslib/lib -ltimestat
