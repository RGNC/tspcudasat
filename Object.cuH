/*
    tspSAT-GPU: Simulating an efficient solution to SAT with Tissue P systems on CUDA 

    This simulator is part of the final degree project entitled:
    "Aceleración de simulaciones de sistemas celulares en soluciones del problema SAT
     usando GPUs" ("Acceleration of cellular systems simulations on solutions to SAT 
     problem using GPUs")
    Jesús Pérez-Carrasco, June 2012, Dpt. Comput. Sci. & A.I. (University of Seville)

    tspSAT-GPU is a subproject of PMCGPU (Parallel simulators for Membrane 
                                        Computing on the GPU)   
 
    Copyright (c) 2012 Jesús Pérez-Carrasco (University of Seville)
                       Miguel Á. Martínez-del-Amor (RGNC, University of Seville)
    
    This file is part of tspSAT-GPU.
  
    tspSAT-GPU is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    tspSAT-GPU is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with tspSAT-GPU.  If not, see <http://www.gnu.org/licenses/>. */
	
#ifndef OBJECT_CUH
#define	OBJECT_CUH

#include "cutil_inline.h"

// Implementación de la clase Object para CUDA

#include "Object.h"

typedef unsigned int Object;

__device__ inline Object d_set_Object1(char variable,short int multiplicidad, short int i,short int j)
{

	Object obj = 0; // Inicialmente, creamos un entero donde todos los bits son ceros para crearnos el nuevo objeto
	obj = (variable) << 24;
	obj = obj | ((multiplicidad) << 16);
	obj = obj | ((i) << 8);
	obj = obj | ((j) & 0xFF);

        return obj;
}

__device__ inline Object d_set_Object2(char variable, int mult2)
{
    
    Object obj = 0;
    obj = (variable) << 24;
    obj = obj | ((mult2) & 0x00FFFFFF);
    
    return obj;
}

__device__ inline Object d_set_Object3(char variable,short int mult, short int i)
{
    
    Object obj = 0;
    obj = (variable) << 24;
    obj = obj | ((mult) << 16);
    obj = obj | ((i) & 0x0000FFFF);
    
    return obj;
}

__device__ inline char d_get_variable(Object obj)
{
	obj = obj >> 24; 

	// Se puede hacer directamente obj = obj & 0xFF000000;
	return (char) obj;
}

__device__ inline Object d_setVar1(Object obj, char var) {
	return obj=(obj&0xFFFFFF)|(var)<<24;
}

__device__ inline short int d_get_mult(Object obj)
{
	obj = (obj & 0x00FF0000) >> 16;
   
	return (short int) obj;
}

__device__ inline short int d_set_mult(Object obj, short int mult) {
	return obj=(obj&0xFF00FFFF)|(mult)<<16;
}

__device__ inline int d_get_extra_mult(Object obj)
{
    obj = (obj & 0x00FFFFFF);
    return (int) obj;
    
}

__device__ inline short int d_set_extra_mult(Object obj, short int mult) {
	return obj=(obj&0xFF000000)|mult;
}

__device__ inline short int d_get_extra_i(Object obj)
{
    obj = (obj & 0x0000FFFF);
    return (int) obj;
    
}

__device__ inline short int d_set_extra_i(Object obj, short int i) {
	return obj=(obj&0xFFFF0000)|i;
}

__device__ inline short int d_get_i(Object obj)
{
	obj = ((obj >> 8) & 0x000000FF);
	return (short int) obj;
}

__device__ inline Object d_setVar2(Object obj, char i) {
	return obj=(obj&0xFFFF00FF)|(i)<<8;
}

__device__ inline short int d_get_j(Object obj)
{
	obj = (obj & 0x000000FF);
        return obj;
}

__device__ inline Object d_setVar3(Object obj, char j) {
	return obj=(obj&0xFFFFFF00)|(j);
}


#endif	/* OBJECT_CUH */

