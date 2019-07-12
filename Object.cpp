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

// Esta clase es la implementación de la clase Object

#include "Object.h"
Object Object1(char variable,short int multiplicidad, short int i,short int j)
{

	Object obj = 0; // Inicialmente, nos creamos un entero donde todos los bits son ceros para crearnos un nuevo objeto
	obj = (variable) << 24;
	obj = obj | ((multiplicidad) << 16);
	obj = obj | ((i) << 8);
	obj = obj | ((j) & 0xFF);

        return obj;
}

Object Object2(char variable, int mult2)
{
    
    Object obj = 0;
    obj = (variable) << 24;
    obj = obj | ((mult2) & 0x00FFFFFF);
    
    return obj;
}

Object Object3(char variable,short int mult, short int i)
{
    
    Object obj = 0;
    obj = (variable) << 24;
    obj = obj | ((mult) << 16);
    obj = obj | ((i) & 0x0000FFFF);
    
    return obj;
}

char get_variable(Object obj)
{
	obj = obj >> 24; //Movemos a la derecha 24 bits

	// Se puede hacer directamente obj = obj & 0xFF000000;
	return (char) obj;
}

short int get_mult(Object obj)
{
	obj = (obj & 0x00FF0000) >> 16;
   
	return (short int) obj;
}

int get_extra_mult(Object obj)
{
    obj = (obj & 0x00FFFFFF);
    return (int) obj;
    
}

short int get_extra_i(Object obj)
{
    obj = (obj & 0x0000FFFF);
    return (int) obj;
    
}

short int get_i(Object obj)
{
	obj = ((obj >> 8) & 0x000000FF);
	return (short int) obj;
}

short int get_j(Object obj)
{
	obj = (obj & 0x000000FF);
        return obj;
}
