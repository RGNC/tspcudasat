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

#ifndef __SATCNF
#define __SATCNF

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "Object.h"



class SATcnf
{
	private:
	unsigned int n;
	unsigned int m;
	unsigned int t;

	Object* cnf;

	public:

		SATcnf()
		{
			n=0;
			m=0;
			t=0;
			cnf = NULL;
		}

		~SATcnf()
		{
			if(cnf != NULL)
				delete cnf;
		}


// Se crea una instancia cnf(en forma normal conjuntiva) a partir del archivo de entrada.

                bool parse(const char* file);

                // Numero de variables en el problema SAT
                unsigned int getN()
                {
                        return n;
	
                }
                // Numero de clausulas en el problema SAT
                unsigned int getM()
                {
			return m;
                }

                //Numero de objetos en la codificacion (del fichero de entrada)
                unsigned int getT()
                {
			return t;
        	}

                //Obtenemos la formula CNF en formato comprimido

                Object* getcnf()
                {
                        return cnf;
                }
};

#endif
