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

#include "SATcnf.h"
#include "GPU_solver.h"
#include "Sequential_hybrid_solver.h"
#include "Sequential_solver_a.h"
#include "Sequential_solver_b.h"
#include "Sequential_solver_c.h"
#include "Sequential_solver_d.h"

#include <string>
#include <iostream>
using namespace std;



int main(int argc, char* argv[])
{

	SATcnf scnf;
	int mode = 5;
        //int sequential_solver=0;
	string input_file;

	bool solution;
	char soltype = '\0';

	while((soltype = getopt(argc,argv,"m:f:s:h:?"))!=-1)
	{

		switch(soltype)
		{

		case 'f':
			input_file = optarg;
			break;	
                case 'm':
                        mode = atoi(optarg);
                        break;
		case 'h':

		case '?':
			cout << "Copyright (C) 2012 J. Pérez-Carrasco, M.A. Martínez-del-Amor" << endl <<
			"This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome to redistribute it under certain conditions" << endl << endl;	            
			
			cout << "The parameters should be the following ones:" <<
					"-f (input file) " <<
					"-m (sequential: 0,1,2(recom.),3; hybrid 4(deprec.); parallel 5) " <<
					"-h (this help);"
			<< endl;
			cout << "TSPCUDASAT Version 1.0" << endl;

			return 0;


		}

	}


	cout << "Reading file..." << endl;

	if(!scnf.parse(input_file.c_str()))

		return -1;

	cout << "without any problem reading the file" << endl;
	cout << "Information of the problem" << endl;

	Object* obj = scnf.getcnf(); // Obtenemos el problema SAT en forma normal conjuntiva

	cout << "Instance size: N =" << scnf.getN() << "M =" << scnf.getM() << endl;

	//Representamos la codificacion de todos los objetos que son entrada del problema SAT

	cout << "objects:" ;

		for(int i=0; i<scnf.getT();i++)
		{

			cout << get_variable(obj[i]) << get_mult(obj[i]) << get_i(obj[i]) << "," << get_j(obj[i]) << "" ;

			cout << endl;

		}

	switch (mode)
	{

		case 0:
			solution = Sequential_solver_a(scnf.getN(),scnf.getM(),scnf.getT(), scnf.getcnf());
			break;
                        
                case 1:
                        solution = Sequential_solver_b(scnf.getN(),scnf.getM(),scnf.getT(), scnf.getcnf());
                        break;
                        
                case 2:
			solution = Sequential_solver_c(scnf.getN(),scnf.getM(),scnf.getT(), scnf.getcnf());
			break;
                        
                case 3:
			solution = Sequential_solver_d(scnf.getN(),scnf.getM(),scnf.getT(), scnf.getcnf());
			break;
                        
                case 4:   
                        solution = Sequential_hybrid_solver(scnf.getN(),scnf.getM(),scnf.getT(), scnf.getcnf());
                        break;
                        
		case 5:
			solution = GPU_solver(scnf.getN(),scnf.getM(),scnf.getT(), scnf.getcnf());
			break;

		default:

			cout << "Wrong type. Use 'h' or '?' for help" << endl;

			return 0;
	}
        
	cout << "The Tissue P System response is:" << solution << endl;

        
	return 0;

}
