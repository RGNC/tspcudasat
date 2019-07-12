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
#include "Object.h"

using namespace std;

bool SATcnf::parse(const char* file)
{
	char string[512];
	char car= 'c';
        char c;
        
	int variable = 0;
	int clausula = 0;

	// Todos los atributos privados, lo inicializamos a cero o vacio.

	char * c_cnf;

	n = 0;
	m = 0;
	t = 0;

	FILE * f = NULL;

	// Leemos el fichero de entrada

	f = fopen(file, "ro"); // Abrimos el archivo de solo lectura

	if(f == NULL)
	{
		perror("Can not open the input:file");
		return false;
	}

	while(!feof(f) && !ferror(f))
	{
		c = fgetc(f);

		if(c == car)  // Obviamos comentarios
		{
			while(c!='\n' && c != EOF)
			{
				c= fgetc(f);
			}

		}

                else if(c=='p') // Pasamos a leer la instancia del problema
		{
			fscanf(f,"%s",&string);

			// Se comprueba si la siguiente palabra es cnf
			if((strcmp(string, "cnf")!=0))
			{
				cerr << "Not a valid input of cnf" << endl;
				fclose(f);
				return false;
			}

			// Leemos n y m
			fscanf (f, "%d""%d", &n, &m);

			if(n<=0 || m<=0)
			{
				cerr << "Invalid value of n or m" << endl;
				fclose(f);
				return false;
			}

			//Inicializamos los datos

			c_cnf = new char[n*m];

			for(int i=0;i<n*m;i++)
			{
				c_cnf[i]='\0';
			}
			clausula++;

			while ((clausula<=m) && (!feof(f) && !ferror(f)))
			{
				fscanf(f,"%d",&variable);

				if (variable == 0)
				{
					clausula++; //Incrementamos el número de clausulas, para leer los siguientes objetos de otras clausulas
                                }
				else if((variable>0)&&(variable <=n))
				{
					c_cnf[((variable-1)* m) + (clausula-1)] = '+';
					t++;  // Incrementamos el contador de objetos
				}
				else if ((variable <0)&&(variable>=((-1)*n)))
				{
					variable *= -1;
					c_cnf[((variable-1)*m) + (clausula-1)] = '-';
					t++; // Incrementamos el contador de objetos
                                        
                                        
				}
			}
		}
	}

	fclose(f);

        cout << "T vale " << t << endl;
	// Inicializamos el array de entrada con las n variables and m clausulas, indicadas en el problema

	int indice = 0;
        int i;
        int j;
	cnf = new Object[t];

	for(j=0;j<m;j++)
	{

		for(i=0;i<n;i++)
		{
			if (c_cnf[i*m + j] == '+')
			{
				cnf[indice++] = Object1('x',0, i+1, j+1);
			}

                        else if(c_cnf[i*m + j] == '-')
                        {
				cnf[indice++] = Object1('y',0, i+1, j+1);
			}
		}
	}
        
        delete c_cnf;
        
	return true;

}





