/*
    tspSAT-GPU: Simulating an efficient solution to SAT with Tissue P systems on CUDA 

    This simulator is part of the final degree project entitled:
    "Aceleraci�n de simulaciones de sistemas celulares en soluciones del problema SAT
     usando GPUs" ("Acceleration of cellular systems simulations on solutions to SAT 
     problem using GPUs")
    Jes�s P�rez-Carrasco, June 2012, Dpt. Comput. Sci. & A.I. (University of Seville)

    tspSAT-GPU is a subproject of PMCGPU (Parallel simulators for Membrane 
                                        Computing on the GPU)   
 
    Copyright (c) 2012 Jes�s P�rez-Carrasco (University of Seville)
                       Miguel �. Mart�nez-del-Amor (RGNC, University of Seville)
    
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

#include <stdio.h>

#include "Sequential_hybrid_solver.h"
#include "Object.h"
#include <iostream>
using namespace std;


void printm2(Object* multiset, unsigned int number_membranes, unsigned int T, unsigned int N) {
    
        int TM2 = 2*N + 3 + 1+ T;
	cout << "Number of membranes: " << number_membranes << endl;

        
	cout << "Multisets: ";
	for (int i=0; i<number_membranes; i++) {
		cout << "|"<< i << "|: ";

                
                for(int c=0;c < (2*N)+2;c++)
                {
                        int o=(i*TM2)+c;
                        cout << get_variable(multiset[o]) << "_" << get_i(multiset[o]) << "," << get_j(multiset[o]) << "*" << get_mult(multiset[o]) << "; ";
                }
                
                cout << get_variable(multiset[i*TM2+2*N+2]) << "_" << get_extra_i(multiset[i*TM2+2*N+2]) << "*" << get_mult(multiset[i*TM2+2*N+2])<< ";";
                cout << get_variable(multiset[i*TM2+2*N+3]) << "*" << get_extra_mult(multiset[i*TM2+2*N+3]) << "; ";
                
		for (int j=0; j<T; j++) {
			int o=(i*TM2)+(2*N)+4+j;
			cout << get_variable(multiset[o]) << "_" << get_i(multiset[o]) << "," << get_j(multiset[o]) << "; ";

		}
		if (i%8==7)
                cout << endl;
	}
	cout << endl;
}

void printm1(Object* multiset) {
    
    cout << "Number of membranes: " << 1 << endl;

            cout << get_variable(multiset[0]) << "_" << get_i(multiset[0]) << "," << get_j(multiset[0]) << "*" << get_mult(multiset[0]) << "; ";
            cout << get_variable(multiset[1]) << "_" << get_i(multiset[1]) << "," << get_j(multiset[1]) << "*" << get_mult(multiset[1]) << "; ";   
            cout << get_variable(multiset[2]) << "_" << get_extra_i(multiset[2]) << "*" << get_mult(multiset[2])<< ";";
            cout << get_variable(multiset[3]) << "_" << get_i(multiset[3]) << "," << get_j(multiset[3]) << "*" << get_mult(multiset[3]) << "; ";
            cout << get_variable(multiset[4]) << "_" << get_i(multiset[4]) << "," << get_j(multiset[4]) << "*" << get_mult(multiset[4]) << "; ";
            
            
        cout << endl;
}

void printcnf(Object* multiset,unsigned int T) {
	
        cout << "T=" << T << endl << "La formula es: " << endl;
	

        for (int o=0; o<T; o++) 
        {
			

                cout << get_variable(multiset[o]) << get_i(multiset[o]) << "," << get_j(multiset[o]) << "*" << get_mult(multiset[o]) << " ";

        }
		
	
	cout << endl;
}




unsigned int division(unsigned int NumMemb,int N,int M,int T, Object* multiset2)
{

	Object obj_actual;
	char variable;
	int TM2 = (2*N) + 3 + 1 + T;
        //cout << "TM2:" << TM2 << endl;


	for(int memb=0;memb<NumMemb; memb++)
	{
		
           
		for(int obj_indice = 0;obj_indice < TM2;obj_indice++)
		{

			obj_actual = multiset2[(memb * TM2)+ obj_indice];
			multiset2[((memb+NumMemb)* TM2) +obj_indice] = obj_actual;

		}
        }
		
                NumMemb = NumMemb * 2;
                //printm2(multiset2,NumMemb,T,N);
                
		return NumMemb;
               
}



void evolution_m2(unsigned int NumMemb,unsigned int T,short int N, Object* multiset2,int mult_a, int e)
{


	int memb = 0;
	Object obj_actual;
	short int i;
	short int j;
	short int multiplicidad;
	char variable = '\0';
	int TM2 = (2*N + 3 + 1 + T);
        int inicioa = (2*N)+3;
    
        
        for(memb;memb<(NumMemb/2);memb++)
         {
		
                        
                obj_actual = multiset2[(memb*TM2)+inicioa];
                variable = get_variable(obj_actual);
		if(variable == 'a')
		{
                     int aux = mult_a - 1;
                     multiset2[(memb*TM2)+inicioa] = Object2('a', aux);
                     i = N-aux;
                     multiset2[(memb*TM2)+(e-1)] = Object1('T',0,i,0); // el contador e se le resta uno porque empezamos el array en 0
                            
                }
                
         }

        for(memb = (NumMemb/2);memb<NumMemb;memb++)
         {
		
                obj_actual = multiset2[(memb*TM2)+inicioa];
                variable = get_variable(obj_actual);
                if(variable == 'a')
		{
                        int aux = mult_a - 1;
                        multiset2[(memb*TM2)+inicioa] = Object2('a', aux);
                        i = N-aux;
                        multiset2[(memb*TM2)+(e-1)] = Object1('F',0,i,0);
				
                }

		
          }
              
          //printm2(multiset2,NumMemb,T,N);
}

void evolution_m1(short int N,short int M,unsigned short T,Object* multiset1)
{
	Object obj_actual;
	char variable;
	short int multiplicidad = N;
	short int i;
        int TM1 = 5;

	for(int pos= 0;pos<TM1;pos++)
	{
                obj_actual = multiset1[pos];
		variable = get_variable(obj_actual);
		if(variable == 'b' || variable == 'c')
		{
			i = get_i(obj_actual);
			multiset1[pos] = Object1(variable, multiplicidad, N+1,0);
                }
                if (variable == 'd')
                {
                    multiset1[pos] = Object3(variable,multiplicidad,N+1);
                }

	}
}


void intercambio (Object* multiset1, Object* multiset2,unsigned int membranas_actuales2,int T,int N,int M)
{
	Object obj_actual;
        char variable = '\0';
        int TM1 = 5;
        int TM2 = (2*N)+3+1+T;
        
        // Actualizamos el multiset1
	for(int pos= 0;pos < TM1;pos++)
	{
                obj_actual = multiset1[pos];
		variable = get_variable(obj_actual);

		if(variable == 'b' || variable == 'c' || variable == 'd')
		{
			multiset1[pos] = Object1('\0', 0, 0, 0);
		}
               
                
        }
		
        multiset1[0] = Object1('f',N,0,0);
        multiset1[1] = Object1('g',N,0,0);

	// Imprimimos multiset1 para ver si se esta actualizando correctamente.

        // Vamos a ver que posee las membranas2 antes del intercambio.
        //printm2(multiset2,membranas_actuales2,T,N);
        
	for(int ind=0; ind<membranas_actuales2; ind++)
	{
		
                for(int pos=2*N;pos<(2*N)+3;pos++)
		{
                        obj_actual = multiset2[(ind*TM2)+pos];
			variable = get_variable(obj_actual);

			if(variable == 'f' || variable == 'g')
			{
				multiset2[(ind*TM2)+pos] = '0';
			}
                        
                        
		}
                
		
                multiset2[(ind*TM2)+(2*N)] = Object1('b', 0, N+1, 0);
		multiset2[(ind*TM2)+(2*N)+1] = Object1('c', 0, N+1, 0);
		multiset2[(ind*TM2)+(2*N)+2] = Object3('d', 0, N+1);
	}
        
        // Comprobamos si se ha hecho correctamente el intercambio.
        //printm2(multiset2,membranas_actuales2,T,N);
}


void generation (int N,int M,int T,Object* multiset1,Object* multiset2)
{

	unsigned int membrana_actual = 1;
	unsigned int e=1;
        unsigned int membranas_actuales2 = 0;
        int mult_a = N;
        // cout << "multiplicidad_a: " << mult_a << endl;
	
        while(e <= N)
	{
		membranas_actuales2 = division(membrana_actual,N,M,T,multiset2);
                //cout << "membranas actuales: " << membranas_actuales2 << endl;
                //printm2(multiset2,membranas_actuales2,T,N);
                
                membrana_actual = membranas_actuales2;
                // cout << "membranas actual: " << membrana_actual << endl;
                
		evolution_m2(membrana_actual,T,N,multiset2,mult_a,e);
                mult_a = mult_a -1;
                //cout << "multiplicidad de a:" << mult_a << endl; 

                e = e+1;
                
	}
        
        // Comprobamos las membranas Etiquetadas con 2
        // printm2(multiset2,membranas_actuales2,T,N);
	
        evolution_m1(N,M,T,multiset1);
        // Comprobamos la membrana Etiquetada con 1
        // printm1(multiset1);
       
        intercambio(multiset1,multiset2,membranas_actuales2,T,N,M);
        

}


void synchronization(Object* multiset2, unsigned int NumMemb,int N,int M, int T)
{
	Object obj_actual;
	char variable;
        int memb;
        int numveces;
        unsigned short mult_t;
        unsigned short mult_f;
	unsigned short i;
	unsigned short j;
        int TM2 = ((2*N)+3+1+T);
        int inicio = N+1;
        int fin = (inicio +((2*N)+M)-1);
        
    
       // Pasar a dos dimensiones
	for (int memb=0;memb<NumMemb;memb++)
	{

              for(int pos = 0;pos < N;pos++)
              {
                        obj_actual = multiset2[(memb*TM2)+pos];
                        variable = get_variable(obj_actual);
			if(variable == 'T')
			{
				i = get_i(obj_actual);
				j = get_j(obj_actual);

				if(j==0)
				{
                                        multiset2[(memb*TM2)+ pos] = Object1('T',0,i,1);
				}
                        }
              
                    
                        if(variable == 'F')
			{
				i = get_i(obj_actual);
				j = get_j(obj_actual);

				if(j==0)
				{
					multiset2[(memb*TM2)+ pos] = Object1('F',0,i,1);
				}
                        }
              }     
        }

        // Comprobacion que se ha pasado a dos dimensiones
        //printm2(multiset2,NumMemb,T,N);
	
        // Evolucionamos una vez obtenidas F y T en dos dimensiones.
        for (int memb=0;memb<NumMemb;memb++)
	{
        
            for (int numveces = 0; numveces<M;numveces++)
            {
                  for(int pos = 0;pos < N;pos++)
                  {
                        obj_actual = multiset2[(memb*TM2)+pos];
                        variable = get_variable(obj_actual);
                        j = get_j(obj_actual);
			
                        if(variable == 'T' && j<=M)
			{
                                i = get_i(obj_actual);
                                
				multiset2[(memb*TM2)+ pos] = Object1('T',0,i,j+1);
                                multiset2[(memb*TM2)+ (pos+N)] = Object1('t',j,i,0);
                        }
              
                    
                        else if(variable == 'F' && j<=M)
			{
		
                                i = get_i(obj_actual);
                               
				multiset2[(memb*TM2)+ pos] = Object1('F',0,i,j+1);
				multiset2[(memb*TM2)+ (pos+N)] = Object1('f',j,i,0);
			}

		  }
               
              }
            
          }
        
        //printm2(multiset2,NumMemb,T,N);
       
        //EVOLUCIONAR las b,c y d
        
          for (int memb=0;memb<NumMemb;memb++)
          {
                for(int cont=inicio;cont<= fin;cont++)
                {
                        multiset2[(memb*TM2)+ (2*N)] = Object1('b',0,cont+1,0);
                        multiset2[(memb*TM2)+ (2*N)+1] = Object1('c',0,N+1,0);
                        multiset2[(memb*TM2)+ (2*N)+2] = Object3('d',0,cont+1);
                }
          }
        
        //printm2(multiset2,NumMemb,T,N);
}


void checking_step(unsigned int NumMemb, unsigned short N, unsigned short M, unsigned short T, Object* multiset2)
{
        Object obj_actual;
        Object obj_aux;
        int memb;
        int contador;
        int pos = 2*N +3 +1;
        int TM2 = pos + T; 
        char variable1;
        char variable2;
        short int mult;
        short int j;
        short int i;
        int inicio = (3*N)+M+1;
        int fin = inicio + (N*M)-1;
        
        
        for(memb=0;memb<NumMemb;memb++)
	{
                for(int cont = pos;cont < TM2;cont++)
		{

                    obj_actual = multiset2[(memb*TM2) + cont];
                    variable1 = get_variable(obj_actual);
                    i = get_i(obj_actual);
                    j = get_j(obj_actual);
                    
                    obj_aux = multiset2[(memb*TM2)+(N-1)+i];
                    variable2 = get_variable(obj_aux);
                    mult = get_mult(obj_aux);
                    
                    if((variable1 == 'x') && (variable2 == 't'))
                    {
                        
                        
                        multiset2[(memb*TM2)+(N-1)+i] = Object1('t',mult-1,i,0); // Actualizamos t
                        multiset2[(memb*TM2)+ cont] = Object1('r',0,j,0);
                    }
                    
                    else if((variable1 == 'y') && (variable2 == 'f'))
                    {
                        
                        
                        multiset2[(memb*TM2)+(N-1)+i] = Object1('f',mult-1,i,0); // Actualizamos f
                        multiset2[(memb*TM2)+ cont] = Object1('r',0,j,0);
                    }
                    
                    
               }
                
                
        }
        
        
        //printm2(multiset2,NumMemb,T,N);
        
        // NOTA
        // Tal vez, podriamos hacer una funcion "reajustar" que nos ponga mas claro, las r 
        // en las membranas2 que poseemos, pero creo que no es necesario.
        
        // Evolucionamos las d  
      for (int memb=0;memb<NumMemb;memb++)
          {
                for(contador=inicio;contador<= fin;contador++)
                {
                    multiset2[(memb*TM2)+(2*N)+2] = Object3('d',0,contador+1);
                }
          }     
          
        cout << "contador vale " << contador << endl;
}



bool checkout_step(long int NumMemb,short int N,short int M,short int T, Object* multiset2)
{
    Object obj_actual;
    Object obj_d;
    int pos = 2*N +3 +1;
    int TM2 = pos + T; 
    int cont;
    char variable;
    short int ir;
    short int id;
    int indice = 3*N + N*M + M ;
    int fin = 3*N + N*M + 2*M + 1;
    int id_max;
    bool result =0; 
    
   
    //Obtengo el id de la primera membrana (supongo que este va a ser el mayor, si no, lo voy modificando)
    obj_d = multiset2[(0*TM2)+(2*N)+2]; 
    id_max = get_extra_i(obj_d);
    
    for(int memb=0;memb<NumMemb;memb++)
    {
        for(cont = pos;cont<TM2;cont++)
        {   
           
           obj_actual = multiset2[(memb*TM2) + cont];
           variable = get_variable(obj_actual);
           
           if(variable == 'r')
           {
               
               ir = get_i(obj_actual);
               
               if(ir <= M)
               {
                   obj_d = multiset2[(memb*TM2)+(2*N)+2];
                   id = get_extra_i(obj_d);
                     
                   if(id == indice + ir)
                   {
                        multiset2[(memb*TM2)+(2*N)+2] = Object3('d', 0, id+1);
                        
                   }   
                    
               }
              
            }
         }
        
        id = get_extra_i (multiset2[(memb*TM2)+(2*N)+2]);
        if (id > id_max)
        {
            id_max = id;
        }
        
    }   
    
     //cout << "id_max es:" << id_max << endl;
     //printm2(multiset2,NumMemb,T,N);
        
   
     // Comprobacion final
     if(id_max == fin)
     {
            result = true;
     }
        
   
    return result;


}

Object* inicializa_multiset1(Object* multiset1)
{
    
        multiset1[0] = Object1('b', 0, 1, 0);
        multiset1[1] = Object1('c', 0, 1, 0);
        multiset1[2] = Object3('d', 0, 1);
        multiset1[3] = Object1('s', 0,1,0);
        multiset1[4] = Object1('n', 0,1,0);
        
        return multiset1;
        
}

Object* inicializa_multiset2(Object* multiset2, int TM2, short int N,Object* cnf, short int T)
{
    
    int iniciocod = (2*N)+3+1;
    int i;
    
    
    for(int pos=0;pos<2*N;pos++)
    {
        multiset2[pos] = Object1('\0', 0, 0, 0);;
    }
    
    multiset2[2*N] = Object1('f', 0, 1, 0);
    multiset2[(2*N)+1] = Object1('g', 0, 1, 0);
    multiset2[(2*N)+2] = 0;
   
   
    multiset2[(2*N)+3] = Object2('a',N);  
        
    
    
    for(int ind = 0;ind < T;ind++)
    {
        multiset2[iniciocod + ind] = cnf[ind];   
    }
    
    return multiset2;
}

bool Sequential_hybrid_solver(int N, int M, int T, Object * cnf) {

	long int NumMemb=(long int)pow(2,N);
	bool solution;
        struct timeval tini, tfin;
	long time;
        int TM1 = 5;
        int TM2 = (2*N)+3+1+T;

        Object* multiset1 = new Object[TM1];
	Object* multiset2 = new Object[TM2*NumMemb];
        
        multiset1 = inicializa_multiset1(multiset1);
        multiset2 = inicializa_multiset2(multiset2,TM2,N,cnf,T);

        // Imprime el cnf inicial.
        printcnf(cnf,T);
        
        // Imprime las membranas con la que estamos trabajando inicialmente 
        //printm1(multiset1);
        //printm2(multiset2,1,T,N);
	

	gettimeofday(&tini, NULL);
	
       
	generation(N,M,T,multiset1,multiset2);
        
        // Comprobamos fase de Generacion
        //printm1(multiset1);
        //printm2(multiset2,NumMemb,T,N);
         
        

	synchronization(multiset2,NumMemb,N,M,T);
        
        //cout << "Comprobamos fase de Sincronizacion" << endl;
        //printm1(multiset1);
        //printm2(multiset2,NumMemb,T,N);
        
         

        checking_step(NumMemb,N,M,T,multiset2);
        
        //cout << "Comprobamos fase de Checking" << endl;
        //printm1(multiset1);
        //printm2(multiset2,NumMemb,T,N);
        
        
       
	solution = checkout_step(NumMemb,N,M,T,multiset2);
        
        /*cout << "Comprobamos fase de Checkout" << endl;
        printm1(multiset1);
        printm2(multiset2,NumMemb,T,N);
         */
        
        gettimeofday(&tfin, NULL);

        time= (tfin.tv_sec - tini.tv_sec)*1000000 + tfin.tv_usec - tini.tv_usec;
        
        cout << endl << "Execution time: " << time/1000.0 << " ms" << endl;

	delete multiset1;
        delete multiset2;

	return solution;
}







