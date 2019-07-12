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

#include <cutil_inline.h>
#include <iostream>
#include <math.h>
#include "Object.h"
#include "GPU_solver.h"

#include "gpu_synchronization_kernel_a.cu"
#include "gpu_synchronization_kernel_b.cu"
#include "gpu_synchronization_kernel_c.cu"
#include "gpu_synchronization_kernel_d.cu"
#include "gpu_checking_kernel.cu"
#include "gpu_checkout_kernel.cu"
#include "gpu_division_evolution_kernel.cu" 

#define MAX_BLOCKS_X 32768
using namespace std;

void gpu_printm1(Object* device_multiset) {
    
    Object * multiset;
    int TM1 = 5;
    multiset = new Object[TM1];
    
    cutilSafeCall(cudaMemcpy(multiset, device_multiset, sizeof(Object)*TM1 , cudaMemcpyDeviceToHost));
    
    cout << "Number of membranes: " << 1 << endl;

            cout << get_variable(multiset[0]) << "_" << get_i(multiset[0]) << "," << get_j(multiset[0]) << "*" << get_mult(multiset[0]) << "; ";
            cout << get_variable(multiset[1]) << "_" << get_i(multiset[1]) << "," << get_j(multiset[1]) << "*" << get_mult(multiset[1]) << "; ";   
            cout << get_variable(multiset[2]) << "_" << get_extra_i(multiset[2]) << "*" << get_mult(multiset[2])<< ";";
            cout << get_variable(multiset[3]) << "_" << get_i(multiset[3]) << "," << get_j(multiset[3]) << "*" << get_mult(multiset[3]) << "; ";
            cout << get_variable(multiset[4]) << "_" << get_i(multiset[4]) << "," << get_j(multiset[4]) << "*" << get_mult(multiset[4]) << "; ";
            
            
        cout << endl;
}

void gpu_printm2(Object* device_multiset, unsigned int number_membranes, unsigned int T, unsigned int N) {
    
        int TM2 = 2*N + 3 + 1+ T;
        Object * multiset;
	multiset = new Object[number_membranes*TM2];
        
        cutilSafeCall(cudaMemcpy(multiset, device_multiset, sizeof(Object) * TM2 * number_membranes, cudaMemcpyDeviceToHost));
        
	cout << "Number of membranes: " << number_membranes << endl;

        cout << "Multisets: ";
	for (int i=0; i<number_membranes; i++) {
		cout << "|"<< i << "|: ";

                
                for(int c=0;c<(2*N)+2;c++)
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

Object* gpu_inicializa_multiset1(Object* multiset1)
{
    
        multiset1[0] = Object1('b', 0, 1, 0);
        multiset1[1] = Object1('c', 0, 1, 0);
        multiset1[2] = Object3('d', 0, 1);
        multiset1[3] = Object1('s', 0,1,0);
        multiset1[4] = Object1('n', 0,1,0);
        
        return multiset1;
        
}

Object* gpu_inicializa_multiset2(Object* multiset2, int TM2, short int N,Object* cnf, short int T)
{
    
    int iniciocod = (2*N)+3+1;
    
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

void gpu_evolution_memb1(short int N,short int M,unsigned short T,Object* multiset1, int e)
{
	Object obj_actual;
	char variable;
        int TM1 = 5;

	for(int pos= 0;pos<TM1;pos++)
	{
                obj_actual = multiset1[pos];
		variable = get_variable(obj_actual);
		if(variable == 'b' || variable == 'c')
		{
			
			multiset1[pos] = Object1(variable,e,e+1 ,0);
                }
                if (variable == 'd')
                {
                    multiset1[pos] = Object3(variable,e,e+1);
                }

	}
}

void gpu_intercambio_memb1 (Object* multiset1, int N)
{
	Object obj_actual;
        char variable = '\0';
        int TM1 = 5;
        
        // Actualizamos el multiset1
	for(int pos= 0;pos<TM1;pos++)
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
}

void ordenacion_seleccion_gpu(unsigned int NumMemb, unsigned short N, unsigned short M, unsigned short T, Object* multiset2)
{
    Object obj_actual;
    Object obj_aux;
    char variable;
    int memb;
    int indice;
    short int i;
    int pos = (2*N)+3+1;
    int TM2 = pos + T;
    int final;
   
    
    
        for(memb=0;memb<NumMemb;memb++)
        {
            final = TM2-1;
            
            for(indice = pos; indice < TM2;indice++)
            {
                obj_actual = multiset2[(memb*TM2)+indice];
                variable = get_variable(obj_actual);
                
                if(variable == 'r')
                {
                    i = get_i(obj_actual);
                    obj_aux = multiset2[(memb*TM2)+(pos-1)+i];
                    multiset2[(memb*TM2)+(pos-1)+i] = obj_actual;
                    multiset2[(memb*TM2)+indice] = obj_aux;
                    
                }
                
             }
            
            
          }
        
 }

/*extern "C"*/ bool GPU_solver(int N, int M, int T, Object * cnf) {

	/* Para el simulador */
	
        Object* d_multiset1;
        Object* d_multiset2;
        Object* d_cnf;
        uint device;
    	cudaDeviceProp deviceProp;
        uint NumMemb = 1;
	uint MaxMemb = (uint) pow(2.0, N);
        uint TM1 = 5;
        uint TM2 = 2*N + 3 + 1 + T;
        uint membporbloquehilos;
        uint numhilos;
        uint bloquehilos;
   	dim3 grid;
    	uint blocksPerRow;
        uint rowsPerGrid;
        bool solution = false;
        bool * d_solution;
        
        // Especifico para el sistema P
         
        // Fase 1: Seleccionar la GPU 
        
        // Selecciona una de las tarjetas graficas que hay en el servidor. cudasetdevice asigna la GPU que le digas a la aplicacion
 
        char * def_dev = getenv("DEFAULT_DEVICE");
        unsigned int dev;
        if (def_dev!=NULL)
                cudaSetDevice(dev= atoi(def_dev));
        else
            // Esto último, devuelve el dispositivo con mas gigaflops que hay en el servidor(en nuestro caso, son iguales)
                cudaSetDevice(dev = cutGetMaxGflopsDeviceId());
      
        // Una vez seleccionada la gpu(tarjeta grafica) obtenemos sus propiedades 
        cutilSafeCall(cudaGetDeviceProperties(&deviceProp, device)); 
        
        // Fase 2: Conseguir las caracteristicas de la GPU para comprobar que son suficientes para poder simular el sistema P
        uint timertotal = 0;
        float TotalTime = 0.0f;

        
        uint MaxDeviceMemb = deviceProp.maxGridSize[0] * deviceProp.maxGridSize[1];
        uint sizeMemb1 = TM1*sizeof(Object);
    	uint deviceGlobalMemb = MaxMemb * TM2 * sizeof(Object);

        cutilCondition(MaxMemb <= MaxDeviceMemb);
        cutilCondition(deviceGlobalMemb <= deviceProp.totalGlobalMem);
        
        // Fase 3: Reservar el maximo de memoria para la GPU
        cutilSafeCall(cudaMalloc((void**)&d_multiset1,sizeMemb1));
        cutilSafeCall(cudaMalloc((void**)&d_solution,sizeof(bool)));
        cutilSafeCall(cudaMalloc((void**)&d_multiset2, deviceGlobalMemb));
        cutilSafeCall(cudaMalloc((void**)&d_cnf, T*sizeof(Object)));

	cutilCheckError(cutCreateTimer(&timertotal));
        cutilCheckError(cutStartTimer(timertotal));
        
        // Inicializo los multisets
        
        Object* multiset1 = new Object[TM1];
	Object* multiset2 = new Object[TM2*NumMemb];
        
        multiset1 = gpu_inicializa_multiset1(multiset1);
        multiset2 = gpu_inicializa_multiset2(multiset2,TM2,N,cnf,T);
        
        
        
        // Fase 4: Creacion y comienzo de tiempos
        
        uint timer = 0;
        float tdiv=0.0f, tsyn=0.0f, tcheck=0.0f, tout=0.0f;
    	cutilCheckError(cutCreateTimer(&timer));
        
        // Fase 5: Copiar datos iniciales
        
        cutilSafeCall(cudaMemcpy(d_cnf, cnf, sizeof(Object)*T, cudaMemcpyHostToDevice));
        cutilSafeCall(cudaMemcpy(d_solution, &solution, sizeof(bool), cudaMemcpyHostToDevice));
        cutilSafeCall(cudaMemcpy(d_multiset1, multiset1, sizeof(Object)*TM1, cudaMemcpyHostToDevice));
        cutilSafeCall(cudaMemcpy(d_multiset2, multiset2, sizeof(Object)*TM2, cudaMemcpyHostToDevice));
        
        //gpu_printm1(d_multiset1);
        //gpu_printm2(d_multiset2, NumMemb,T,N);
        
        // Fase 6: LLamada a los kernel
        
        grid = dim3(1);
        
    // FASE DE GENERACION
        
        int e = 1;
        int mult_a = N;
        int indice;
        
        
        while(e <= N)
        {
            
            mult_a = mult_a -1;
            indice = N - (mult_a);
            
            cutilCheckError(cutResetTimer(timer));
            cutilCheckError(cutStartTimer(timer));
            gpu_division_evolution_kernel<<<grid,TM2>>>(NumMemb,d_multiset2,e, mult_a, N, indice);
            cutilCheckMsg("Kernel execution failed");
            cudaDeviceSynchronize();
            cutilCheckError(cutStopTimer(timer));
          
            tdiv+= cutGetTimerValue(timer);
                    
            NumMemb <<= 1; // Esto lo que hace es multiplicar por 2, y la operacion es menos costosa (Operaciones Atómicas)
            
            // Configuracion de Parametros
            // Si hay mas de 65536 memb, usamos la segunda dimension
        	if (NumMemb <= MAX_BLOCKS_X) 
                {
            		// We can use a 1D Grid
            		blocksPerRow = NumMemb;
            		rowsPerGrid  = 1;
        	} 
                
                else {
			// We need to use 2D Grid
            		//blocksPerRow = MAX_BLOCKS_X;
            		//rowsPerGrid = numMemb/MAX_BLOCKS_X;
			blocksPerRow = rowsPerGrid = (uint) sqrt(NumMemb);

            		while ((blocksPerRow * rowsPerGrid) < NumMemb)
                		blocksPerRow++;
        	}
            
                grid = dim3(blocksPerRow, rowsPerGrid);
                //cout << "NumMembranas=" << NumMemb << ", maxgridx=" << deviceProp.maxGridSize[0] << ", blocksx=" << blocksPerRow << ", y="<<rowsPerGrid<<endl;
               
                gpu_evolution_memb1(N,M,T,multiset1,e);
              
                //gpu_evolution_m2<<<grid,TM2>>>(NumMemb,T,N,multiset2,mult_a,e);
                
                e = e+1;
            
        }
        //cutilSafeCall(cudaMemcpy(d_multiset1, multiset1, sizeof(Object)*TM1, cudaMemcpyHostToDevice));
        //gpu_printm1(d_multiset1);
        
        // Hacemos la fase de intercambio
        gpu_intercambio_memb1(multiset1, N);
        cutilSafeCall(cudaMemcpy(d_multiset1, multiset1, sizeof(Object)*TM1, cudaMemcpyHostToDevice));
        
        //cout << "Estados de las membranas al acabar la fase de Generacion" << endl;
        //gpu_printm1(d_multiset1);
        //gpu_printm2(d_multiset2, NumMemb,T,N);
        
      
        // FASE DE SINCRONIZACION
        
        //cout <<"Fase de Sincronizacion" << endl;
        /*cout <<"Opcion A1" << endl;
        numhilos = TM2;
        cutilCheckError(cutStartTimer(timer)); 
        gpu_synchronization_kernel_a<<<grid,numhilos>>>(NumMemb,d_multiset2,N,M,TM2);
        cutilCheckMsg("Kernel execution failed");
        cudaDeviceSynchronize();
        cutilCheckError(cutStopTimer(timer));
          
        tsyn= cutGetTimerValue(timer);
              
        //cout << "Estados de las membranas al acabar la fase de Sincronizacion" << endl;
        //gpu_printm1(d_multiset1);
        //gpu_printm2(d_multiset2, NumMemb,T,N);
        */
        
        /*cout <<"Opcion A2" << endl;
        numhilos = N;
        cutilCheckError(cutStartTimer(timer)); 
        gpu_synchronization_kernel_a<<<grid,numhilos>>>(NumMemb,d_multiset2,N,M,TM2);
        cutilCheckMsg("Kernel execution failed");
        cudaDeviceSynchronize();
        cutilCheckError(cutStopTimer(timer));
          
        tsyn= cutGetTimerValue(timer);
              
        //cout << "Estados de las membranas al acabar la fase de Sincronizacion" << endl;
        //gpu_printm1(d_multiset1);
        //gpu_printm2(d_multiset2, NumMemb,T,N);
        */
        
        /*cout <<"Opcion B" << endl;
        numhilos = N;
        cutilCheckError(cutResetTimer(timer));
        cutilCheckError(cutStartTimer(timer)); 
        gpu_synchronization_kernel_b<<<grid,numhilos,((2*N)+3)*sizeof(uint)>>>(NumMemb,d_multiset2,N,M,TM2);
        cutilCheckMsg("Kernel execution failed");
        cudaDeviceSynchronize();
        cutilCheckError(cutStopTimer(timer));
          
        tsyn= cutGetTimerValue(timer);
              
        //cout << "Estados de las membranas al acabar la fase de Sincronizacion" << endl;
        //gpu_printm1(d_multiset1);
        //gpu_printm2(d_multiset2, NumMemb,T,N);
        */
        
        /*cout <<"Opcion C" << endl;
        membporbloquehilos = 256/N; //
        numhilos = membporbloquehilos*N;
        bloquehilos = NumMemb/membporbloquehilos+1;
        
        //cout << NumMemb << " " << membporbloquehilos<< " " << bloquehilos << " " << blocksPerRow << " " << rowsPerGrid << endl;
        
        //grid = dim3(1);
        
        // Configuracion de Parametros
        // Si hay mas de 65536 memb, usamos la seguna dimension
        	if (bloquehilos <= MAX_BLOCKS_X) 
                {
            		// We can use a 1D Grid
            		blocksPerRow = bloquehilos;
            		rowsPerGrid  = 1;
        	} 
                
                else {
			// We need to use 2D Grid
            		//blocksPerRow = MAX_BLOCKS_X;
            		//rowsPerGrid = numMemb/MAX_BLOCKS_X;
			blocksPerRow = rowsPerGrid = (uint) sqrt(bloquehilos);

            		while ((blocksPerRow * rowsPerGrid) < bloquehilos)
                		blocksPerRow++;
        	}
            
                grid = dim3(blocksPerRow, rowsPerGrid);
        
        //cout << NumMemb << " " << membporbloquehilos<< " " << bloquehilos << " " << blocksPerRow << " " << rowsPerGrid << endl;
        cutilCheckError(cutResetTimer(timer));
        cutilCheckError(cutStartTimer(timer));        
        gpu_synchronization_kernel_c<<<grid,numhilos>>>(membporbloquehilos,NumMemb,d_multiset2,N,M,TM2);
        cutilCheckMsg("Kernel execution failed");
        cudaDeviceSynchronize();
        cutilCheckError(cutStopTimer(timer));
          
        tsyn= cutGetTimerValue(timer);
        
        //cout << "Estados de las membranas al acabar la fase de Sincronizacion" << endl;
        //gpu_printm1(d_multiset1);
        //gpu_printm2(d_multiset2, NumMemb,T,N);
       */

        
        //cout <<"Opcion D" << endl;
        membporbloquehilos = 256/N; //
        numhilos = membporbloquehilos*N;
        bloquehilos = NumMemb/membporbloquehilos+1;
        
        //cout << NumMemb << " " << membporbloquehilos<< " " << bloquehilos << " " << blocksPerRow << " " << rowsPerGrid << endl;
        
        //grid = dim3(1);
        
        // Configuracion de Parametros
        // Si hay mas de 65536 memb, usamos la seguna dimension
        	if (bloquehilos <= MAX_BLOCKS_X) 
                {
            		// We can use a 1D Grid
            		blocksPerRow = bloquehilos;
            		rowsPerGrid  = 1;
        	} 
                
                else {
			// We need to use 2D Grid
            		//blocksPerRow = MAX_BLOCKS_X;
            		//rowsPerGrid = numMemb/MAX_BLOCKS_X;
			blocksPerRow = rowsPerGrid = (uint) sqrt(bloquehilos);

            		while ((blocksPerRow * rowsPerGrid) < bloquehilos)
                		blocksPerRow++;
        	}
            
                dim3 grid_sync = dim3(blocksPerRow, rowsPerGrid);
        
        //cout << NumMemb << " " << membporbloquehilos<< " " << bloquehilos << " " << blocksPerRow << " " << rowsPerGrid << endl;
        cutilCheckError(cutResetTimer(timer));
        cutilCheckError(cutStartTimer(timer));        
        gpu_synchronization_kernel_d<<<grid_sync,numhilos,((2*N)+3)*sizeof(Object)*membporbloquehilos>>>(membporbloquehilos,NumMemb,d_multiset2,N,M,TM2);
        cutilCheckMsg("Kernel execution failed");
        cudaDeviceSynchronize();
        cutilCheckError(cutStopTimer(timer));
          
        tsyn= cutGetTimerValue(timer);
        
        //cout << "Estados de las membranas al acabar la fase de Sincronizacion" << endl;
        //gpu_printm1(d_multiset1);
        //gpu_printm2(d_multiset2, NumMemb,T,N);
       
        
        // FASE DE CHEQUEO      
        
        //cout <<"Fase de Chequeo" << endl;
        
        numhilos = T;
        cutilCheckError(cutResetTimer(timer));
        cutilCheckError(cutStartTimer(timer)); 
        gpu_checking_kernel<<<grid,numhilos,(N+T)*sizeof(Object)>>>(NumMemb,d_multiset2,N,M,TM2,T);
        cutilCheckMsg("Kernel execution failed");
        cudaDeviceSynchronize();
        cutilCheckError(cutStopTimer(timer));
          
        tcheck= cutGetTimerValue(timer);
              
        //cout << "Estados de las membranas al acabar la fase de Chequeo" << endl;
        //gpu_printm1(d_multiset1);
        //gpu_printm2(d_multiset2, NumMemb,T,N);
     
        
        // FASE DE SALIDA      
        
        //cout <<"Fase de Salida" << endl;
        
        numhilos = T;
        cutilCheckError(cutResetTimer(timer));
        cutilCheckError(cutStartTimer(timer)); 
        gpu_checkout_kernel<<<grid,numhilos,(T)*sizeof(Object)>>>(NumMemb,d_multiset2,N,M,TM2,T,d_solution);
        cutilCheckMsg("Kernel execution failed");
        cudaDeviceSynchronize();
        cutilCheckError(cutStopTimer(timer));
          
        tout= cutGetTimerValue(timer);
              
        //cout << "d_solution vale:" << d_solution << endl;
        //cout << "Estados de las membranas al acabar la fase de Salida" << endl;
        cutilSafeCall(cudaMemcpy(&solution, d_solution, sizeof(bool), cudaMemcpyDeviceToHost));
        //gpu_printm1(d_multiset1);
        //gpu_printm2(d_multiset2, NumMemb,T,N);
        
        cutilCheckError(cutStopTimer(timertotal));
        TotalTime = cutGetTimerValue(timertotal);
        
        //cout << "Tiempo: division=" << tdiv << "ms, sincronización=" << tsyn << "ms, chequeo=" << tcheck << "ms, salida="<< tout <<"ms" << endl;
        //cout << "sincronización=" << tsyn << endl;
        cout << "Execution time: " << TotalTime << " ms" << endl;
        cutilCheckError(cutDeleteTimer(timer));
        cutilCheckError(cutDeleteTimer(timertotal));
        
        return solution;
        
}       
