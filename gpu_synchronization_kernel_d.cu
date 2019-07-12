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

#include "Object.cuH"

__global__ static void gpu_synchronization_kernel_d(uint membporbloquehilos,uint NumMemb, Object* multiset2, uint N, uint M, uint TM2) {
    const uint bid = (blockIdx.y * gridDim.x) + blockIdx.x;
    const uint tid = threadIdx.x;
    const uint idm = (membporbloquehilos * bid)+(tid / N); // Me idedntifica la membrana dentro del bloque de hilo correspondiente
    const uint sidm = idm - membporbloquehilos*bid;
    const uint sSize = (2 * N) + 3;
    const uint i = tid % N; // Para que todos los indices sean iguales en todas las membranas

    __shared__ int pos1;
    __shared__ int pos2;

    uint inicio = N + 1;
    uint fin = inicio + (N + M) - 1;
    uint pos_b = 2 * N;
    uint pos_d = (2 * N) + 2;


    // Nos creamos el array en memoria compartida

    extern __shared__ uint sm2[]; // Array Dinámico

    if (idm >= NumMemb) return;

    // Inicializamos el array

    sm2[(sidm * sSize) + i] = multiset2[(idm * TM2) + i];
    sm2[(sidm * sSize) + N + i] = 0;

    if ((N >= 3) && (i < 3)) {
        sm2[(sidm * sSize)+(2 * N) + i] = multiset2[(idm * TM2)+((2 * N) + i)];
    }
    else if ((N < 3) && (i == 0)) {
        for (int ind = 0; ind < 3; ind++) {
            sm2[(sidm * sSize)+(2 * N) + ind] = multiset2[(idm * TM2)+((2 * N) + ind)];
        }
    }

    __syncthreads();


    // Cuerpo


    if (tid == 0) //Para inicializar la variable compartida, por ejemplo seleccionamos este hilo(podríamos haber seleccionado otro cualquiera)
    {
        pos1 = 0;
        pos2 = 0;

    }

    __syncthreads();


    while(inicio <= fin)
    {
        
        
       if(pos2<pos1 && i>=pos2 && i<pos1)
        {
            if((tid < N) && (d_get_j(sm2[(sidm*sSize)+i])==M+1))
            {
                pos2 = pos2 +1;
                
            }
            
            if((d_get_variable(sm2[(sidm*sSize)+i]) == 'T') && (d_get_j(sm2[(sidm*sSize)+i])<=M))
           {
                sm2[(sidm*sSize)+i+N] = d_set_Object1('t', d_get_j(sm2[(sidm*sSize)+i]), d_get_i(sm2[(sidm*sSize)+i]), 0);
                sm2[(sidm*sSize)+i] = d_set_Object1('T', 0,d_get_i(sm2[(sidm*sSize)+i]),(d_get_j(sm2[(sidm*sSize)+i]))+1); 
               
           }
        
           if((d_get_variable(sm2[(sidm*sSize)+i]) == 'F') && (d_get_j(sm2[(sidm*sSize)+i])<=M))
           {
                sm2[(sidm*sSize)+i+N] = d_set_Object1('f', d_get_j(sm2[(sidm*sSize)+i]), d_get_i(sm2[(sidm*sSize)+i]), 0);
                sm2[(sidm*sSize)+i] = d_set_Object1('F', 0,d_get_i(sm2[(sidm*sSize)+i]),(d_get_j(sm2[(sidm*sSize)+i]))+1); 
                
           }
        
        }
        
      __syncthreads();
         
         
        
        if(pos1<N && i==pos1)
        {
            if(d_get_variable(sm2[(sidm*sSize)+i]) == 'T')
            {
                sm2[(sidm*sSize)+i] = d_set_Object1('T', 0,d_get_i(sm2[(sidm*sSize)+i]),1);
            }
        
            if(d_get_variable(sm2[(sidm*sSize)+i]) == 'F')
            {
                sm2[(sidm*sSize)+i] = d_set_Object1('F', 0,d_get_i(sm2[(sidm*sSize)+i]),1);
            }
        
        }
        
      __syncthreads();
      
    if(tid==0)
    {
          pos1 = pos1 + 1;
        
    } 
      
    
    if(i == 0)
    {
        sm2[(sidm*sSize)+pos_b] = d_set_Object1('b', 0, inicio+1, 0);
    }
    
    if(i == 0)
    {
        sm2[(sidm*sSize)+pos_d] = d_set_Object3('d', 0, inicio+1);
    }
     
    
          inicio = inicio+1;
    
          
  }
   


    // Reestructuracion del array

    multiset2[(idm * TM2) + i] = sm2[(sidm * sSize) + i];
    multiset2[(idm * TM2)+(N + i)] = sm2[(sidm * sSize) + N + i];

    if ((N >= 3) && (i < 3)) {
        multiset2[(idm * TM2)+(2 * N) + i] = sm2[(sidm * sSize)+(2 * N) + i];

    } else if ((N < 3) && (i == 0)) {
        for (int ind = 0; ind < 3; ind++) {
            multiset2[(idm * TM2)+(2 * N) + ind] = sm2[(sidm * sSize)+(2 * N) + ind];
        }

    }

    __syncthreads();

} 
