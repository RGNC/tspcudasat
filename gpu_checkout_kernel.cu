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


__global__ static void gpu_checkout_kernel(uint NumMemb, Object* multiset2, uint N, uint M, uint TM2, uint T, bool *d_solution)
{
    const uint bid = (blockIdx.y * gridDim.x) + blockIdx.x;
    const uint tid = threadIdx.x;
    
    __shared__ int id;
    __shared__ bool enc;
    int cont;
    uint final = (2*N)+(N*M)+(2*M)+1;
    
     // Nos creamos el array en memoria compartida
    
    extern __shared__ uint sm2[];
    
    
    if (bid >= NumMemb) return;
    
    // Inicializamos la variable compartida
     if(tid==0)
    {
        id = d_get_extra_i(multiset2[(bid*TM2)+(2*N)+2]);
        enc=false;
        
    }
    
   
   __syncthreads();
    
    // Inicializamos el array
    
    // Primeramente, hacemos una copia en memoria compartida de cnf
    if(tid < T)
    {    
       sm2[tid] = multiset2[(bid*TM2)+(2*N)+4+tid];
    }
     
    
    // Implementacion de la funcion
    
    // Aqui se lanzan T hilos simultaneamente por cada membrana         
    for(cont=1;cont<=M;cont++) 
    {        
            
            if((d_get_variable(sm2[tid]) == 'r') && (d_get_i(sm2[tid])== cont))
            {
                enc=true;
                
            }
            
            __syncthreads();
            
            if (!enc)
                break;
            else if((enc) && (tid==0))
            {
                id++;
                enc=false;
            }
            
            __syncthreads();
    }
        
     if(id == final)
     {
         *d_solution = true;
     }
      __syncthreads();   
    
     // Reestructuracion del array
    if(tid < T)
    {    
       multiset2[(bid*TM2)+(2*N)+4+tid]=sm2[tid];
    }
     
     
     // En la ultima posicion, es decir en la posicion de la memoria compartida (T), almacenamos d
     if(tid==0)
     {
        multiset2[(bid*TM2)+(2*N)+2] = d_set_Object3('d',0,id);
     }
    
    
}
    
