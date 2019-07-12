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

//#include "Object.h"
#include "Object.cuH"


__global__ static void gpu_division_evolution_kernel(uint NumMemb,Object* multiset2,int e, int mult_a, int N, int indice)

{
    
    const uint bid = (blockIdx.y * gridDim.x) + blockIdx.x;
    const uint tid = threadIdx.x;
    const uint numBlocks = NumMemb;
    const uint blockSize = blockDim.x; //numero de hilos totales que hay
    const uint inicio_a = (2*N)+3;
    const uint pos_d = (2*N)+2;
    
    char var[] = {'b','c'};
    
    
    if (bid >= NumMemb)
        return;
    
    multiset2[(blockSize*(numBlocks+bid))+tid] = multiset2[bid*blockSize+tid];
    
     
    __syncthreads();
       
    if (tid == inicio_a)
    {
           multiset2[(blockSize*bid)+tid] = d_set_Object2('a',mult_a);
           multiset2[(blockSize*(numBlocks+bid))+tid] = d_set_Object2('a',mult_a);
    }
       
        __syncthreads();
       
      
       
       if(tid == e-1)
       {
           multiset2[(blockSize*(numBlocks+bid))+tid] = d_set_Object1('F', 0, indice, 0);
           multiset2[bid*blockSize+tid] = d_set_Object1('T', 0, indice, 0);
           
       }
       
        __syncthreads();
        
        if(e==N && tid<2)
        {
            
            multiset2[(bid*blockSize)+(2*N)+tid]=d_set_Object1(var[tid],0,N+1,0);
            multiset2[(blockSize*(numBlocks+bid))+(2*N)+tid] = d_set_Object1(var[tid],0,N+1,0);
        }
       
        if(e==N && tid==pos_d)
        {
             multiset2[(bid*blockSize)+tid]= d_set_Object3('d', 0, N+1);
             multiset2[(blockSize*(numBlocks+bid))+tid] = d_set_Object3('d', 0, N+1);
        
        }
       
     
}
