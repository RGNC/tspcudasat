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

__global__ static void gpu_synchronization_kernel_c(uint membporbloquehilos,uint NumMemb, Object* multiset2, uint N, uint M, uint TM2)
{
    const uint bid = (blockIdx.y * gridDim.x) + blockIdx.x;
    const uint tid = threadIdx.x;
    const uint idMemb = (membporbloquehilos*bid)+(tid/N); // Me idedntifica la membrana dentro del bloque de hilo correspondiente
    const uint i = tid%N; // Para que todos los indices sean iguales en todas las membranas
    const uint blockSize = TM2;
    __shared__ int pos1;
    __shared__ int pos2;
    
    if (idMemb>=NumMemb) return;
    
    uint inicio = N+1;
    uint fin =inicio + (N + M)-1;
    uint pos_b = 2*N;
    uint pos_d = (2*N)+2;
   
   if(tid == 0)  //Para inicializar la variable compartida, por ejemplo seleccionamos este hilo(podríamos haber seleccionado otro cualquiera)
    {
        pos1=0;
        pos2=0;
        
    }
    
    __syncthreads();
    
   
    while(inicio <= fin)
    {
        
        
       if(pos2<pos1 && i>=pos2 && i<pos1)
        {
            if((tid < N) && (d_get_j(multiset2[(idMemb*blockSize)+i])==M+1))
            {
                pos2 = pos2 +1;
                
            }
            
            if((d_get_variable(multiset2[(idMemb*blockSize)+i]) == 'T') && (d_get_j(multiset2[(idMemb*blockSize)+i])<=M))
           {
                multiset2[(idMemb*blockSize)+i+N] = d_set_Object1('t', d_get_j(multiset2[(idMemb*blockSize)+i]), d_get_i(multiset2[(idMemb*blockSize)+i]), 0);
                multiset2[(idMemb*blockSize)+i] = d_set_Object1('T', 0,d_get_i(multiset2[(idMemb*blockSize)+i]),(d_get_j(multiset2[(idMemb*blockSize)+i]))+1); 
               
           }
        
           if((d_get_variable(multiset2[(idMemb*blockSize)+i]) == 'F') && (d_get_j(multiset2[(idMemb*blockSize)+i])<=M))
           {
                multiset2[(idMemb*blockSize)+i+N] = d_set_Object1('f', d_get_j(multiset2[(idMemb*blockSize)+i]), d_get_i(multiset2[(idMemb*blockSize)+i]), 0);
                multiset2[(idMemb*blockSize)+i] = d_set_Object1('F', 0,d_get_i(multiset2[(idMemb*blockSize)+i]),(d_get_j(multiset2[(idMemb*blockSize)+i]))+1); 
                
           }
        
        }
        
      __syncthreads();
         
         
        
        if(pos1<N && i==pos1)
        {
            if(d_get_variable(multiset2[(idMemb*blockSize)+i]) == 'T')
            {
                multiset2[(idMemb*blockSize)+i] = d_set_Object1('T', 0,d_get_i(multiset2[(idMemb*blockSize)+i]),1);
            }
        
            if(d_get_variable(multiset2[(idMemb*blockSize)+i]) == 'F')
            {
                multiset2[(idMemb*blockSize)+i] = d_set_Object1('F', 0,d_get_i(multiset2[(idMemb*blockSize)+i]),1);
            }
        
        }
        
      __syncthreads();
      
    if(tid==0)
    {
          pos1 = pos1 + 1;
        
    } 
      
    
    if(i == 0)
    {
        multiset2[(idMemb*blockSize)+pos_b] = d_set_Object1('b', 0, inicio+1, 0);
    }
    
    if(i == 0)
    {
        multiset2[(idMemb*blockSize)+pos_d] = d_set_Object3('d', 0, inicio+1);
    }
     
    
          inicio = inicio+1;
    
          
  }
   

}
