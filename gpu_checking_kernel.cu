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


__global__ static void gpu_checking_kernel(uint NumMemb, Object* multiset2, uint N, uint M, uint TM2, uint T)
{
    const uint bid = (blockIdx.y * gridDim.x) + blockIdx.x;
    const uint tid = threadIdx.x;
    
    
    __shared__ uint evol_d;
    Object obj_cnf;
    Object obj_tf;
    int c;
    uint inicio = (2*N)+M+1;
    uint fin = (inicio)+((N*M)-1);
    
     // Nos creamos el array en memoria compartida
    
    extern __shared__ uint sm2[];
    
    
    if (bid >= NumMemb) return;
    
    // Inicializamos la variable compartida
    if(tid==0)
    {
        evol_d = d_get_i(multiset2[(bid*TM2)+(2*N)+2]);
    }
    
     __syncthreads();
    
    // Inicializamos el array
    
    if(tid < N)
    {    
        sm2[tid] = multiset2[(bid*TM2)+N+tid];
    }
    
     sm2[N+tid] = multiset2[(bid*TM2)+(2*N)+4+tid];
    
    // Implementacion de la funcion
    // Esto es totalmente secuencial;solo trabaja un objeto(solo trabaja un solo hilo a la vez por el modelo(pos la b, que solo tenemos una para gastar))
    for(c=0;c<T;c++)
    {
        
        if(tid==c)
        {
            evol_d++;
            Object obj_cnf = sm2[N+tid];
            Object obj_tf = sm2[d_get_i(obj_cnf)-1];
            
            if(((d_get_variable(obj_cnf) == 'x') && (d_get_variable(obj_tf) == 't')) || ((d_get_variable(obj_cnf) == 'y') && (d_get_variable(obj_tf) == 'f')))
            {
                sm2[N+tid] = d_set_Object1('r', 0, d_get_j(obj_cnf), 0);
                sm2[d_get_i(obj_cnf)-1] = d_set_Object1(d_get_variable(obj_tf), d_get_mult(obj_tf)-1, d_get_i(obj_tf), 0);
                
            }
        }
        
        __syncthreads();
        
    }
    
    //cout << "C vale" << c << endl; 
     
    if(tid==0)
    {
        while(evol_d <= fin)
        {
            evol_d++;
        }
        
    }
    
    __syncthreads();
    
    // Reestructuracion del array
    
    if(tid < N)
    {    
        multiset2[(bid*TM2)+N+tid] = sm2[tid];
    }
    multiset2[(bid*TM2)+(2*N)+4+tid]=sm2[N+tid];
    
    if(tid==0)
    {    
        multiset2[(bid*TM2)+(2*N)+2] = d_set_Object3('d', 0, evol_d);
    }
}
    
    
    
    
    
    
