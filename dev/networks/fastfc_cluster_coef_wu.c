//
//    Copyright (c) 2014-2015 Juan Garcia-Prieto Cuesta
//
//    This file is part of Fast Functional Connectivity (FastFC)
//
//    FastFC is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    FastFC is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with FastFC.  If not, see <http://www.gnu.org/licenses/>.
//
//    ------------------------------------------ 
//    Contact:
//    Juan Garcia-Prieto <juangpc@gmail.com>
//    ------------------------------------------
//
//    Please consider helping by citing our research.
// 
//    J. Garcia-Prieto, E. Pereda
//
#include <mex.h>
#include <math.h>
#include <omp.h>

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    int *K,node_i,node_j,z;
    double *W,*W2,*C,cyc3,w2_aux;
    
    const int n_nodes=(int)mxGetM(prhs[0]);
    W=(double*)mxGetPr(prhs[0]);
    K=(int*)mxMalloc(n_nodes*sizeof(int));
    W2=(double*)mxCalloc(n_nodes*n_nodes,sizeof(double));
    plhs[0]=mxCreateDoubleMatrix(n_nodes,1,mxREAL);
    C=(double*)mxGetPr(plhs[0]);
    
    omp_set_num_threads(omp_get_num_procs());
    #pragma omp parallel
    {
        #pragma omp for private(node_i,node_j,z) reduction(+:w2_aux)
        for(node_i=0;node_i<n_nodes;node_i++)
        {
            K[node_i]=n_nodes;
            for(node_j=0;node_j<node_i;node_j++)
                K[node_i]-=(int)(W[node_i*n_nodes+node_j]==0.);
            for(;node_j<n_nodes;node_j++)
            {
                K[node_i]-=(int)(W[node_i*n_nodes+node_j]==0.);
                w2_aux=0.;
                for(z=0;z<n_nodes;z++)
                    w2_aux+=pow(W[node_i*n_nodes+z],1./3.)*
                            pow(W[node_j*n_nodes+z],1./3.);
                W2[node_i*n_nodes+node_j]=w2_aux;
                W2[node_j*n_nodes+node_i]=w2_aux;
            }
        }
    }   
    
    #pragma omp parallel
    {
        #pragma omp for private(node_i,z) reduction(+:cyc3)
        for(node_i=0;node_i<n_nodes;node_i++)
        {
            cyc3=0.;
            for(z=0;z<n_nodes;z++)
                cyc3+=W2[node_i*n_nodes+z]*pow(W[node_i*n_nodes+z],1./3.);
            if(cyc3==0.)
                C[node_i]=0.;
            else
                C[node_i]=cyc3/(K[node_i]*(K[node_i]-1));
        }
    }

    mxFree(K);
    mxFree(W2);
}