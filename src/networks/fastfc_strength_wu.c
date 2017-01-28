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
#include <omp.h>

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    int node_i,node_j;
    double *W,*S,s_aux;
    const int n_nodes=(int)mxGetM(prhs[0]);
    
    plhs[0]=mxCreateDoubleMatrix(1,n_nodes,mxREAL);
    W=(double*)mxGetPr(prhs[0]);
    S=(double*)mxGetPr(plhs[0]);
    
    #pragma omp parallel
    {
        #pragma omp for private(node_i,node_j) reduction(+:s_aux)
        for(node_i=0;node_i<n_nodes;node_i++)
        {
            s_aux=0;
            for(node_j=0;node_j<n_nodes;node_j++)
            {
                s_aux+=W[node_i*n_nodes+node_j];
            }
            S[node_i]=s_aux;
        }
    }
}