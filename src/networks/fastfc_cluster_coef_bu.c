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
//	 -------------------------------------------------------
//   Please consider citing our work:
//  
//   Efficient computation of functional brain networks: towards real-time functional connectivity
//	 Frontiers in Neuroinformatics (2017) García-Prieto Juan, Bajo Ricardo, Pereda Ernesto
//
//   -------------------------------------------------------
//   Contact: Juan Garcia-Prieto    juangpc (at) gmail.com
//   -------------------------------------------------------
//
#include <mex.h>
#include <omp.h>

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
 
    int *n_neis,*n_triangs,*neis,link_i,link_j,nei_i,nei_j;
    double *A,*C;
    
    const int n_links=(int)mxGetM(prhs[0]);
    
    A=(double*)mxGetPr(prhs[0]);
    
    neis=(int*)mxMalloc(n_links*(n_links-1)*sizeof(int));
    n_neis=(int*)mxCalloc(n_links,sizeof(int));
    n_triangs=(int*)mxCalloc(n_links,sizeof(int));
    
    plhs[0]=mxCreateDoubleMatrix(n_links,1,mxREAL);
    C=(double*)mxGetPr(plhs[0]);
    
    omp_set_num_threads(omp_get_num_procs());
    #pragma omp parallel
    {
        #pragma omp for private(link_i,link_j)
        for(link_i=0;link_i<n_links;link_i++)
            for(link_j=0;link_j<n_links;link_j++)
                if(((int)A[link_i*n_links+link_j]) && (link_i!=link_j))
                {
                    neis[link_i*(n_links-1)+n_neis[link_i]]=link_j;
                    n_neis[link_i]++;
                }
        
        #pragma omp for private(link_i,nei_i,nei_j)
        for(link_i=0;link_i<n_links;link_i++)
            for(nei_i=0;nei_i<n_neis[link_i]-1;nei_i++)
                for(nei_j=nei_i+1;nei_j<n_neis[link_i];nei_j++)
                    if((int)A[neis[link_i*(n_links-1)+nei_i]*n_links+neis[link_i*(n_links-1)+nei_j]])
                        n_triangs[link_i]++;

        #pragma omp for private(link_i)
        for(link_i=0;link_i<n_links;link_i++)
        {
           if (n_neis[link_i]>1)
           {
                C[link_i]=(double)n_triangs[link_i]/(n_neis[link_i]*(n_neis[link_i]-1)*.5);
            }else{
                C[link_i]=0;
            }
        }
    }
    
	mxFree(neis);
	mxFree(n_neis);
	mxFree(n_triangs);    

}
