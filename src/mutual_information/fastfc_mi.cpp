/*
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

// Copyright 2009 Alexander Kraskov, Harald Stoegbauer, Peter Grassberger
//-----------------------------------------------------------------------------------------
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should receive a copy of the GNU General Public License
// along with this program.  See also <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------------------------- 
// Contacts:
//
// Harald Stoegbauer <h.stoegbauer@gmail.com>
// Alexander Kraskov <alexander.kraskov@gmail.com>
//-----------------------------------------------------------------------------------------
// Please reference
// 
// A. Kraskov, H. Stogbauer, and P. Grassberger,
// Estimating mutual information.
// Phys. Rev. E 69 (6) 066138, 2004
//
// in your published research.
*/
#define _CRT_SECURE_NO_WARNINGS

#include <mex.h>
#include <matrix.h>
#include <omp.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "miutils.h"

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    int N,edim,tau,K,n_sensors,sensor_i,sensor_j,i,d,n_threads;
    const int dim=2;
    double addnoise=-1;
    double *data,**x;//*min,*max;
    double *psi;
//     double s,me;
    double *out_mi,*mi,mi_value;
    
    
    data=(double*)mxGetPr(prhs[0]);
    N=(int)mxGetM(prhs[0]);
    n_sensors=(int)mxGetN(prhs[0]);
	edim=(int)mxGetScalar(prhs[1]);
    tau=(int)mxGetScalar(prhs[2]);
    K=(int)mxGetScalar(prhs[3]);
    const int n_indexes=n_sensors*n_sensors;
    n_threads=omp_get_num_procs();
            
            
    if (nrhs<4)
        mexErrMsgIdAndTxt("MATLAB:trial_mi_milca","Not enough input arguments\n");
    if (nrhs==6)
        addnoise=mxGetScalar(prhs[4]);
    if (nrhs>6)
        mexErrMsgIdAndTxt("MATLAB:trial_mi_milca","Too many input arguments\n");
    
    x=(double**)mxMalloc(dim*sizeof(double*));

    for (d=0;d<dim;d++)
        x[d]=(double*)mxMalloc(N*sizeof(double));
  
    psi=(double*)mxMalloc((N+1)*sizeof(double));
    psi[1]=-(double).57721566490153;
    for (i=1;i<N;i++) 
        psi[i+1]=psi[i]+1/(double)i;
    
    mi=(double*)mxMalloc(n_indexes*sizeof(double));
    
    plhs[0]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    out_mi=(double*)mxGetPr(plhs[0]);
    
    omp_set_num_threads(n_threads);
    
    for(sensor_i=0;sensor_i<n_sensors;++sensor_i)
    {
        #pragma omp parallel 
        {
            #pragma omp for private(sensor_j,i,mi_value)
            for(sensor_j=sensor_i;sensor_j<n_sensors;++sensor_j)
            {
                for(i=0;i<N;++i)
                {
                    x[0][i]=data[i+sensor_i*N];
                    x[1][i]=data[i+sensor_j*N];
                }
                redr_embed(x,(int)dim, edim, tau,(int)N,K,psi,&mi_value);
                mi[sensor_j+(n_sensors*sensor_i)]=mi_value;
                mi[sensor_i+(n_sensors*sensor_j)]=mi_value;
            }
        }
    }
    
//     #pragma omp parallel for private(i) shared(mi,out_mi) num_threads(n_threads)
    for(i=0;i<n_indexes;++i)
        out_mi[i]=mi[i];

    
    for (d=0;d<dim;d++) 
        mxFree(x[d]); 
    
    mxFree(mi);
    mxFree(x);
    mxFree(psi);
}