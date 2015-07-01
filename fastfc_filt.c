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
#include <stdio.h>

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{

    // y=trial_filtfilt(data, 
    //
    // data = trial. sensors=columns
    // filt = filter denominator.
    //
    
    mwSize n_samples,n_sensors,l_x_p,l_filt,l_filt_p;
    double *x,*filt,*y;
    
    x=(double*)mxGetPr(prhs[0]);
    n_samples=(mwSize)mxGetM(prhs[0]);
    n_sensors=(mwSize)mxGetN(prhs[0]);
    
    filt=(double*)mxGetPr(prhs[1]);
    l_filt=(mwSize)mxGetN(prhs[1]);
    
    l_x_p=n_samples+2*l_filt;
    l_filt_p=l_x_p
        
    
    
    
    plhs[0]=mxCreateDoubleMatrix(n_samples,n_sensors,mxReal);
    y=(double*)mxGetPr(plhs[0]);
    

}