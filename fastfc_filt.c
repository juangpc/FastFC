
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