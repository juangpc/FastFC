#include <mex.h>
#include <stdio.h>
#include <fftw3.h>


void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    
    mwSize i,aux_i;
    float *b,*x;
    double *x_d,*b_d,*y;
    fftwf_complex *fft_b,*fft_x;
    fftwf_plan p_r2c,p_c2r;
    
    const mwSize l_b=(mwSize)mxGetN(prhs[0]);
    const mwSize n_samples=(mwSize)mxGetM(prhs[1]);
    const mwSize n_sensors=(mwSize)mxGetN(prhs[1]);
    
    const mwSize l_p=n_samples+2*(l_b-1);
    const int is_even=(l_p%2)?0:1;

    b_d=(double*)mxGetPr(prhs[0]);    
    x_d=(double*)mxGetPr(prhs[1]);
    
    b=(float*)mxMalloc(l_p*sizeof(float));
    x=(float*)mxMalloc(l_p*sizeof(float));
    fft_x=(fftwf_complex*)fftwf_malloc(l_p*sizeof(fftwf_complex));
    fft_b=(fftwf_complex*)fftwf_malloc(l_p*sizeof(fftwf_complex));
    
    p_r2c=fftwf_plan_dft_r2c_1d(l_p,b,fft_b,FFTW_MEASURE);
    p_c2r=fftwf_plan_dft_c2r_1d(l_p,fft_b,b,FFTW_MEASURE);
    
    //copy b with normal padding
    for(i=0;i<l_b;i++)
        b[i]=(float)b_d[i];
    for(;i<l_p;i++)
        b[i]=(float)0.;
    //execute fft of b
    fftwf_execute_dft_r2c(p_r2c,b,fft_b);
    
    //calculate filter mask and store in b (as it is Real)
    for(i=0;i<l_p/2+1;i++)
        b[i]=fft_b[i][0]*fft_b[i][0]+fft_b[i][1]*fft_b[i][1];
    
    //no need for fft_b any more
    fftwf_free(fft_b);

    //copy sensor with circular padding
    for(i=0;i<l_b-2;i++)
        x[i]=(2.*(float)x_d[0])-(float)x_d[l_b-2-i];
    for(i=0;i<n_samples;i++)
        x[l_b-2+i]=(float)x_d[i];
    aux_i=1;
    for(;i<l_p;i++)
        x[l_b-2+i]=(2.*(float)x_d[n_samples-1])-(float)x_d[n_samples-1-aux_i++];
    //calc fft of sensor
    fftwf_execute_dft_r2c(p_r2c,x,fft_x);
    
    //filter in freq domain
    for(i=0;i<l_p/2+1;i++)
    {
        fft_x[i][0]*=b[i];
        fft_x[i][1]*=b[i];
    }
    
    fftwf_execute_dft_c2r(p_c2r,fft_x,x);

    plhs[0]=mxCreateDoubleMatrix(n_samples,n_sensors,mxREAL);
    y=(double*)mxGetPr(plhs[0]);
    
    aux_i=0;
    for(i=l_b-2;i<l_p-l_b;i++)
        y[aux_i++]=(double)(x[i]/(float)l_p);

    fftwf_free(fft_x);
    fftwf_destroy_plan(p_r2c);
    fftwf_destroy_plan(p_c2r);
    mxFree(x);
    mxFree(b);
}