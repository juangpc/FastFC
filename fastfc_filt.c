#include <stdio.h>
#include <mex.h>
#include <fftw3.h>

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    unsigned int n_sensors,sensor_i,N,aux_i;
    mwSize n_samples,sample_i,l_p;
    double *x_d,*b_d,*y;
    float *si,*b_p;
    fftwf_complex *fft_b_p,*fft_si;
    fftwf_plan p_r2c, p_c2r;
    
    x_d=(double*)mxGetPr(prhs[0]);
    n_samples=(mwSize)mxGetM(prhs[0]);
    n_sensors=(mwSize)mxGetN(prhs[0]);
    
    b_d=(double*)mxGetPr(prhs[1]);
    N=(int)mxGetN(prhs[1]); //filter order is N-1
    
    l_p=n_samples+2*N;
    
    b_p=(float*)mxMalloc(l_p*sizeof(float));
    fft_b_p=(fftwf_complex*)fftwf_malloc(l_p*sizeof(fftwf_complex));
    
    si=(float*)mxMalloc(l_p*sizeof(float));
    fft_si=(fftwf_complex*)fftwf_malloc(l_p*sizeof(fftwf_complex));

    p_r2c=fftwf_plan_dft_r2c_1d(l_p,b_p,fft_b_p,FFTW_MEASURE);
    p_c2r=fftwf_plan_dft_c2r_1d(l_p,fft_b_p,b_p,FFTW_MEASURE);
   
    for(sample_i=0;sample_i<N;sample_i++)
        b_p[sample_i]=(float)b_d[sample_i];
    for(;sample_i<l_p;sample_i++)
        b_p[sample_i]=0.;
    
    fftwf_execute_dft_r2c(p_r2c,b_p,fft_b_p);
    
    ///////
    
    for(sample_i=0;sample_i<l_p;sample_i++)
        b_p[sample_i]=fft_b_p[sample_i][0]*fft_b_p[sample_i][0]-fft_b_p[sample_i][1]*fft_b_p[sample_i][1];
    
    fftwf_free(fft_b_p);
    
    for(sample_i=0;sample_i<N-1;sample_i++)
        si[sample_i]=2.*(float)x_d[0]-(float)x_d[N-1-sample_i];
    for(sample_i=0;sample_i<n_samples;sample_i++)
        si[sample_i+N-1]=(float)x_d[sample_i];
    aux_i=1;
    for(;sample_i<l_p;sample_i++)
        si[sample_i]=2.*(float)x_d[n_samples-1]-(float)x_d[n_samples-1-aux_i++];
    
    fftwf_execute_dft_r2c(p_r2c,si,fft_si);
    
    for(sample_i=0;sample_i<l_p;sample_i++)
        fft_si[sample_i][0]=fft_si[sample_i][0]*b_p[sample_i];
    
    fftwf_execute_dft(p_r2c,fft_si,si);
    
    plhs[0]=mxCreateDoubleMatrix(n_samples,n_sensors,mxREAL);
    y=(double*)mxGetPr(plhs[0]);
    
    aux_i=0;
    for(sample_i=N-1;sample_i<N-1+n_samples;sample_i++)
        y[aux_i++]=b_p[sample_i]*fft_si[sample_i][0];
                
    fftwf_free(fft_si);
    fftwf_destroy_plan(p_r2c);
    fftwf_destroy_plan(p_c2r);
    mxFree(si);
    mxFree(b_p);
}