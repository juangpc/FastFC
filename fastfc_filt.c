#include <mex.h>
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
    const unsigned int n_sensors=(unsigned int)mxGetN(prhs[1]);
    const short int mode=(short int)mxGetScalar(prhs[2]);
            
    const mwSize l_p=n_samples+2*(l_b-1);
    const int is_even=(l_p%2)?0:1;

    b_d=(double*)mxGetPr(prhs[0]);    
    x_d=(double*)mxGetPr(prhs[1]);
    
    b=(float*)malloc(l_p*sizeof(float));
    x=(float*)malloc(l_p*sizeof(float));
    fft_x=(fftwf_complex*)fftwf_malloc(l_p*sizeof(fftwf_complex));
    fft_b=(fftwf_complex*)fftwf_malloc(l_p*sizeof(fftwf_complex));
    
    if(mode==2)
    {
        p_r2c=fftwf_plan_dft_r2c_1d(l_p,b,fft_b,FFTW_EXHAUSTIVE);
        p_c2r=fftwf_plan_dft_c2r_1d(l_p,fft_b,b,FFTW_EXHAUSTIVE);
    }
    else
    {
        p_r2c=fftwf_plan_dft_r2c_1d(l_p,b,fft_b,FFTW_ESTIMATE);
        p_c2r=fftwf_plan_dft_c2r_1d(l_p,fft_b,b,FFTW_ESTIMATE);
    }                                           

    //copy b with normal padding
    for(i=0;i<l_b;i++)
        b[i]=(float)b_d[i];
    for(;i<l_p;i++)
        b[i]=(float)0.;
    //execute fft of b
    fftwf_execute_dft_r2c(p_r2c,b,fft_b);
    
    //copy sensor with circular padding
    for(i=0;i<l_b-2;i++)
        x[i]=(2.*(float)x_d[0])-(float)x_d[l_b-2-i];
    for(;i<l_p-l_b;i++)
        x[i]=(float)x_d[i-l_b+2];
    for(;i<l_p;i++)
        x[i]=(2.*x[l_p-l_b-1])-x[2*(l_p-l_b-1)-i];
    
    //calc fft of sensor
    fftwf_execute_dft_r2c(p_r2c,x,fft_x);
    
    //calculate filter mask and 
    //filter in freq domain 
    for(i=0;i<l_p/2+1;i++)
    {
        fft_x[i][0]*=fft_b[i][0]*fft_b[i][0]+fft_b[i][1]*fft_b[i][1];
        fft_x[i][1]*=fft_b[i][0]*fft_b[i][0]+fft_b[i][1]*fft_b[i][1];
    }
    
    fftwf_execute_dft_c2r(p_c2r,fft_x,x);
    plhs[0]=mxCreateDoubleMatrix(n_samples,n_sensors,mxREAL);
    y=(double*)mxGetPr(plhs[0]);
    
    aux_i=0;
    for(i=l_b-2;i<l_p-l_b;i++)
        y[aux_i++]=(double)(x[i]/l_p);

    mxFree(b);
    mxFree(x);
    fftwf_destroy_plan(p_c2r);
    fftwf_destroy_plan(p_r2c);
    fftwf_free(fft_b);    
    fftwf_free(fft_x);

}