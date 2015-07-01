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
#include <fftw3.h>
#include <omp.h>

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    int sensor_i,thread_i;
    mwSize i;
    float *b,**x;
    double *x_d,*b_d,*y;
    fftwf_complex **fft_x,*fft_b;
    fftwf_plan p_r2c,p_c2r;
    
    const mwSize l_b=(mwSize)mxGetN(prhs[0]);
    const mwSize n_samples=(mwSize)mxGetM(prhs[1]);
    const int n_sensors=(int)mxGetN(prhs[1]);
    const int mode=(int)mxGetScalar(prhs[2]);
    const int n_threads=(int)mxGetScalar(prhs[3]);
    const mwSize l_p=n_samples+2*(l_b-1);
    const float n_factor=(float)1./l_p;

    b_d=(double*)mxGetPr(prhs[0]);    
    x_d=(double*)mxGetPr(prhs[1]);
    
    b=(float*)mxMalloc(l_p*sizeof(float));
    x=(float**)mxMalloc(n_threads*sizeof(float*));
    fft_x=(fftwf_complex**)mxMalloc(n_threads*sizeof(fftwf_complex*));
    fft_b=(fftwf_complex*)fftwf_malloc(l_p*sizeof(fftwf_complex));
    
    plhs[0]=mxCreateDoubleMatrix(n_samples,n_sensors,mxREAL);
    y=(double*)mxGetPr(plhs[0]);

    for(thread_i=0;thread_i<n_threads;thread_i++)
    {
        x[thread_i]=(float*)mxMalloc(l_p*sizeof(float));
        fft_x[thread_i]=(fftwf_complex*)fftwf_malloc(l_p*sizeof(fftwf_complex));
    }
    
    if(mode==3)
    {
        p_r2c=fftwf_plan_dft_r2c_1d((int)l_p,b,fft_b,FFTW_EXHAUSTIVE);
        p_c2r=fftwf_plan_dft_c2r_1d((int)l_p,fft_b,b,FFTW_EXHAUSTIVE);
    }
    else if(mode==2)
    {
        p_r2c=fftwf_plan_dft_r2c_1d((int)l_p,b,fft_b,FFTW_MEASURE);
        p_c2r=fftwf_plan_dft_c2r_1d((int)l_p,fft_b,b,FFTW_MEASURE);
    }
    else
    {
        p_r2c=fftwf_plan_dft_r2c_1d((int)l_p,b,fft_b,FFTW_ESTIMATE);
        p_c2r=fftwf_plan_dft_c2r_1d((int)l_p,fft_b,b,FFTW_ESTIMATE);
    }                                           

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
	
	fftwf_free(fft_b);  

    omp_set_num_threads(n_threads);
    #pragma omp parallel
    {
        #pragma omp for private(sensor_i,i,thread_i)
        for(sensor_i=0;sensor_i<n_sensors;sensor_i++)
        {
            thread_i=omp_get_thread_num();

            //copy sensor with circular padding
            for(i=0;i<l_b-2;i++)
                x[thread_i][i]=((float)2.*(float)x_d[sensor_i*n_samples])-
                        (float)x_d[sensor_i*n_samples+l_b-2-i];
            for(;i<l_p-l_b;i++)
                x[thread_i][i]=(float)x_d[sensor_i*n_samples+i-l_b+2];
            for(;i<l_p;i++)
                x[thread_i][i]=((float)2.*(float)x[thread_i][l_p-l_b-1])-
                        (float)x[thread_i][2*(l_p-l_b-1)-i];

            //calc fft of sensor
            fftwf_execute_dft_r2c(p_r2c,x[thread_i],fft_x[thread_i]);

            //filter in freq domain 
            for(i=0;i<l_p/2+1;i++)
            {
                fft_x[thread_i][i][0]*=b[i];
                fft_x[thread_i][i][1]*=b[i];
            }

            fftwf_execute_dft_c2r(p_c2r,fft_x[thread_i],x[thread_i]);

            for(i=0;i<n_samples;i++)
                y[sensor_i*n_samples+i]=(double)(x[thread_i][i+l_b-2]*n_factor);
        }
    }
    
    fftwf_destroy_plan(p_r2c);
	fftwf_destroy_plan(p_c2r);
    for(thread_i=0;thread_i<n_threads;thread_i++)
    {   
        mxFree(x[thread_i]);
        fftwf_free(fft_x[thread_i]);
    }
    mxFree(fft_x);    
    mxFree(x);
    mxFree(b);    

}