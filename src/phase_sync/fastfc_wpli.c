//
//    Copyright (c) 2014-2016 Juan Garcia-Prieto Cuesta
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
#include <fftw3.h>

#define _2PI 6.283185307179586

void calc_hamm_win(float*w,int wl)
{
    int i;
    const int wlh=(wl%2)? wl/2:(wl+1)/2;
    const float k=(float)_2PI/(wl-1);
    for(i=0;i<wlh;i++)
    {
        w[i]=(float)0.54-(float)0.46*cosf(k*i);
        w[wl-i-1]=w[i];
    }
}

float norm_squared(float*vec,int l)
{
    int i;
    float dist_sq=0;
    for(i=0;i<l;i++)
        dist_sq += vec[i]*vec[i];
    return dist_sq;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    int first_samp_in_win,win_i,sens_i,sens_j,win_samp_i,f_i,nii;
    float *x_f,*win,*pxy,*pxx,k_coh,win_mean,sum_,sum_abs;
    double *x_d,*wpli,*imC;
    fftwf_complex *h;
    fftwf_plan ftplan;
    
    const int n_samp=(int)mxGetM(prhs[0]);
    const int n_sens=(int)mxGetN(prhs[0]);
	const int win_l=(n_samp<5)? n_samp:(int)(n_samp*0.2222);
	const int n_overlap=(int)(n_samp<5)? 0:(win_l/2);
	const int n_wins=(int)(n_samp<5)? 1:(n_samp-n_overlap)/(win_l-n_overlap);
	const int is_even=(n_samp%2)?0:1;
    const int n_indices=(int)(n_sens*n_sens-n_sens)/2;
	const int f_select=(win_l%2)?((win_l+1)/2):(win_l/2+1);
    const int n_threads=omp_get_num_procs();
    
    x_d=(double*)mxGetPr(prhs[0]);
    win=(float*)mxMalloc(win_l*sizeof(float));
    
    x_f=(float*)fftwf_malloc(win_l*n_sens*sizeof(float));
    h=(fftwf_complex*)fftwf_malloc(win_l*n_sens*sizeof(fftwf_complex));
    ftplan=fftwf_plan_dft_r2c_1d(win_l,x_f,h,FFTW_MEASURE);
    
    pxy=(float*)mxCalloc(f_select*n_indices,sizeof(float));
    pxx=(float*)mxCalloc(n_sens,sizeof(float));
    
    calc_hamm_win(win,win_l);
    k_coh=((float)1.)/(sqrtf((float)f_select)*n_wins*norm_squared(win,win_l));
    first_samp_in_win=0;
    
    for(win_i=0;win_i<n_wins;win_i++)
    {
        #pragma omp parallel
        {
        
            #pragma omp for private(sens_i,win_mean,win_samp_i)
            for(sens_i=0;sens_i<n_sens;sens_i++)
            {
                win_mean=0;
                for(win_samp_i=0;win_samp_i<win_l;win_samp_i++)
                    win_mean+=(float)x_d[sens_i*n_samp+first_samp_in_win+win_samp_i];

                win_mean/=win_l;

                for(win_samp_i=0;win_samp_i<win_l;win_samp_i++)
                    x_f[sens_i*win_l+win_samp_i]=win[win_samp_i]*
                            ((float)x_d[sens_i*n_samp+first_samp_in_win+win_samp_i]-win_mean);

                fftwf_execute_dft_r2c(ftplan,x_f+sens_i*win_l,h+sens_i*win_l);
            }

            #pragma omp for private(sens_i,sens_j,f_i,nii)
            for(sens_i=0;sens_i<n_sens;sens_i++)
            {
                for(f_i=0;f_i<f_select;f_i++)
                    pxx[sens_i] +=
                            h[sens_i*win_l+f_i][0]*h[sens_i*win_l+f_i][0] +
                            h[sens_i*win_l+f_i][1]*h[sens_i*win_l+f_i][1];

                nii=(sens_i+1)*(sens_i+2)/2;
                for(sens_j=sens_i+1;sens_j<n_sens;sens_j++)
                    for(f_i=0;f_i<f_select;f_i++)
                        pxy[(sens_i*n_sens+sens_j-nii)*f_select+f_i] +=
                                h[sens_i*win_l+f_i][0]*h[sens_j*win_l+f_i][1] -
                                h[sens_i*win_l+f_i][1]*h[sens_j*win_l+f_i][0];
            }
        }
        
        first_samp_in_win+=win_l-n_overlap;
    }
    
    
    mxFree(win);
    fftwf_free(x_f);
    fftwf_free(h);
    fftwf_destroy_plan(ftplan);
    
    plhs[0]=mxCreateDoubleMatrix(n_sens,n_sens,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(n_sens,n_sens,mxREAL);
    wpli=(double*)mxGetPr(plhs[0]);
    imC=(double*)mxGetPr(plhs[1]);
    
    sum_=(float)0.;
    sum_abs=(float)0.;

#pragma omp parallel
    {
#pragma omp for private(sens_i,sens_j,nii,f_i) reduction(+:sum_,sum_abs)
        for(sens_i=0;sens_i<n_sens;sens_i++)
        {
            nii=(sens_i+1)*(sens_i+2)/2;
            wpli[sens_i*n_sens+sens_i]=0.;
            imC[sens_i*n_sens+sens_i]=0.;
            
            for(sens_j=sens_i+1;sens_j<n_sens;sens_j++)
            {
                sum_=(float)0.;
                sum_abs=(float)0.;
                for(f_i=0;f_i<f_select;f_i++)
                {
                    sum_ +=          pxy[(sens_i*n_sens+sens_j-nii)*f_select+f_i];
                    sum_abs += fabsf(pxy[(sens_i*n_sens+sens_j-nii)*f_select+f_i]);
                }
                
                wpli[sens_i*n_sens+sens_j]=(double)(fabsf(sum_)/sum_abs);
                imC[sens_i*n_sens+sens_j]=(double) k_coh*sum_/sqrtf(pxx[sens_i] + pxx[sens_j]);
                
                wpli[sens_j*n_sens+sens_i]=wpli[sens_i*n_sens+sens_j];
                imC[sens_j*n_sens+sens_i]=imC[sens_i*n_sens+sens_j];
            }
        }
    }
    mxFree(pxx);
    mxFree(pxy);
}

