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
#include <math.h>
#include <omp.h>

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{   
    int i,s_i,s_j,th_i;
    float *b,*x,*adj_plv,*pval_plv,*adj_pli,*adj_wpli,phi,sin_phi,plv_i;
    float sum_cos=0;
    float sum_sin=0;
    float sum_sin_sign=0;
    float sum_abs_sin=0;
    
    double *x_d,*b_d,*out_plv,*out_pli,*out_wpli,*out_pval;
    fftwf_complex *fft_b,**h,**a;
    fftwf_plan p_r2c_1,p_r2c_2,p_c2r,p_c2c;
    
    const int n_samples=(int)mxGetM(prhs[0]);
    const int n_sensors=(int)mxGetN(prhs[0]);
    const int init_sample=(int)mxGetScalar(prhs[1]);
	const int l_b=(int)mxGetN(prhs[2]);
	const int mode=(int)mxGetScalar(prhs[3]);
	const int n_threads=1;//omp_get_num_procs();
	
    const int n_indexes=n_sensors*n_sensors;
	const int n_samples_eff=n_samples-2*init_sample;
    const int l_p=n_samples+2*(l_b-1);
    const int is_even=(n_samples%2)?0:1;
    const int n_f=n_samples/2+((is_even)?0:1);
    const float factor_n_l_p=(float)1./l_p;
    const float factor_n_ns=(float)1./n_samples;
    const float factor_n_eff=(float)1./n_samples_eff;
            
    x_d=(double*)mxGetPr(prhs[0]);
    b_d=(double*)mxGetPr(prhs[2]);
    
    b=(float*)fftwf_malloc(l_p*sizeof(float));
    fft_b=(fftwf_complex*)fftwf_malloc(l_p*sizeof(fftwf_complex));
    
    x=(float*)mxMalloc(l_p*n_sensors*sizeof(float));
    h=(fftwf_complex**)mxMalloc(n_threads*sizeof(fftwf_complex*));
    a=(fftwf_complex**)mxMalloc(n_threads*sizeof(fftwf_complex*));

    for(th_i=0;th_i<n_threads;th_i++)
    {
        h[th_i]=(fftwf_complex*)fftwf_malloc(l_p*sizeof(fftwf_complex));
        a[th_i]=(fftwf_complex*)fftwf_malloc(n_samples*sizeof(fftwf_complex));
    }
    
    if(mode==3)
    {
        p_r2c_1=fftwf_plan_dft_r2c_1d(l_p,b,fft_b,FFTW_EXHAUSTIVE);
        p_r2c_2=fftwf_plan_dft_r2c_1d(n_samples,b,fft_b,FFTW_EXHAUSTIVE);
        p_c2r  =fftwf_plan_dft_c2r_1d(l_p,*h,x,FFTW_EXHAUSTIVE);
        p_c2c  =fftwf_plan_dft_1d(n_samples,*h,*a,FFTW_BACKWARD,FFTW_EXHAUSTIVE);
    }else if(mode==2)
    {
        p_r2c_1=fftwf_plan_dft_r2c_1d(l_p,b,fft_b,FFTW_PATIENT);
        p_r2c_2=fftwf_plan_dft_r2c_1d(n_samples,x,*h,FFTW_PATIENT);
        p_c2r  =fftwf_plan_dft_c2r_1d(l_p,*h,x,FFTW_PATIENT);
        p_c2c  =fftwf_plan_dft_1d(n_samples,*h,*a,FFTW_BACKWARD,FFTW_PATIENT);
    }else if(mode==1)
    {
        p_r2c_1=fftwf_plan_dft_r2c_1d(l_p,b,fft_b,FFTW_MEASURE);
        p_r2c_2=fftwf_plan_dft_r2c_1d(n_samples,b,fft_b,FFTW_MEASURE);
        p_c2r  =fftwf_plan_dft_c2r_1d(l_p,*h,x,FFTW_MEASURE);
        p_c2c  =fftwf_plan_dft_1d(n_samples,*h,*a,FFTW_BACKWARD,FFTW_MEASURE);
    }else
    {
        p_r2c_1=fftwf_plan_dft_r2c_1d(l_p,b,fft_b,FFTW_ESTIMATE);
        p_r2c_2=fftwf_plan_dft_r2c_1d(n_samples,b,fft_b,FFTW_ESTIMATE);
        p_c2r  =fftwf_plan_dft_c2r_1d(l_p,*h,x,FFTW_ESTIMATE);
        p_c2c  =fftwf_plan_dft_1d(n_samples,*h,*a,FFTW_BACKWARD,FFTW_ESTIMATE);
    }        
    
    for(i=0;i<l_b;i++)
        b[i]=(float)b_d[i];
    for(;i<l_p;i++)
        b[i]=(float)0.;
    
    fftwf_execute_dft_r2c(p_r2c_1,b,fft_b);
        
    for(i=0;i<l_p/2+1;i++)
        b[i]=fft_b[i][0]*fft_b[i][0]+fft_b[i][1]*fft_b[i][1];
    
    fftwf_free(fft_b);
    
	omp_set_num_threads(n_threads);
    #pragma omp parallel private(i,s_i,th_i) shared(n_sensors,n_samples,l_b,l_p,n_f)
    {
        #pragma omp for 
        for(s_i=0;s_i<n_sensors;s_i++)
        {
            th_i=omp_get_thread_num();
            
            for(i=0;i<l_b-2;i++)
                x[s_i*l_p+i]=((float)2.*(float)x_d[s_i*n_samples])-(float)x_d[s_i*n_samples+l_b-2-i];
            for(;i<l_p-l_b;i++)
                x[s_i*l_p+i]=(float)x_d[s_i*n_samples+i-l_b+2];
            for(;i<l_p;i++)
                x[s_i*l_p+i]=((float)2.*x[s_i*l_p+l_p-l_b-1])-x[s_i*l_p+2*(l_p-l_b-1)-i];
            
            fftwf_execute_dft_r2c(p_r2c_1,&x[s_i*l_p],h[th_i]);
            
            for(i=0;i<l_p/2+1;i++)
            {
                h[th_i][i][0]*=b[i];
                h[th_i][i][1]*=b[i];
            }
            fftwf_execute_dft_c2r(p_c2r,h[th_i],&x[s_i*l_p]);
            
            for(i=0;i<l_p;i++)
                x[s_i*l_p+i]*=factor_n_l_p;

            //el error está en esta función. Esto peta.
            fftwf_execute_dft_r2c(p_r2c_2,&x[(s_i*l_p)+l_b-2],h[th_i]);
            
            for(i=1;i<n_f;i++)
            {
                h[th_i][i][0]*=2;
                h[th_i][i][1]*=2;
            }
            if(is_even)
                i++;
            for(;i<n_samples;i++)
            {
                h[th_i][i][0]*=(float)0.;
                h[th_i][i][1]*=(float)0.;
            }
            
            fftwf_execute_dft(p_c2c,h[th_i],a[th_i]);
            
            for(i=0;i<n_samples_eff;i++)
                x[s_i*l_p+i]=atan2f(
                        factor_n_ns*a[th_i][i+init_sample][1],
                        factor_n_ns*a[th_i][i+init_sample][0]);
        }
    }
    
    fftwf_destroy_plan(p_r2c_1);
    fftwf_destroy_plan(p_r2c_2);
    fftwf_destroy_plan(p_c2c);
    fftwf_destroy_plan(p_c2r);
    
    for(th_i=0;th_i<n_threads;th_i++)
    {
        fftwf_free(h[th_i]);
        fftwf_free(a[th_i]);
    }
    mxFree(h);
    mxFree(a);
    fftwf_free(b);
    
    adj_plv =(float*)mxMalloc(n_indexes*sizeof(float));
    pval_plv=(float*)mxMalloc(n_indexes*sizeof(float));
    adj_pli =(float*)mxMalloc(n_indexes*sizeof(float));
    adj_wpli=(float*)mxMalloc(n_indexes*sizeof(float));
        
    for(s_i=0;s_i<n_sensors;++s_i)
    {
        adj_plv[s_i+(n_sensors*s_i)]=1.0;
        adj_pli[s_i+(n_sensors*s_i)]=1.0;
        adj_wpli[s_i+(n_sensors*s_i)]=1.0;
        pval_plv[s_i+(n_sensors*s_i)]=0.0;
        
        #pragma omp parallel
        {
            #pragma omp for private(s_j,i,phi,sin_phi,plv_i) reduction(+:sum_cos,sum_sin,sum_sin_sign,sum_abs_sin) //private(sensor_j,i,phi,sin_phi,x_f_i,x_f_j,sum_cos,sum_sin,sum_sin_sign,sum_abs_sin)
            for(s_j=s_i+1;s_j<n_sensors;++s_j)
            {
                sum_cos=0;
                sum_sin=0;
                sum_sin_sign=0;
                sum_abs_sin=0;
                
                for(i=0;i<n_samples_eff;++i)
                {
                    phi=x[s_i*l_p+i]-x[s_j*l_p+i];
                    sum_cos+=cosf(phi);
                    sin_phi=sinf(phi);
                    sum_sin+=sin_phi;
                    sum_sin_sign+=(float)(sin_phi>0.0)-(sin_phi<0.0);
                    sum_abs_sin+=fabsf(sin_phi);
                }
                sum_cos*=factor_n_eff;
                sum_sin*=factor_n_eff;
                
                plv_i=sqrtf(sum_cos*sum_cos+sum_sin*sum_sin);
                adj_pli[s_j+(n_sensors*s_i)]=fabsf(sum_sin_sign/n_samples_eff);
                adj_wpli[s_j+(n_sensors*s_i)]=(fabsf(sum_sin)*n_samples_eff)/sum_abs_sin;
                adj_plv[s_j+(n_sensors*s_i)]=plv_i;
                pval_plv[s_j+(n_sensors*s_i)]=expf(
                        sqrtf((float)(1+4*n_samples_eff)+
                        (float)4.0*n_samples_eff*n_samples_eff*((float)1.0-plv_i*plv_i))
                        -(1+2*n_samples_eff));
                
                adj_plv[s_i+(n_sensors*s_j)]=plv_i;
                adj_pli[s_i+(n_sensors*s_j)]=adj_pli[s_j+(n_sensors*s_i)];
                adj_wpli[s_i+(n_sensors*s_j)]=adj_wpli[s_j+(n_sensors*s_i)];
                pval_plv[s_i+(n_sensors*s_j)]=pval_plv[s_j+(n_sensors*s_i)];
                
            }
        }
    }

    mxFree(x);
    
    plhs[0]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    plhs[3]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    
    
    #pragma omp parallel sections private(i) 
    {
        #pragma omp section
        {
            out_plv=(double*)mxGetPr(plhs[0]);
            for(i=0;i<n_indexes;++i)
                out_plv[i]=(double)adj_plv[i];
        }
        #pragma omp section 
        {
            out_pval=(double*)mxGetPr(plhs[1]);
            for(i=0;i<n_indexes;++i)
                out_pval[i]=(double)pval_plv[i];
        }

        #pragma omp section
        {
            out_pli=(double*)mxGetPr(plhs[2]);
            for(i=0;i<n_indexes;++i)
                out_pli[i]=(double)adj_pli[i]; 
        }
        #pragma omp section 
        {
            out_wpli=(double*)mxGetPr(plhs[3]);
            for(i=0;i<n_indexes;++i)
                out_wpli[i]=(double)adj_wpli[i];
        }
    }
    
    mxFree(adj_plv);
    mxFree(pval_plv);
    mxFree(adj_pli);
    mxFree(adj_wpli);
    
}