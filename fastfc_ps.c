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
    int i,f_i,sample_i,sensor_i,first_samp_sensor_i,sensor_j;
    float *x_f,*x_f_i,*x_f_j,*adj_plv,*adj_pli,*adj_wpli,*pval_plv,plv_i,
            sum_cos,sum_sin,sum_sin_sign,phi,sin_phi,sum_abs_sin;
    double *x_d,*out_plv,*out_pli,*out_wpli,*out_pval;
    fftwf_complex *a,*h;
    fftwf_plan p_r2c, p_c2c;
    
    const int n_samples=(int)mxGetM(prhs[0]);
    const int n_sensors=(int)mxGetN(prhs[0]);
    const int n_elements=n_samples*n_sensors;
    const int is_even=(n_samples%2)? 0:1;
    const int n_f=n_samples/2+((is_even)? 0:1);
    const float factor_n=(float)1.0/n_samples;
    const int n_indexes=n_sensors*n_sensors;
    const int init_sample=(int)mxGetScalar(prhs[1]);
    const int n_threads=omp_get_num_procs();
    const int last_sample=n_samples-init_sample;
    const int n_samples_eff=n_samples-2*init_sample;
            
    x_d=(double*)mxGetPr(prhs[0]);
    
    x_f=(float*)mxMalloc(n_elements*sizeof(float));
    h=(fftwf_complex*)fftwf_malloc(n_elements*sizeof(fftwf_complex));
    a=(fftwf_complex*)fftwf_malloc(n_elements*sizeof(fftwf_complex));
    
    p_r2c=fftwf_plan_dft_r2c_1d(n_samples,x_f,h,FFTW_EXHAUSTIVE);
    p_c2c=fftwf_plan_dft_1d(n_samples,h,a,FFTW_BACKWARD,FFTW_EXHAUSTIVE);
    
    omp_set_num_threads(n_threads);
    #pragma omp parallel
    {
        #pragma omp for private(sensor_i,sample_i,first_samp_sensor_i,f_i)
        for(sensor_i=0;sensor_i<n_sensors;++sensor_i)
        {
            first_samp_sensor_i=sensor_i*n_samples;
            for(sample_i=0;sample_i<n_samples;++sample_i)
                x_f[sample_i+first_samp_sensor_i]=(float)x_d[sample_i+first_samp_sensor_i];
            
            fftwf_execute_dft_r2c(p_r2c,&x_f[first_samp_sensor_i],&h[first_samp_sensor_i]);
            
            for(f_i=1;f_i<n_f;++f_i)
            {
                h[f_i+first_samp_sensor_i][0]*=2;
                h[f_i+first_samp_sensor_i][1]*=2;
            }
            
            if(is_even)
                ++f_i;
            
            for(;f_i<n_samples;++f_i)
            {
                h[f_i+first_samp_sensor_i][0]=0;
                h[f_i+first_samp_sensor_i][1]=0;
            }
            
            fftwf_execute_dft(p_c2c,&h[first_samp_sensor_i],&a[first_samp_sensor_i]);
            
            for(sample_i=0;sample_i<n_samples;++sample_i)
                x_f[sample_i+first_samp_sensor_i]=atan2f(
                        factor_n*a[sample_i+first_samp_sensor_i][1],
                        factor_n*a[sample_i+first_samp_sensor_i][0]);
        }
    }
    
    fftwf_free(h);
    fftwf_free(a);
    fftwf_destroy_plan(p_r2c);
    fftwf_destroy_plan(p_c2c);
    
    adj_plv =(float*)mxMalloc(n_indexes*sizeof(float));
    pval_plv=(float*)mxMalloc(n_indexes*sizeof(float));
    adj_pli =(float*)mxMalloc(n_indexes*sizeof(float));
    adj_wpli=(float*)mxMalloc(n_indexes*sizeof(float));
        
    for(sensor_i=0;sensor_i<n_sensors;++sensor_i)
    {
        adj_plv[sensor_i+(n_sensors*sensor_i)]=1.0;
        adj_pli[sensor_i+(n_sensors*sensor_i)]=1.0;
        adj_wpli[sensor_i+(n_sensors*sensor_i)]=1.0;
        pval_plv[sensor_i+(n_sensors*sensor_i)]=0.0;
        
        #pragma omp parallel
        {
            #pragma omp for private(sensor_j,i,phi,sin_phi,x_f_i,x_f_j,plv_i) reduction(+:sum_cos,sum_sin,sum_sin_sign,sum_abs_sin) //private(sensor_j,i,phi,sin_phi,x_f_i,x_f_j,sum_cos,sum_sin,sum_sin_sign,sum_abs_sin)
            for(sensor_j=sensor_i+1;sensor_j<n_sensors;++sensor_j)
            {
                x_f_i=x_f+(sensor_i*n_samples);
                x_f_j=x_f+(sensor_j*n_samples);
                sum_cos=0;
                sum_sin=0;
                sum_sin_sign=0;
                sum_abs_sin=0;
                
                for(i=init_sample;i<last_sample;++i)
                {
                    phi=x_f_i[i]-x_f_j[i];
                    sum_cos+=cosf(phi);
                    sin_phi=sinf(phi);
                    sum_sin+=sin_phi;
                    sum_sin_sign+=(float)(sin_phi>0.0)-(sin_phi<0.0);
                    sum_abs_sin+=fabsf(sin_phi);
                }
                sum_cos/=n_samples_eff;
                sum_sin/=n_samples_eff;
                
                plv_i=sqrtf(sum_cos*sum_cos+sum_sin*sum_sin);
                adj_pli[sensor_j+(n_sensors*sensor_i)]=fabsf(sum_sin_sign/n_samples_eff);
                adj_wpli[sensor_j+(n_sensors*sensor_i)]=(fabsf(sum_sin)*n_samples_eff)/sum_abs_sin;
                adj_plv[sensor_j+(n_sensors*sensor_i)]=plv_i;
                pval_plv[sensor_j+(n_sensors*sensor_i)]=
 expf(sqrtf((float)(1+4*n_samples_eff)+4.0*(float)n_samples_eff*(float)n_samples_eff*(1.0-plv_i*plv_i))-(float)(1+2*n_samples_eff));
                
                adj_plv[sensor_i+(n_sensors*sensor_j)]=plv_i;
                adj_pli[sensor_i+(n_sensors*sensor_j)]=adj_pli[sensor_j+(n_sensors*sensor_i)];
                adj_wpli[sensor_i+(n_sensors*sensor_j)]=adj_wpli[sensor_j+(n_sensors*sensor_i)];
                pval_plv[sensor_i+(n_sensors*sensor_j)]=pval_plv[sensor_j+(n_sensors*sensor_i)];
                
            }
        }
    }

    mxFree(x_f);
    
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
    mxFree(adj_pli);
    mxFree(adj_wpli);
}
    
    
    
    