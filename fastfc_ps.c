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

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    char *wisdom_str;
    int n_samples,n_sensors,n_elements,wisdom_str_length,n_f,ii,i,
            sensor_i,sensor_j,n_indexes;
    float factor_n,*x_f,*x_f_i,*x_f_j,*adj_plv,*adj_pli,*adj_wpli,
            sum_cos,sum_sin,sum_sin_sign,phi,sin_phi,sum_abs_sin;
    double *x_d,*out_plv,*out_pli,*out_wpli;
    fftwf_complex *a,*h,*h_i;
    fftwf_plan p_r2c, p_c2c;
    
    x_d=(double*)mxGetPr(prhs[0]);
    n_samples=(int)mxGetM(prhs[0]);
    n_sensors=(int)mxGetN(prhs[0]);
    n_elements=n_samples*n_sensors;
    factor_n=(float)1.0/n_samples;

    x_f=(float*)mxMalloc(n_elements*sizeof(float));
    i=n_elements;
    while(i--)
        x_f[i]=(float)x_d[i];
    
    wisdom_str_length=(int)mxGetN(prhs[1])*sizeof(mxChar)+1;
    wisdom_str=(char*)mxMalloc(wisdom_str_length);
    mxGetString(prhs[1], wisdom_str, (mwSize)wisdom_str_length);
    fftwf_import_wisdom_from_string(wisdom_str);

    h=(fftwf_complex*)fftwf_malloc(n_samples*n_sensors*sizeof(fftwf_complex));
    a=(fftwf_complex*)fftwf_malloc(n_samples*n_sensors*sizeof(fftwf_complex));
    
    p_r2c=fftwf_plan_many_dft_r2c(1,&n_samples,n_sensors,
            x_f,&n_samples,1,n_samples,
            h,&n_samples,1,n_samples,
            FFTW_ESTIMATE);
    p_c2c=fftwf_plan_many_dft(1,&n_samples,n_sensors,
            h,&n_samples,1,n_samples,
            a,&n_samples,1,n_samples,
            FFTW_BACKWARD,FFTW_ESTIMATE);
    
    fftwf_execute(p_r2c);

    for(sensor_i=0;sensor_i<n_sensors;++sensor_i)
    {
        h_i=h+(sensor_i*n_samples);
        if(n_samples%2)
        {
            n_f=n_samples/2+1;
            for(ii=1;ii<n_f;++ii)
            {
                h_i[ii][0]*=2;
                h_i[ii][1]*=2;
            }
        }else{
            n_f=n_samples/2;
            for(ii=1;ii<n_f;++ii)
            {
                h_i[ii][0]*=2;
                h_i[ii][1]*=2;
            }
            ++ii;
        }
        for(;ii<n_samples;++ii)
        {
            h_i[ii][0]=0;
            h_i[ii][1]=0;
        }
    }
    
    fftwf_execute(p_c2c);
    
    mxFree(wisdom_str);
    fftwf_free(h);
    fftwf_destroy_plan(p_r2c);
    fftwf_destroy_plan(p_c2c);

    for(i=0;i<n_elements;++i)
        x_f[i]=atan2f(factor_n*a[i][1],factor_n*a[i][0]);
    
    fftwf_free(a);
    
    adj_plv =(float*)mxMalloc(n_sensors*n_sensors*sizeof(float));
    adj_pli =(float*)mxMalloc(n_sensors*n_sensors*sizeof(float));
    adj_wpli=(float*)mxMalloc(n_sensors*n_sensors*sizeof(float));
    
    for(sensor_i=0;sensor_i<n_sensors;++sensor_i)
    {
        adj_plv[sensor_i+(n_sensors*sensor_i)]=1.0;
        adj_pli[sensor_i+(n_sensors*sensor_i)]=1.0;
        adj_wpli[sensor_i+(n_sensors*sensor_i)]=1.0;

        for(sensor_j=sensor_i+1;sensor_j<n_sensors;++sensor_j)
        {
            x_f_i=x_f+(sensor_i*n_samples);
            x_f_j=x_f+(sensor_j*n_samples);
            sum_cos=0;
            sum_sin=0;
            sum_sin_sign=0;
            sum_abs_sin=0;
            
            for(i=0;i<n_samples;++i)
            {
               phi=x_f_i[i]-x_f_j[i];
               sum_cos+=cosf(phi);
               sin_phi=sinf(phi);
               sum_sin+=sin_phi;
               sum_sin_sign+=(float)(sin_phi>0.0)-(sin_phi<0.0);
               sum_abs_sin+=fabs(sin_phi);
            }
            sum_cos/=n_samples;
            sum_sin/=n_samples;
                        
            adj_plv[sensor_j+(n_sensors*sensor_i)]=
                    sqrt(sum_cos*sum_cos+sum_sin*sum_sin);
            adj_pli[sensor_j+(n_sensors*sensor_i)]=fabs(sum_sin_sign/n_samples);
            adj_wpli[sensor_j+(n_sensors*sensor_i)]=fabs(sum_sin)/sum_abs_sin*n_samples;
            
            adj_plv[sensor_i+(n_sensors*sensor_j)]=adj_plv[sensor_j+(n_sensors*sensor_i)];
            adj_pli[sensor_i+(n_sensors*sensor_j)]=adj_pli[sensor_j+(n_sensors*sensor_i)];
            adj_wpli[sensor_i+(n_sensors*sensor_j)]=adj_wpli[sensor_j+(n_sensors*sensor_i)];
        }
    }
    
    mxFree(x_f);
    plhs[0]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);

    n_indexes=n_sensors*n_sensors;

    out_plv=(double*)mxGetPr(plhs[0]);
    for(i=0;i<n_indexes;++i)
        out_plv[i]=(double)adj_plv[i];

    out_pli=(double*)mxGetPr(plhs[1]);
    for(i=0;i<n_indexes;++i)
        out_pli[i]=(double)adj_pli[i];

    out_wpli=(double*)mxGetPr(plhs[2]);
    for(i=0;i<n_indexes;++i)
        out_wpli[i]=(double)adj_wpli[i];

    mxFree(adj_plv);
    mxFree(adj_pli);
    mxFree(adj_wpli);
}