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
//	 -------------------------------------------------------
//   Please consider citing our work:
//  
//   Efficient computation of functional brain networks: towards real-time functional connectivity
//	 Frontiers in Neuroinformatics (2017) García-Prieto Juan, Bajo Ricardo, Pereda Ernesto
//
//   -------------------------------------------------------
//   Contact: Juan Garcia-Prieto    juangpc (at) gmail.com
//   -------------------------------------------------------
//
#include <mex.h>
#include <omp.h>
#include <math.h>

struct neighbor{
    int index;
    float dist;
};

int comp(const void *a,const void *b)
{
    return ((*(struct neighbor *)a).dist<(*(struct neighbor *)b).dist)?
        -1:
        ((*(struct neighbor *)a).dist>(*(struct neighbor *)b).dist);
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    
    //call
    // shml(trial(by columns),emb_d,tau, k , w_emb , states_eff_step,
    //            (0)          (1)  (2) (3)   (4)       (5)          
    
    int state_eff_i,state_i,state_j,tau_i,n_states,n_states_eff,n_indexes,
            actual_state,k_in,index_i,sensor_i,sensor_j;
    float *mean_d_nn,*mean_d,*adj_s,*adj_h,*adj_m,*adj_l,
            dist_diff,dist_sum_1,dist_sum_2,cum_dist_sum_1,
            cum_sum_s,cum_sum_h,cum_sum_m,cum_sum_l,
            g_i,g_i_g_ik,mean_rank_dist,mean_cond_dist;
    double *data,*data_i,*out_s,*out_h,*out_m,*out_l;
    struct neighbor *emb_dist,*emb_dist_i,*emb_dist_ord,*emb_dist_ord_i,
            *knn,*knn_i,this_neighbor;
    
    const int n_samples=(int)mxGetM(prhs[0]);
    const int n_sensors=(int)mxGetN(prhs[0]);
    const int emb_d=(int)mxGetScalar(prhs[1]);
    const int tau=(int)mxGetScalar(prhs[2]);
    const int k=(int)mxGetScalar(prhs[3]);
    const int w_emb=(int)mxGetScalar(prhs[4]);
    const int states_eff_step=(int)mxGetScalar(prhs[5]);
    const int n_threads=omp_get_num_procs();
    
    n_states=n_samples-(tau*(emb_d-1));
    n_states_eff=n_states/states_eff_step;
    n_indexes=n_sensors*n_sensors;
    g_i=((float)n_states)/2;
    g_i_g_ik=g_i-(((float)k+1)/2);
    
    data=(double*)mxGetPr(prhs[0]);
    
    emb_dist=(struct neighbor*)mxMalloc(n_states_eff*n_states*n_sensors*sizeof(struct neighbor));
    emb_dist_ord=(struct neighbor*)mxMalloc(n_states_eff*n_states*n_sensors*sizeof(struct neighbor));
    mean_d=(float*)mxMalloc(n_states_eff*n_sensors*sizeof(float));
    mean_d_nn=(float*)mxMalloc(n_states_eff*n_sensors*sizeof(float));
    knn=(struct neighbor*)mxMalloc(n_states_eff*k*n_sensors*sizeof(struct neighbor));

    omp_set_num_threads(n_threads);
    
    for(sensor_i=0;sensor_i<n_sensors;sensor_i++)
    {
        data_i=data+sensor_i*n_samples;
        emb_dist_i=emb_dist+sensor_i*n_states_eff*n_states;
        emb_dist_ord_i=emb_dist_ord+sensor_i*n_states_eff*n_states;
        knn_i=knn+sensor_i*n_states_eff*k;
        
        #pragma omp parallel
        {
            #pragma omp for private(state_eff_i,state_i,state_j,tau_i,dist_diff,k_in,actual_state,this_neighbor) reduction(+:cum_dist_sum_1,dist_sum_1,dist_sum_2)
            for(state_eff_i=0;state_eff_i<n_states_eff;state_eff_i++)
            {
                state_j=state_eff_i*states_eff_step;
                cum_dist_sum_1=0.0;
                for(state_i=0;state_i<n_states;state_i++)
                {
                    dist_sum_1=0.0;
                    for(tau_i=0;tau_i<emb_d;tau_i++)
                    {
                        dist_diff=(float)data_i[state_i+tau_i*tau]-(float)data_i[state_j+tau_i*tau];
                        dist_sum_1+=dist_diff*dist_diff;
                    }
                    cum_dist_sum_1+=dist_sum_1;
                    emb_dist_i[state_i+state_eff_i*n_states].index=state_i;
                    emb_dist_ord_i[state_i+state_eff_i*n_states].index=state_i;
                    emb_dist_i[state_i+state_eff_i*n_states].dist=dist_sum_1;
                    emb_dist_ord_i[state_i+state_eff_i*n_states].dist=dist_sum_1;
                }
                qsort(emb_dist_ord_i+state_eff_i*n_states,n_states,sizeof(struct neighbor),comp);
                mean_d[state_eff_i+sensor_i*n_states_eff]=cum_dist_sum_1/((float)n_states-1);
                
                state_i=1;
                k_in=0;
                actual_state=emb_dist_ord_i[state_eff_i*n_states].index;
                dist_sum_2=0.0;
                while(k_in<k && state_i<n_states)
                {
                    this_neighbor=emb_dist_ord_i[state_i+state_eff_i*n_states];
                    if(abs(this_neighbor.index-actual_state)>w_emb)
                    {
                        knn_i[k_in+state_eff_i*k].index=this_neighbor.index;
                        knn_i[k_in+state_eff_i*k].dist=this_neighbor.dist;
                        dist_sum_2+=this_neighbor.dist;
                        k_in++;
                    }
                    state_i++;
                }
                mean_d_nn[state_eff_i+sensor_i*n_states_eff]=dist_sum_2/(float)k;
            }
        }
    }
    
    adj_s=(float*)mxMalloc(n_indexes*sizeof(float));
    adj_h=(float*)mxMalloc(n_indexes*sizeof(float));
    adj_m=(float*)mxMalloc(n_indexes*sizeof(float));
    adj_l=(float*)mxMalloc(n_indexes*sizeof(float));
    
    for(sensor_i=0;sensor_i<n_sensors;sensor_i++)
    {
        #pragma omp parallel
        {
            #pragma omp for private(sensor_j,state_i,state_eff_i,tau_i) reduction(+:cum_sum_s,cum_sum_h,cum_sum_m,cum_sum_l,mean_cond_dist,mean_rank_dist)
            for(sensor_j=0;sensor_j<n_sensors;sensor_j++)
            {
                cum_sum_s=0;
                cum_sum_h=0;
                cum_sum_m=0;
                cum_sum_l=0;
                for(state_eff_i=0;state_eff_i<n_states_eff;state_eff_i++)
                {
                    mean_cond_dist=0.0;
                    mean_rank_dist=0.0;
                    for(tau_i=0;tau_i<k;tau_i++)
                    {
                        mean_cond_dist+=emb_dist[knn[tau_i+state_eff_i*k+sensor_j*n_states_eff*k].index+
                                state_eff_i*n_states+sensor_i*n_states_eff*n_states].dist;
                        state_i=0;
                        while(knn[tau_i+state_eff_i*k+sensor_j*n_states_eff*k].index !=
                                emb_dist_ord[state_i+state_eff_i*n_states+sensor_i*n_states*n_states_eff].index)
                            state_i++;
                        mean_rank_dist+=(float)state_i;
                    }
                    mean_cond_dist/=(float)k;
                    cum_sum_s+=mean_d_nn[state_eff_i+sensor_i*n_states_eff]/mean_cond_dist;
                    cum_sum_h+=logf(mean_d[state_eff_i+sensor_i*n_states_eff]/mean_cond_dist);
                    cum_sum_m+=(mean_d[state_eff_i+sensor_i*n_states_eff]-mean_cond_dist)/
                            (mean_d[state_eff_i+sensor_i*n_states_eff]-mean_d_nn[state_eff_i+sensor_i*n_states_eff]);
                    cum_sum_l+=(g_i-(mean_rank_dist/(float)k));
                }

                adj_s[sensor_j+sensor_i*n_sensors]=cum_sum_s/(float)n_states_eff;
                adj_h[sensor_j+sensor_i*n_sensors]=cum_sum_h/(float)n_states_eff;
                adj_m[sensor_j+sensor_i*n_sensors]=cum_sum_m/(float)n_states_eff;
                adj_l[sensor_j+sensor_i*n_sensors]=cum_sum_l/((float)n_states_eff*g_i_g_ik);
            }
        }
    }
    
    mxFree(emb_dist);
    mxFree(emb_dist_ord);
    mxFree(mean_d);
    mxFree(mean_d_nn);
    mxFree(knn);
    
    plhs[0]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    plhs[3]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    
    #pragma omp parallel sections private(index_i)
    {
        #pragma omp section
        {
        out_s=(double*)mxGetPr(plhs[0]);
        for(index_i=0;index_i<n_indexes;index_i++)
            out_s[index_i]=(double)adj_s[index_i];
        }
        
        #pragma omp section
        {
            out_h=(double*)mxGetPr(plhs[1]);
            for(index_i=0;index_i<n_indexes;index_i++)
                out_h[index_i]=(double)adj_h[index_i];
        }
        
        #pragma omp section
        {
            out_m=(double*)mxGetPr(plhs[2]);
            for(index_i=0;index_i<n_indexes;index_i++)
                out_m[index_i]=(double)adj_m[index_i];
        }
        
        #pragma omp section
        {
            out_l=(double*)mxGetPr(plhs[3]);
            for(index_i=0;index_i<n_indexes;index_i++)
                out_l[index_i]=(double)adj_l[index_i];
        }
        
    }
    
    mxFree(adj_s);
    mxFree(adj_h);
    mxFree(adj_m);
    mxFree(adj_l);
    
}

