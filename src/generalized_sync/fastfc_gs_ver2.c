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
#include <matrix.h>
#include <math.h>
#include <omp.h>

#define INF (float)mxGetInf();

typedef struct neighbour{
    int index;
    float dist;
}neighbour;

int cmp(const void *a,const void *b)
{
    return (((neighbour*)a)->dist < ((neighbour*)b)->dist)?
        -1:(((neighbour*)a)->dist > ((neighbour*)b)->dist);
}

void init_knn(const int n_sensors,const int n_states,const int k,
        neighbour * k_nn)
{
    int sensor_i,state_i,k_i;
    for(sensor_i=0;sensor_i<n_sensors;sensor_i++)
        for(state_i=0;state_i<n_states;state_i++)
            for(k_i=0;k_i<k;k_i++)
                k_nn[sensor_i*n_states*k+state_i*k+k_i].dist=(float)0.;
}

void push2knn(int * k_in,const int k,const int i,const int W,
        const float d,const int j,neighbour * k_nn)
{
    int k_ii,k_i=0;
    if(abs(i-j)>W)
    {
        if(*k_in<k)
        {
            while((k_nn[k_i].dist>d)&(abs(i-j)>W))
                k_i++;
            for(k_ii=*k_in;k_ii>k_i;k_ii--)
                k_nn[k_ii]=k_nn[k_ii-1];
            
            k_nn[k_i].dist=d;
            k_nn[k_i].index=j;
            *k_in+=1;
        } else if(k_nn[k_i].dist>d) {
            k_i++;
            while((k_i<k)&(k_nn[k_i].dist>d))
                k_i++;
            for(k_ii=0;k_ii<k_i-1;k_ii++)
                k_nn[k_ii]=k_nn[k_ii+1];
            k_nn[k_i-1].dist=d;
            k_nn[k_i-1].index=j;
        }
    }
}

void build_emb_dist(const int m,const int tau,const int n_states,
        const int k,const int W,
        neighbour * emb_d,
        float * r,float * r_k,
        neighbour * k_nn,const double * X)
{
    int state_i,state_j,k_in,d,k_i;
    float dist,cum_dist,cum_dist_k;
    
    k_in=0;
    
    for(state_i=0;state_i<n_states;state_i++)
    {
        cum_dist=0;
        for(state_j=0;state_j<state_i;state_j++)
        {
            dist=0;
            for(d=0;d<m;d++)
                dist+=((float)X[state_i+d*tau]-(float)X[state_j+d*tau])*
                        ((float)X[state_i+d*tau]-(float)X[state_j+d*tau]);
            
            cum_dist+=dist;
            emb_d[state_i*n_states+state_j].index=state_j;
            emb_d[state_i*n_states+state_j].dist=dist;
            
            push2knn(&k_in,k,state_i,W,dist,state_j,k_nn+state_i*k);
        }
        emb_d[state_i*n_states+state_j].index=state_j;
        emb_d[state_i*n_states+state_j].dist=INF;
        for(state_j++;state_j<n_states;state_j++)
        {
            dist=0;
            for(d=0;d<m;d++)
                dist+=((float)X[state_i+d*tau]-(float)X[state_j+d*tau])*
                        ((float)X[state_i+d*tau]-(float)X[state_j+d*tau]);
            
            cum_dist+=dist;
            emb_d[state_i*n_states+state_j].index=state_j;
            emb_d[state_i*n_states+state_j].dist=dist;
            
            push2knn(&k_in,k,state_i,W,dist,state_j,k_nn+state_i*k);
        }
		k_in=0;
        r[state_i]=cum_dist/(float)(n_states-1);
        
        cum_dist_k=0.;
        for(k_i=0;k_i<k;k_i++)
            cum_dist_k+=k_nn[state_i*k+k_i].dist;
        r_k[state_i]=cum_dist_k/(float)k;
    }
}

void build_SHM_indexes(const int sens_i,const int sens_j,const int n_sensors,
        const int n_states,const int k,
        neighbour * emb_d,
        neighbour * knn,
        float * r_ik,
        float * r_i,
        float * S,float * H, float * M)
{
    int state_i,k_i;
    float dist,cum_S,cum_H,cum_M;
    
    cum_S=0.;
    cum_H=0.;
    cum_M=0.;

    for(state_i=0;state_i<n_states;state_i++)
    {
        dist=0.;
        for(k_i=0;k_i<k;k_i++)
            dist+=emb_d[state_i*n_states+knn[state_i*k+k_i].index].dist;
        dist/=(float)k;
        cum_S+=r_ik[state_i]/dist;
        cum_H+=logf(r_i[state_i]/dist);
        cum_M+=(r_i[state_i]-dist)/(r_i[state_i]-r_ik[state_i]);
    }
    
    S[sens_i*n_sensors+sens_j]=cum_S/(float)n_states;
    H[sens_i*n_sensors+sens_j]=cum_H/(float)n_states;
    M[sens_i*n_sensors+sens_j]=cum_M/(float)n_states;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    int sens_i,sens_j,index_i;
    float *r_i,*r_ik,*S,*H,*M;
    double *X_d,*S_d,*H_d,*M_d;
    neighbour *emb_dist,*k_nn;
    
    const int n_samples=(int)mxGetM(prhs[0]);
    const int n_sensors=(int)mxGetN(prhs[0]);
    const int k=(int)mxGetScalar(prhs[1]);
    const int m=(int)mxGetScalar(prhs[2]);
    const int tau=(int)mxGetScalar(prhs[3]);
    const int W=(int)mxGetScalar(prhs[4]);
    const int n_states=n_samples-(m-1)*tau;
    const int n_indexes=n_sensors*n_sensors;
    
    X_d=(double*)mxGetPr(prhs[0]);
    
    emb_dist=(neighbour*)mxMalloc(n_states*n_states*n_sensors*sizeof(neighbour));
    k_nn=(neighbour*)mxMalloc(n_states*k*n_sensors*sizeof(neighbour));
    r_i=(float*)mxMalloc(n_states*n_sensors*sizeof(float));
    r_ik=(float*)mxMalloc(n_states*n_sensors*sizeof(float));
    
    init_knn(n_sensors,n_states,k,k_nn);

    omp_set_num_threads(omp_get_num_procs());
    #pragma omp parallel private(sens_i) shared(m,tau,n_states,k,W)
    {
        #pragma omp for
        for(sens_i=0;sens_i<n_sensors;sens_i++)
            build_emb_dist(m,tau,n_states,k,W,
                    emb_dist+sens_i*n_states*n_states,
                    r_i+sens_i*n_states,
                    r_ik+sens_i*n_states,
                    k_nn+sens_i*n_states*k,
                    X_d+sens_i*n_samples);
    }

    S=(float*)mxMalloc(n_indexes*sizeof(float));
    H=(float*)mxMalloc(n_indexes*sizeof(float));
    M=(float*)mxMalloc(n_indexes*sizeof(float));

    
    for(sens_i=0;sens_i<n_sensors;sens_i++)
    {
        #pragma omp parallel private(sens_j) shared(sens_i,n_sensors,n_states,k)
        {
        #pragma omp for
        for(sens_j=0;sens_j<n_sensors;sens_j++)
        {
            build_SHM_indexes(sens_i,sens_j,n_sensors,n_states,k,
                    emb_dist+sens_i*n_states*n_states,
                    k_nn+sens_j*k*n_states,
                    r_ik+sens_i*n_states,
                    r_i+sens_i*n_states,
                    S,H,M);
        }
        }
    }

    mxFree(emb_dist);
    mxFree(k_nn);
    mxFree(r_i);
    mxFree(r_ik);
    
    plhs[0]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(n_sensors,n_sensors,mxREAL);

    #pragma omp parallel sections private(index_i)
    {
        #pragma omp section
        {
        S_d=(double*)mxGetPr(plhs[0]);
        for(index_i=0;index_i<n_indexes;index_i++)
            S_d[index_i]=(double)S[index_i];
        }
        
        #pragma omp section
        {
            H_d=(double*)mxGetPr(plhs[1]);
            for(index_i=0;index_i<n_indexes;index_i++)
                H_d[index_i]=(double)H[index_i];
        }
        
        #pragma omp section
        {
            M_d=(double*)mxGetPr(plhs[2]);
            for(index_i=0;index_i<n_indexes;index_i++)
                M_d[index_i]=(double)M[index_i];
        }
        
    }
    
    mxFree(S);
    mxFree(H);
    mxFree(M);
    
}
