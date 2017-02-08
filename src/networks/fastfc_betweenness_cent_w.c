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
#include <matrix.h>

#define INF (float)mxGetInf();

typedef struct _node{
    int inode;
    int val;
    struct _node* prev;
    struct _node* next;
}node;

node * flush_node(int n,int i,node * list,node * list_init)
{
    if(i==list->val)
        return list->next;
    if(i==n-1)
        return list;
    else{   
        list_init[i].prev->next=list_init[i].next;
        list_init[i].next->prev=list_init[i].prev;
        return list;
	}
}

node * init_nodes(int n,int inode,node * list)
{
    int i=0;
    list->inode=inode;
    list->val=i;
    list->prev=NULL;
    list->next=&list[1];
    for(i=1;i<n-1;i++)
    {   
        list[i].inode=inode;
        list[i].val=i;
        list[i].prev=&list[i-1];
        list[i].next=&list[i+1];
    }
    list[n-1].inode=inode;
    list[n-1].val=n-1;
    list[n-1].prev=&list[n-2];
    list[n-1].next=NULL;
	return list;
}

void init_path(int n,int *p)
{
    int i;
    for(i=0;i<n;i++)
        p[i*n]=0;
}

void push_node2path(int n,int prev_n,int next_n,int* p)
{
    p[prev_n*n+(p[prev_n*n]++)+1]=next_n;
}

void count_path(int n, int n_i, int * C, int * p, int * list)
{
    int l_list, i, k, i_old, nodes2find,nodes_found;
    
    list[0]=n_i;
    l_list=1;
    nodes_found=1;
    nodes2find=n-nodes_found;
    i=n_i;
    
    while(nodes2find)
    {
        if(p[i*n])
        {
            list[l_list++]=p[i*n+p[i*n]];
             i_old=i;
             nodes2find--;
             i=p[i*n+p[i*n]];
        }else{
            for(k=1;k<l_list-1;k++)
                C[list[k]]+=1;
            
            p[i_old*n]--;
            l_list=1;
            nodes_found++;
            nodes2find=n-nodes_found;
            i=n_i;
        }
    }
    for(k=1;k<l_list-1;k++)
        C[list[k]]+=l_list-1-k;
}
        
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    int **itinerary,**path,*Betw,*B,*B_i,i,n_i,inode,min_node,l_nodes,min_node_prev,th_i;
    float *D,*D_i,min_d;
    double *L,*D_d,*B_d,*Betw_d;
    node **nodes_init,*nodes,*node_i;
    
    const int n=(int)mxGetM(prhs[0]);
    const int n_elements=n*n;
    const int n_threads=omp_get_num_procs();
    
    nodes_init=(node**)mxMalloc(n_threads*sizeof(node*));
    path=(int**)mxMalloc(n_threads*sizeof(int*));
    itinerary=(int**)mxMalloc(n_threads*sizeof(int*));
    for(th_i=0;th_i<n_threads;th_i++)
    {
        nodes_init[th_i]=(node*)mxMalloc(n*sizeof(node));
        path[th_i]=(int*)mxMalloc(n_elements*sizeof(int));
        itinerary[th_i]=(int*)mxMalloc(n*sizeof(int));
    }
    
    L=(double*)mxGetPr(prhs[0]);
    D=(float*)mxMalloc(n_elements*sizeof(float));
    B=(int*)mxCalloc(n_elements,sizeof(int));
    Betw=(int*)mxCalloc(n,sizeof(int));
    
    for(i=0;i<n_elements;i++)
        D[i]=INF;
    for(i=0;i<n;i++)
        D[i*n+i]=0.;
    
    omp_set_num_threads(n_threads);
    #pragma omp parallel private(n_i,i,l_nodes,inode,th_i,node_i,nodes,D_i,B_i,min_d,min_node,min_node_prev) shared(n,path,nodes_init,Betw,itinerary)
    {
        #pragma omp for 
        for(n_i=0;n_i<n;n_i++)
        {
            l_nodes=n-1;
            inode=n_i;
            th_i=omp_get_thread_num();
            init_path(n,path[th_i]);
            nodes=init_nodes(n,inode,nodes_init[th_i]); 
            nodes=flush_node(n,inode,nodes,nodes_init[th_i]);

            D_i=D+n_i*n;
            B_i=B+n_i*n;
            while(l_nodes)
            {
                node_i=nodes;
                for(i=0;i<l_nodes;i++)
                {
                    if(D_i[inode]+(float)L[node_i->val*n+inode]<D_i[node_i->val])
                    {
                        D_i[node_i->val]=D_i[inode]+(float)L[node_i->val*n+inode];
                        B_i[node_i->val]=B_i[inode]+1;
                        node_i->inode=inode;
                    }
                    node_i=node_i->next;
                }
                node_i=nodes;
                min_d=D_i[node_i->val];
                min_node=node_i->val;
                min_node_prev=node_i->inode;
                for(i=0;i<l_nodes;i++)
                {
                    if(D_i[node_i->val]<min_d)
                    {
                        min_d=D_i[node_i->val];
                        min_node=node_i->val;
                        min_node_prev=node_i->inode;
                    }
                    node_i=node_i->next;
                }
                inode=min_node;
                nodes=flush_node(n,min_node,nodes,nodes_init[th_i]);
                push_node2path(n,min_node_prev,inode,path[th_i]);
                l_nodes--;
            }
            
            count_path(n,n_i,Betw,path[th_i],itinerary[th_i]);
        }
    }

    plhs[0]=mxCreateDoubleMatrix(n,n,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(n,n,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(n,1,mxREAL);
        
    #pragma omp parallel sections private(i)
    {
        #pragma omp section
        {
            D_d=(double*)mxGetPr(plhs[0]);
            for(i=0;i<n_elements;++i)
                D_d[i]=(double)D[i];
        }

        #pragma omp section
        {
            B_d=(double*)mxGetPr(plhs[1]);
            for(i=0;i<n_elements;++i)
                B_d[i]=(double)B[i];
        }

        #pragma omp section
        {
            Betw_d=(double*)mxGetPr(plhs[2]);
            for(i=0;i<n;++i)
                Betw_d[i]=(double)Betw[i];
        }

    }
    
    for(th_i=0;th_i<n_threads;th_i++)
    {
        mxFree(nodes_init[th_i]);
        mxFree(path[th_i]);
        mxFree(itinerary[th_i]);
    }
    mxFree(nodes_init);
    mxFree(path);
    mxFree(itinerary);
    mxFree(D);
    mxFree(B);
    mxFree(Betw);

}
