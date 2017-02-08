#include <mex.h>
#include <math.h>
// 
// __inline double fastabs(double x){
//     unsigned long int *x_int=&x;
//     unsigned long int n=*x_int&0x7fffffff;
//     return (double)n;
//   }

// 
// inline double fastabs(double x){
//     float a =(float)x;
// 
//     unsigned int b = *reinterpret_cast<unsigned int*>(&a);
//     return (double) (b&0x7fffffff) ;
// }
// 

//  inline double fastabs(double x){
//         float x_f=(float)x;
//         unsigned int=
//         return (double&) (((unsigned int&)x_f)&0x7fffffff) ;
//    }
// 

// double fastabs(double x){
//     double x_abs=fabs(x);
//     double *ptr_x;
//     ptr_x=(double*)mxMalloc(2*sizeof(double));
//     *ptr_x=x;
//     *(ptr_x+1)=0x7fffffffffffffff;
//     unsigned long int x_int=(((unsigned long int&)x)&0x7fff ffff ffff ffff);
//     //return (double&) x_int ;
// 	return *(ptr_x+1);
// }

// double fastabs(double x)
// {   
//     double mask=9.2233720368547758e18;
//     unsigned long int xx=((unsigned long int&)x&(unsigned long int&)mask);
//     return (double&)xx;
// }

double fastabs(float x)
{
    unsigned int int_x=((unsigned int&)x&0x7fffffff);
    return (double)(float&)int_x;
}
    



void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
//     unsigned int int_in;
    double x;
    double in=mxGetScalar(prhs[0]);

    x=fastabs(in);
    
    plhs[0]=mxCreateDoubleScalar(x);

}
    
    
    
//  
//     int_in=*(long int*)&in;
//     p_in=(long int*)&in;
// //     reinterpret_cast
// //  int_in_1=*p_in;
// //  int_in_2=3;
// //  
// //  int_in_1=int_in_1*int_in_2;
// 	mexPrintf("\nsize of int: %d\n",sizeof(int));
// 	mexPrintf("\nsize of double: %d\n",sizeof(double));
// 	mexPrintf("\nsize of long int: %d\n",sizeof(long int));
// 	mexPrintf("\nsize of long long int: %d\n",sizeof(long long int));
// 	mexPrintf("\nsize of float: %d\n",sizeof(float));
// 	mexPrintf("\nsize of pointer to int: %d\n",sizeof(int*));
// 	mexPrintf("\nsize of pointer to double: %d\n",sizeof(double));
// 	mexPrintf("\nsize of pointer to long int: %d\n",sizeof(long int*));
// 	mexPrintf("\nsize of pointer to float: %d\n",sizeof(float*));
// 

// 	 mexPrintf("\nint_in: %f\n",int_in);
// 	 mexPrintf("\np_in: %f\n\n",*p_in);
 
                                      
