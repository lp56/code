#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_complex.h>
int index(int *x, int *y, int N);
void xhop_l1(int i, gsl_matrix *M, int N, double H[2]);
void xhop_l2(int i, gsl_matrix *M, int N, double H[2], int s);
void xhop_r1(int i, gsl_matrix *M, int N, double H[2]);
void xhop_r2(int i, gsl_matrix *M, int N, double H[2], int s);
void yhop_l1(int i, gsl_matrix *M, int N, double H[2]);
void yhop_l2(int i, gsl_matrix *M, int N, double H[2], int s);
void yhop_r1(int i, gsl_matrix *M, int N, double H[2]);
void yhop_r2(int i, gsl_matrix *M, int N, double H[2], int s);
void onsite(int i, gsl_matrix *M, int N, double U);
void set(gsl_matrix *M, double H[2], int N, int S,double t, int s, double U);

void alpha(double U, double H[2], double t, double alp[21], double *ec){

printf("running alpha\n");

int N=7;
int S=N*(N+1)/2;
int s=0;

gsl_eigen_symmv_workspace *W=gsl_eigen_symmv_alloc(S);
gsl_vector *evals = gsl_vector_alloc(S);
gsl_matrix *evecs = gsl_matrix_alloc(S,S);
gsl_matrix *M=gsl_matrix_alloc(S,S);
set(M,H,N,S,t,s,U);
gsl_eigen_symmv(M,evals,evecs,W);
gsl_eigen_symmv_sort(evals,evecs,GSL_EIGEN_SORT_VAL_ASC);




for (int i=0; i<S; i++){
		alp[i]=gsl_matrix_get(evecs,i,0);					
		}

*ec=gsl_vector_get(evals,0);


printf("alpha ran ok\n");

}

