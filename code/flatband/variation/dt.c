#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_permutation.h>

int main(){

int *signum;
int si=-1;

gsl_permutation *p=gsl_permutation_calloc(2);

gsl_matrix *M= gsl_matrix_alloc(2,2);

gsl_matrix_set(M,0,0,2);

gsl_matrix_set(M,0,1,2);

gsl_matrix_set(M,1,0,3);

gsl_matrix_set(M,1,1,1);

gsl_linalg_LU_decomp(M,p,signum);

printf("%g\n",gsl_linalg_LU_det(M,si));

gsl_matrix_free(M);
gsl_permutation_free(p);

return 0;
}


