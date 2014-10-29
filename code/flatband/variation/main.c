#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_eigen.h>

void mset(int N, gsl_matrix *M, gsl_matrix *Mb);
void msolve(gsl_matrix *M, gsl_matrix *Mb,gsl_vector_complex *al, gsl_vector *be, gsl_matrix_complex *evec,gsl_eigen_genv_workspace *W);
void print(int p,int N, gsl_vector_complex *al, gsl_vector *be, gsl_matrix_complex *evec);

int main(){




int N=2; /*number of particles*/
int p=1; /* 1: print coefficients and eigenvalues    0: print eigenvales only*/




int L=4*N-2;

gsl_matrix *M= gsl_matrix_alloc(L/2,L/2), *Mb= gsl_matrix_alloc(L/2,L/2);
gsl_vector_complex *al = gsl_vector_complex_alloc(L/2); gsl_vector *be = gsl_vector_alloc(L/2);
gsl_matrix_complex *evec = gsl_matrix_complex_alloc(L/2,L/2);
gsl_eigen_genv_workspace *W=gsl_eigen_genv_alloc(L/2);


mset(L,M,Mb);

msolve(M,Mb,al,be,evec,W);

print(p,N,al,be,evec);

gsl_matrix_free(M); gsl_matrix_free(Mb); gsl_vector_free(be); gsl_vector_complex_free(al); gsl_matrix_complex_free(evec); gsl_eigen_genv_free(W);
 
return 0;
}

