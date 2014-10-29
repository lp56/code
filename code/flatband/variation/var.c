#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<complex.h>
#include<math.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_permutation.h>
#include<stdlib.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_eigen.h>

int delta(int L, int i, int j){
	i=i%L; j=j%L;
	if(i==j) return 1;
	else return 0;
}	

int step(int i){
	if(i>=0) return 1;
	else return 0;
}

double n(int L,double a, double b,  double c, int i, int j){
	return step(abs(i%L-j%L)-1)*(a*sqrt(2)-b)*(c*sqrt(2)-b)+(delta(L,j,i+1)+delta(L,j,i-1))*a*c+
	delta(L,i,j)*(a*a+b*b+c*c);
}


void mset(int L, gsl_matrix *M, gsl_matrix *Mb){
	double a=sqrt(2)*(sqrt(17)-3)/(sqrt(17)-5),c=1, A=a*sqrt(2)+1, B=sqrt(2)*2, C=1+sqrt(2)*a;
	for(int i=0;i<L/2;i++){
		for(int j=0;j<L/2;j++){
			gsl_matrix_set(M,i,j,n(L,1,a,c,i,j));
			gsl_matrix_set(Mb,i,j,n(L,A,B,C,i,j));
		}
	}
}

void msolve( gsl_matrix *M, gsl_matrix *Mb,gsl_vector_complex *al, gsl_vector *be, gsl_matrix_complex *evec, gsl_eigen_genv_workspace *W){
	
	gsl_eigen_genv(Mb,M,al,be,evec,W);
}

void print(int p,int N, gsl_vector_complex *al, gsl_vector *be, gsl_matrix_complex *evec){
	
	FILE *t;
	char *m = "w";
	t=fopen("en",m);

	int L=4*N-2;
	for(int k=0;k<L/2;k++){
		double eva= GSL_REAL(gsl_vector_complex_get(al,k))/gsl_vector_get(be,k);
		fprintf(t,"E-E0 (%i)= %g\n\n",k,eva);
		if(p==1){
			for(int i=0;i<L/2;i++){
				double evr = GSL_REAL(gsl_matrix_complex_get(evec,i,k)),
				evi = GSL_IMAG(gsl_matrix_complex_get(evec,i,k));
				fprintf(t,"%g+i%g\n",evr,evi);
		}
			}
		fprintf(t,"\n");
	}	
}
