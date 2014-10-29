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
void alpha(double U, double H[2], double t, double alp[21], double *ec);
int delta(int i, int j, int N);
int mod(int i);
double jev(double alp[21], int N, int h);
double jevr(double alp[21], int N, int h);
double jevl(double alp[21], int N, int h);
void m1set(int n, double U, double H[2], double t, gsl_matrix *M, double *ec);
void m2set(int n, double U, double H[2], double t, gsl_matrix *M, double *ec);

int main(){

int N = 5; /*number of particles*/
double U = 10000;
double H[] = {sqrt(2),1};
double t = sqrt(2);
/*double ev = -2;*/
double a=0; double *ec=&a;

int n=2*(N-2)+3;

gsl_matrix *m1= gsl_matrix_alloc(n,n);
gsl_matrix *m2= gsl_matrix_alloc(n,n);

	double alp[21];

	alpha(U,H,t,alp,ec);

	int h=0;	
	for(int i=0; i<n;i++){
		for(int j=0; j<n; j++){
			h=j-i;
			if(h<0) h+=n;
			gsl_matrix_set(m1,i,j,jev(alp,n,h));
	/*		printf("%8g", gsl_matrix_get(m1,i,j));*/
		}
/*	printf("\n");*/}


	alpha(U,H,t,alp,ec);

	h=0;
	for(int i=0; i<n;i++){
		for(int j=0; j<n; j++){
			h=j-i;
			if(h<0) h+=n;
			gsl_matrix_set(m2,i,j,jevr(alp,n,h)+jevl(alp,n,h));
	/*		printf("%10g", gsl_matrix_get(m2,i,j));*/
		}
/*	printf("\n")*/;
	}


gsl_vector_complex *al = gsl_vector_complex_alloc(n);
gsl_vector *be = gsl_vector_alloc(n);
gsl_matrix_complex *evec = gsl_matrix_complex_alloc(n,n);
gsl_eigen_genv_workspace *w = gsl_eigen_genv_alloc(n);
gsl_eigen_genv(m1,m2,al,be,evec,w);

FILE *tea;
char *m="w";
tea=fopen("egg",m);



for(int i=0;i<n;i++){
	fprintf(tea,"%g\n",GSL_REAL(gsl_vector_complex_get(al,i))/gsl_vector_get(be,i));
}

fprintf(tea,"\n%g\n",*ec);

gsl_eigen_genv_free(w);
gsl_matrix_free(m1);
gsl_matrix_free(m2);
gsl_vector_complex_free(al);
gsl_vector_free(be);
gsl_matrix_complex_free(evec);
return 0;
}
