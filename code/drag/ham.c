#include<gsl/gsl_complex.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<complex.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<gsl/gsl_complex_math.h>

struct params{
	int Nx;
	double Dx;
	double a1x;
	double a1y;
	int Nt;
	double vx;
	double vy;
	double sig;
	double om;
	double m;
	double g11;
	double g22;
	double g;
	double a;
	double Dt;
};

double dlt(int x, int y){
	if(x==y) return 1;
	else return 0;
}

double Ax(complex double *psi, int x, int y, struct params p){
	return p.a1x*conj(*(psi+p.Nx*x+y))**(psi+p.Nx*x+y);
}
double Ay(complex double *psi, int x, int y, struct params p){
	return p.a1y*conj(*(psi+p.Nx*x+y))**(psi+p.Nx*x+y);
}

complex double jx(complex double *psi, int x, int y, struct params p){
	return
	-I*(*(psi+p.Nx*x+y)*(conj(*(psi+p.Nx*(x+1)+y)-*(psi+p.Nx*(x-1)+y))/p.Dx/2+I*Ax(psi,x,y,p))*conj(*(psi+p.Nx+y))
	-conj(*(psi+p.Nx*x+y))*((*(psi+p.Nx*(x+1)+y)-*(psi+p.Nx*(x-1)+y))/p.Dx/2-I*Ax(psi,x,y,p))**(psi+p.Nx+y))/2/p.m;
}

complex double jy(complex double *psi, int x, int y, struct params p){
	return
	-I*(*(psi+p.Nx*x+y)*(conj(*(psi+p.Nx*x+y+1)-*(psi+p.Nx*x+y-1))/p.Dx/2+I*Ay(psi,x,y,p))*conj(*(psi+p.Nx+y))
	-conj(*(psi+p.Nx*x+y))*((*(psi+p.Nx*x+y+1)-*(psi+p.Nx*x+y-1))/p.Dx/2-I*Ay(psi,x,y,p))**(psi+p.Nx+y))/2/p.m;
}

complex double ihm(int i, int j, complex double *psi, double t, struct params p){

	return
	-I/p.Dx/p.Dx*(dlt(j,i+p.Nx)+dlt(j,i-p.Nx)-2*dlt(i,j))-I/p.Dx/p.Dx*(dlt(i,j+1)+dlt(i,j-1)-2*dlt(i,j))/2/p.m
	-Ax(psi,i,j,p)/2/p.Dx*(dlt(j,i+p.Nx)-dlt(j,i-p.Nx))-Ay(psi,i,j,p)/2/p.Dx*(dlt(j,i+1)-dlt(j,i-1))
	-dlt(i,j)*((Ax(psi,i+1,j,p)-Ax(psi,i-1,j,p))/2/p.Dx+(Ay(psi,i,j+1,p)-Ay(psi,i,j-1,p))/2/p.Dx
	+Ax(psi,i+1,j,p)*Ax(psi,i+1,j,p)+Ay(psi,i+1,j,p)*Ay(psi,i+1,j,p)) /* gauge covariant kinetic term*/
	
	+dlt(i,j)*(p.a1x*jx(psi,i,j,p)+p.a1y*jy(psi,i,j,p) /* current nonlinearity */

	+(pow(p.a1x*8*p.om/(p.g11-p.g22),2)+pow(p.a1y*8*p.om/(p.g11-p.g22),2))/2/p.m /* scalar field */

	+p.g*conj(*(psi+p.Nx*i+j))**(psi+p.Nx*i+j) /* density nonlinearity */

	+p.a*(cos(-(i*p.Dx-p.vx*t)*(i*p.Dx-p.vx*t)/p.sig/p.sig/2)+I*sin(-(i*p.Dx-p.vx*t)*(i*p.Dx-p.vx*t)/p.sig/p.sig/2)
	*cos(-(j*p.Dx-p.vy*t)*(j*p.Dx-p.vy*t)/p.sig/p.sig/2)+I*sin(-(j*p.Dx-p.vy*t)*(j*p.Dx-p.vy*t)/p.sig/p.sig/2))); /* impurity */

}


void ham(complex double * psi, double t, struct params p){

	gsl_vector_complex *X = gsl_vector_complex_alloc(p.Nx*p.Nx);
	gsl_vector_complex *b = gsl_vector_complex_alloc(p.Nx*p.Nx);
	gsl_matrix_complex *hm = gsl_matrix_complex_alloc(p.Nx*p.Nx,p.Nx*p.Nx);
	gsl_matrix_complex *hmp = gsl_matrix_complex_alloc(p.Nx*p.Nx,p.Nx*p.Nx);


	for(int i=0;i<p.Nx*p.Nx;i++){
	for(int j=0;j<p.Nx*p.Nx;j++){
		gsl_matrix_complex_set(hm,i,j,
		gsl_complex_rect(creal(ihm(i,j,psi,t,p)),cimag(ihm(i,j,psi,t,p))));

		gsl_matrix_complex_set(hmp,i,j,
		gsl_complex_rect(creal(dlt(i,j)+ihm(i,j,psi,t,p)),cimag(ihm(i,j,psi,t,p))));
	}}

	for(int i=0;i<p.Nx;i++){
	for(int j=0;j<p.Nx;j++){
		gsl_vector_complex_set(b,p.Nx*i+j,gsl_complex_rect(creal(*(psi+p.Nx*i+j)),cimag(*(psi+p.Nx*i+j))));
	}}

	gsl_blas_cgemv(CblasNoTrans,gsl_complex_rect(-1/2*p.Dt,0),b,gsl_complex_rect(1,0),b);

	gsl_permutation *p = gsl_permutation_alloc(p.Nx*p.Nx,p.Nx*p.Nx);
	
	int sgn;
	
	gsl_linalg_complex_LU_decomp(hmp,p,&sgn);
	gsl_linalg_complex_LU_solve(hmp,p,b,x);

	
	for(int i=0;i<p.Nx;i++){
	for(int j=0;j<p.Nx;j++){
		*(psi+p.Nx*i+j)=gsl_vector_get(x,p.Nx*i+j)
	}}

	gsl_vector_complex_free(X);
	gsl_vector_complex_free(b);
	gsl_matrix_complex_free(hm);
	gsl_matrix_complex_free(hmp);
}
