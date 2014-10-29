#include<stdio.h>
#include<complex.h>
#include<stdlib.h>

void ham(complex double *psi, double t, struct params p);

struct params p;

int main(){

	p.Nx = 1000; /*Number of sites in 1D  */
	p.Dx = 0.0001; /* Position space step size */
	p.a1x = 1;
	p.a1y = 1;
	p.Nt = 1000;
	p.vx = 1; /*impurity velocity */
	p.vy = 0;
	p.sig = 1; /* impurity width */
	p.om = 1; /* rabi freq: characterizes scalar term strength */
	p.m = 1; /* mass */
	p.g11 = 1;
	p.g22 = 1;
	p.g = 1; /* strength of density nonlinearity */
	p.a = 1; /*impurity strength */
	p.Dt = 0.001 /* time step */

	complex double *psi = malloc(p.Nx*p.Nx*sizeof(complex double));
	
	for(int i=0;i<p.Nx*p.Nx;i++){
		*(psi+i)=1./p.Nx/p.Dx;
	}

	for(double i=0; i<p.Nt;i++){
		double t=p.Dx*p.Nx/p.vx*1./Nt*i;
		ham(psi,t,p);
		for(int j=0;j<p.Nx;j++){
		for(int k=0;k<p.Nx;k++){
			printf("%g\n",pow(creal(*(psi+p.Nx*j+k),2))+pow(cimag(*(psi+p.Nx*j+k),2)));
		}}
		printf("\n\n");
	}

	free(psi);

	return 0;
}
