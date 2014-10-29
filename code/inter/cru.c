#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_eigen.h>
#include<math.h>
#include<complex.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>

#define PI 3.141592

int main(){
	
	FILE *fi,*fil;
	fi=fopen("state.dat","w");
	fil=fopen("energy","w");
	

	int st =0;/* print which state */
	int L = 80;
	double t=-1;

	double phase[L/4];

/*for(int k=0;k<=30;k++){*/
	
	for(int i=0;i<L/8;i++){
		phase[i]=PI/2;
		}


	for(int i=L/8;i<L/4;i++){
		phase[i]=-PI/2;
		}

/*	phase[10]=-PI/2;
	phase[11]=-PI/2;
	phase[12]=-PI/2;
	phase[13]=-PI/2;
	phase[14]=-PI/2;
	phase[15]=-PI/2; */

	gsl_matrix_complex *H=gsl_matrix_complex_alloc(L,L);
	gsl_vector *eval=gsl_vector_alloc(L);
	gsl_matrix_complex *evec=gsl_matrix_complex_alloc(L,L);
	gsl_eigen_hermv_workspace *W = gsl_eigen_hermv_alloc(L);
	

	for(int i=0;i<L-2;i++){

		int j=i/4; 

		if(i%2==1){
			gsl_matrix_complex_set(H,i,i+1,gsl_complex_polar(t,0));
			gsl_matrix_complex_set(H,i+1,i,gsl_complex_polar(t,0));
			gsl_matrix_complex_set(H,i,i+2,gsl_complex_polar(t,phase[j]));
			gsl_matrix_complex_set(H,i+2,i,gsl_complex_polar(t,-phase[j]));
		}
		if(i%2==0){
			gsl_matrix_complex_set(H,i,i+3,gsl_complex_polar(t,0));
			gsl_matrix_complex_set(H,i+3,i,gsl_complex_polar(t,0));
			gsl_matrix_complex_set(H,i,i+2,gsl_complex_polar(t,-phase[j]));
			gsl_matrix_complex_set(H,i+2,i,gsl_complex_polar(t,phase[j]));
		}
		
	}


	gsl_eigen_hermv(H,eval,evec,W);
	gsl_eigen_hermv_sort(eval,evec,GSL_EIGEN_SORT_VAL_ASC);


	for(int i=0;i<L;i++){
		fprintf(fi,"%i	%.20g	%.20g\n",i,GSL_REAL(gsl_matrix_complex_get(evec,i,st)),GSL_IMAG(gsl_matrix_complex_get(evec,i,st)));
	}
		fprintf(fi,"\n\n");
	for(int i=0;i<L;i++){
		fprintf(fil,"%.20g\n",gsl_vector_get(eval,i));
	}	

	return 0;
}
