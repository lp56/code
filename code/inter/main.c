#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_eigen.h>
#include<math.h>
#include<complex.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
int main(){
	
	FILE *fi,*fil;
	fi=fopen("state","w");
	fil=fopen("energy","w");
	

	int st = 127;/* print which state */
	int L_1 = 100, L_c = 50, L_2 = 100; /* L_c = length of crutez ladder */ int L=L_1+L_2+L_c;
	gsl_complex t1= gsl_complex_rect(-1,0);
	gsl_complex t2= gsl_complex_rect(-1,0); /* hoppings */
	gsl_complex tc= gsl_complex_rect(-1,0);
	gsl_complex ga11 = gsl_complex_rect(-0.1,0);
	gsl_complex ga12 = gsl_complex_rect(-0.1,0);
	gsl_complex ga21 = gsl_complex_rect(-0.1,0);/*couplings between chains and cl*/
	gsl_complex ga22 = gsl_complex_rect(-0.1,0);
	gsl_complex la = gsl_complex_rect(0,0); /* energy offset */
	int pbc = 1; /* 1 = periodic, 0 = open. NOTE, FOR PBC = 1 MUST HAVE t1=t2. */

	gsl_matrix_complex *H=gsl_matrix_complex_alloc(L,L);
	gsl_vector *eval=gsl_vector_alloc(L);
	gsl_matrix_complex *evec=gsl_matrix_complex_alloc(L,L);
	gsl_eigen_hermv_workspace *W = gsl_eigen_hermv_alloc(L);

	for(int i=0;i<L_1;i++){
		for(int j=0;j<L_1;j++){
			if(i==j+1||i==j-1) gsl_matrix_complex_set(H,i,j,t1);
			if(i==j) gsl_matrix_complex_set(H,i,j,la);
		}
	}

	if(L_1!=0){
	gsl_matrix_complex_set(H,L_1-1,L_1,ga11);
	gsl_matrix_complex_set(H,L_1,L_1-1,ga11);
	gsl_matrix_complex_set(H,L_1-1,L_1+1,ga12);
	gsl_matrix_complex_set(H,L_1+1,L_1-1,ga12);
	}

	for(int i=L_1;i<L_1+L_c-2;i++){
		if(i%2==1){
			gsl_matrix_complex_set(H,i,i+1,tc);
			gsl_matrix_complex_set(H,i+1,i,tc);
		}
		if(i%2==0){
			gsl_matrix_complex_set(H,i,i+3,tc);
			gsl_matrix_complex_set(H,i+3,i,tc);
		}
		gsl_matrix_complex_set(H,i,i+2,gsl_complex_mul(gsl_complex_rect(0,pow(-1,i+1)),tc));
		gsl_matrix_complex_set(H,i+2,i,gsl_complex_mul(gsl_complex_rect(0,pow(-1,i)),tc));
	}

	if(L_2!=0){
	gsl_matrix_complex_set(H,L_1+L_c-1,L_1+L_c,ga21);
	gsl_matrix_complex_set(H,L_1+L_c,L_1+L_c-1,ga21);
	gsl_matrix_complex_set(H,L_1+L_c-2,L_1+L_c,ga22);
	gsl_matrix_complex_set(H,L_1+L_c,L_1+L_c-2,ga22);
	}

	for(int i=L_1+L_c;i<L;i++){
		for(int j=L_1+L_c;j<L;j++){
			if(i==j+1||i==j-1) gsl_matrix_complex_set(H,i,j,t2);
			if(i==j) gsl_matrix_complex_set(H,i,j,la);
		}
	}

	if(pbc==1){
		gsl_matrix_complex_set(H,0,L-1,t1);
		gsl_matrix_complex_set(H,L-1,0,t1);
	}

	gsl_eigen_hermv(H,eval,evec,W);
	gsl_eigen_hermv_sort(eval,evec,GSL_EIGEN_SORT_VAL_ASC);


	for(int i=0;i<L;i++){
		fprintf(fi,"%g	%g\n",GSL_REAL(gsl_matrix_complex_get(evec,i,st)),GSL_IMAG(gsl_matrix_complex_get(evec,i,st)));
		fprintf(fil,"%g\n",gsl_vector_get(eval,i));
	}

	return 0;
}
