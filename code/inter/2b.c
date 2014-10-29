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
	fi=fopen("state.dat","w");
	fil=fopen("energy","w");
	

	int st = 0;/* print which state */
	int L_c = 20, L_s=200; /* L_c = length of crutez ladder */ int L=L_c+2*L_s;
	gsl_complex t1= gsl_complex_rect(-sqrt(2)+0.01,0);
	gsl_complex tc= gsl_complex_rect(-1,0);
	int pbc = 1; /* 1 = periodic, 0 = open. NOTE, FOR PBC = 1 MUST HAVE t1=t2. */
	gsl_complex la= gsl_complex_rect(0,0);


		gsl_matrix_complex *H=gsl_matrix_complex_alloc(L,L);
		gsl_vector *eval=gsl_vector_alloc(L);
		gsl_matrix_complex *evec=gsl_matrix_complex_alloc(L,L);
		gsl_eigen_hermv_workspace *W = gsl_eigen_hermv_alloc(L);
		
		for(int i=0;i<2*L_s;i++){
			gsl_matrix_complex_set(H,i,i+2,t1);
			gsl_matrix_complex_set(H,i+2,i,t1);
			gsl_matrix_complex_set(H,i,i,la);

		}

		for(int i=2*L_s;i<L_c+2*L_s-2;i++){
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


		gsl_eigen_hermv(H,eval,evec,W);
		gsl_eigen_hermv_sort(eval,evec,GSL_EIGEN_SORT_VAL_ASC);
		double td=GSL_REAL(t1), tcr=GSL_REAL(tc);

		for(int i=0;i<L;i++){
			fprintf(fi,"%i	%.20g	%.20g\n",i,GSL_REAL(gsl_matrix_complex_get(evec,i,st)),GSL_IMAG(gsl_matrix_complex_get(evec,i,st)));
		}
		for(int i=0;i<L;i++){
			fprintf(fil,"%.20g\n",gsl_vector_get(eval,i)-4*tcr/sqrt(4-pow(td,2)/pow(tcr,2)));
		}
	return 0;
}
