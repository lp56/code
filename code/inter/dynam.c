#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_eigen.h>
#include<math.h>
#include<complex.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

int main(){
	
	FILE *fi;
	fi=fopen("dynam.dat","w");
	
	int st = 0;/* print which state */
	int L_c = 10, L_s=20; /* L_c = length of crutez ladder */ int L=L_c+4*L_s;
	gsl_complex t1= gsl_complex_rect(-1,0);
	gsl_complex t2= gsl_complex_rect(-1,0); /* hoppings */
	gsl_complex tc= gsl_complex_rect(-1,0);
	int pbc = 1; /* 1 = periodic, 0 = open. NOTE, FOR PBC = 1 MUST HAVE t1=t2. */
	gsl_complex la=gsl_complex_rect(3,0);
	
	gsl_matrix_complex *H=gsl_matrix_complex_alloc(L,L);
	gsl_matrix_complex *H_0=gsl_matrix_complex_alloc(L,L);
	gsl_vector *eval=gsl_vector_alloc(L);
	gsl_matrix_complex *evec=gsl_matrix_complex_alloc(L,L);
	gsl_eigen_hermv_workspace *W = gsl_eigen_hermv_alloc(L);
	
	
	
	for(int i=0;i<2*L_s;i++){
		gsl_matrix_complex_set(H,i,i+2,t1);
		gsl_matrix_complex_set(H,i,i+2,t1);
		gsl_matrix_complex_set(H,i+2,i,t1);
		if(i!=2*L_s-1 && i!=2*L_s-2){
			gsl_matrix_complex_set(H_0,i,i,la);
			gsl_matrix_complex_set(H_0,i+2,i,t1);
			gsl_matrix_complex_set(H_0,i,i,la);
		}
	
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
	
	for(int i=L_c+2*L_s-2;i<L-2;i++){
		gsl_matrix_complex_set(H,i,i+2,t2);
		gsl_matrix_complex_set(H,i+2,i,t2);
		if(i>L_c+2*L_s-1){ 
			gsl_matrix_complex_set(H,i,i,la);
			gsl_matrix_complex_set(H_0,i,i,la);
			gsl_matrix_complex_set(H_0,i,i+2,t2);
			gsl_matrix_complex_set(H_0,i+2,i,t2);
		}
	}
	
	if(pbc==1){
		gsl_matrix_complex_set(H,0,L-2,t1);
		gsl_matrix_complex_set(H,L-2,0,t1);
		gsl_matrix_complex_set(H,L-1,1,t2);
		gsl_matrix_complex_set(H,1,L-1,t2);
		gsl_matrix_complex_set(H_0,0,L-2,t1);
		gsl_matrix_complex_set(H_0,L-2,0,t1);
		gsl_matrix_complex_set(H_0,L-1,1,t2);
		gsl_matrix_complex_set(H_0,1,L-1,t2);
	}
	
	gsl_eigen_hermv(H,eval,evec,W);
	gsl_eigen_hermv_sort(eval,evec,GSL_EIGEN_SORT_VAL_ASC);


	double re[L*L], im[L*L];

	for(int i=0;i<L;i++){for(int j=0;j<L;j++){
		re[i+L*j]=GSL_REAL(gsl_matrix_complex_get(H_0,i,j));
		im[i+L*j]=GSL_IMAG(gsl_matrix_complex_get(H_0,i,j));
		}}
	
	int func(double t, const double y[], double f[], void *params){
			for(int i=0;i<L;i++){
				for(int j=0;j<L;j++){
					f[i]=re[i+L*j]*y[j]-im[i+L*j]*y[j+L];
					f[i+L]=im[i+L*j]*y[j]+re[i+L*j]*y[j+L];
				}
			}
		return GSL_SUCCESS;
	}
	
	
	
	gsl_odeiv2_system sys = {func,NULL,2*L,NULL};
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,0.01,0.001,0.001);
	
	
	double t=0.0, tf=100.0;
	
	double y[2*L];
	for(int i=0;i<L;i++){
		y[i]=GSL_REAL(gsl_matrix_complex_get(evec,i,0));
		y[i+L]=GSL_IMAG(gsl_matrix_complex_get(evec,i,0));
	}
	
	y[3]=0.2;
	y[5]=0.7;
	
	for(int i=1;i<=50;i++){
		double ti=i*tf/100.0;
		int status = gsl_odeiv2_driver_apply(d,&t,ti,y);
	
		if(status!=GSL_SUCCESS){
			printf("error number %d, bitch!\n", status);
			break;
		}
		for(int j=0;j<L;j++)fprintf(fi,"%i	%e	%e\n",j,y[j],y[j+L]);
	
		fprintf(fi,"\n\n");
	
	}
	
	gsl_odeiv2_driver_free;
	
	return 0;
	
}
