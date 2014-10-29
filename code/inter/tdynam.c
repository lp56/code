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
		gsl_matrix_complex_set(H,0,L-2,t1);
		gsl_matrix_complex_set(H,L-2,0,t1);
		gsl_matrix_complex_set(H,L-1,1,t2);
		gsl_matrix_complex_set(H_0,1,L-1,t2);
		gsl_matrix_complex_set(H_0,L-2,0,t1);
		gsl_matrix_complex_set(H_0,L-1,1,t2);
		gsl_matrix_complex_set(H_0,1,L-1,t2);
	}

	gsl_eigen_hermv(H,eval,evec,W);
	gsl_eigen_hermv_sort(eval,evec,GSL_EIGEN_SORT_VAL_ASC);

	int x[]={1,-1};

	int func(double t, const double y[], double f[], void *params){
		for(int i=0;i<2;i++){
			f[i]=x[i]*y[(i+1)%2];
		}
		return GSL_SUCCESS;
	}	


	gsl_odeiv2_system sys = {func,NULL,2,NULL};
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,1e-6,1e-6,0.0);


	double t=0.0, tf=100.0;
	
	double y[2];
	
	y[0]=0;
	y[1]=1;

	for(int i=1;i<=100;i++){
		double ti=i*tf/1000.0;
		int status = gsl_odeiv2_driver_apply(d,&t,ti,y);
	
		if(status!=GSL_SUCCESS){
			printf("error number %d, bitch!\n", status);
			break;
		}
		for(int j=0;j<L;j++) fprintf(fi,"%e	%e	%e\n",t,y[0],y[1]);


	}


	gsl_odeiv2_driver_free;
	
	return 0;

}
