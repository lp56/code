#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_complex.h>


int index(int *x, int *y, int *z, int N){
	int n=0;
	for(int i=0; i<*x;i++)
		n+=0.5*(pow((N-i),2)+(N-i));
	for(int i=0;i<*y-*x;i++)
		n+=(N-i-*x);
	for(int i=0;i<*z-*y; i++)
		n++;
return n;
}


void sort(int *x, int *y, int *z){
	int xyz[]={0,0,0};
	if(*x>=*y && *x>=*z){
		if(*y>*z){
			xyz[0]=*z;
			xyz[1]=*y;
			xyz[2]=*x;
		}
		else{
			xyz[0]=*y;
			xyz[1]=*z;
			xyz[2]=*x;
		}
	}
	else if(*y>=*x && *y>=*z){
		if(*z>*x){
			xyz[0]=*x;
			xyz[1]=*z;
			xyz[2]=*y;
		}
		else{
			xyz[0]=*z;
			xyz[1]=*x;
			xyz[2]=*y;
		}
	}
	else if(*z>=*y && *z>=*x){
		if(*x>*y){
			xyz[0]=*y;
			xyz[1]=*x;
			xyz[2]=*z;
		}
		else{
			xyz[0]=*x;
			xyz[1]=*y;
			xyz[2]=*z;
		}
	}
	*x=xyz[0];
	*y=xyz[1];
	*z=xyz[2];
}
void findxyz(int i, int *x, int *y, int *z, int N){
	int q=0;
	int j=(pow(N,2)+N)*0.5;
	while (i>j-1){ q++; j+=(pow(N-q,2)+(N-q))*0.5;}
	*x=q;
	int k=N-q;
	while (i>j-(pow(N-*x,2)+(N-*x))*0.5+k-1){q++; k+=(N-q);}
	*y=q;
	*z=i-j+(pow(N-*x,2)+(N-*x))*0.5-k+N-q+*y;
}


void xhop_l1(int i, gsl_matrix *M, int N, double t, double H[2]){
	int a=0,b=0,c=0,d=0;
	int *x=&a, *y=&b, *z=&c, *xn=&d;
	findxyz(i,x,y,z,N);
	if(*x==0) *xn=N-1;
	else *xn=*x-1;
	if(*x==*y && *y==*z)
		H[0]*=sqrt(2)/sqrt(6);
	else if(*x == *y || *x == *z){
		if(*xn != *y && *xn != *z)
			H[0]*=1/sqrt(2);
	}
	if(*xn==*y && *y==*z)
		H[0]*=sqrt(6)/sqrt(2);
	else if(*xn == *y || *xn == *z){
		if(*x != *y && *x != *z)
			H[0]*=sqrt(2);
	}
	sort(xn,y,z);
	gsl_matrix_set(M,index(xn,y,z,N),i,H[0]+gsl_matrix_get(M,index(xn,y,z,N),i));
}
void xhop_l2(int i, gsl_matrix *M, int N, double t, double H[2]){
	int a=0,b=0,c=0,d=0;
	int *x=&a, *y=&b, *z=&c, *xn=&d;
	findxyz(i,x,y,z,N);
	if(*x%2==0){
		if(*x==0) *xn=N-2;
		else *xn=*x-2;
		if(*x==*y && *y==*z)
			H[1]*=sqrt(2)/sqrt(6);
		else if(*x == *y || *x == *z){
			if(*xn != *y && *xn != *z)
				H[1]*=1/sqrt(2);
		}
		if(*xn==*y && *y==*z)
			H[1]*=sqrt(6)/sqrt(2);
		else if(*xn == *y || *xn == *z){
			if(*x != *y && *x != *z)
				H[1]*=sqrt(2);
		}
		sort(xn,y,z);
		gsl_matrix_set(M,index(xn,y,z,N),i,H[1]+gsl_matrix_get(M,index(xn,y,z,N),i));
	}
}
void xhop_r1(int i, gsl_matrix *M, int N, double t, double H[2]){
	int a=0,b=0,c=0,d=0;
	int *x=&a, *y=&b, *z=&c, *xn=&d;
	findxyz(i,x,y,z,N);
	if(N%2==0){
		if(*x==N-1) *xn=0;
		else *xn=*x+1;
	}
	else{
		if(*x==N-1) *xn=0;
		else *xn=*x+1;
	}
	if(*x==*y && *y==*z)
		H[0]*=sqrt(2)/sqrt(6);
	else if(*x == *y || *x == *z){
		if(*xn != *y && *xn != *z)
			H[0]*=1/sqrt(2);
		}
	if(*xn==*y && *y==*z)
		H[0]*=sqrt(6)/sqrt(2);
	else if(*xn == *y || *xn == *z){
		if(*x != *y && *x != *z)
			H[0]*=sqrt(2);
		}
	sort(xn,y,z);
	gsl_matrix_set(M,index(xn,y,z,N),i,H[0]+gsl_matrix_get(M,index(xn,y,z,N),i));
}
void xhop_r2(int i, gsl_matrix *M, int N, double t, double H[2]){
	int a=0,b=0,c=0,d=0;
	int *x=&a, *y=&b, *z=&c, *xn=&d;
	findxyz(i,x,y,z,N);
	if(*x%2==0){
		if(N%2==0){
			if(*x==N-1) *xn=1;
			else if(*x==N-2) *xn=0;
			else *xn=*x+2;
		}
		else{
			if(*x==N-1) *xn=1;
			else *xn=*x+2;
		}
		if(*x==*y && *y==*z)
			H[1]*=sqrt(2)/sqrt(6);
		else if(*x == *y || *x == *z){
			if(*xn != *y && *xn != *z)
				H[1]*=1/sqrt(2);
		}
		if(*xn==*y && *y==*z)
			H[1]*=sqrt(6)/sqrt(2);
		else if(*xn == *y || *xn == *z){
			if(*x != *y && *x != *z)
				H[1]*=sqrt(2);
		}
		sort(xn,y,z);
		gsl_matrix_set(M,index(xn,y,z,N),i,H[1]+gsl_matrix_get(M,index(xn,y,z,N),i));
	}
}
void yhop_l1(int i, gsl_matrix *M, int N, double t, double H[2]){
	int a=0,b=0,c=0,d=0;
	int *x=&a, *y=&b, *z=&c, *yn=&d;
	findxyz(i,x,y,z,N);
	if(*y==0) *yn=N-1;
	else *yn=*y-1;
	if(*x==*y && *y==*z)
		H[0]*=sqrt(2)/sqrt(6);
	else if(*x == *y || *y == *z){
		if(*yn != *x && *yn != *z)
			H[0]*=1/sqrt(2);
	}
	if(*yn==*x && *x==*z)
		H[0]*=sqrt(6)/sqrt(2);
	else if(*yn == *x || *yn == *z){
		if(*y != *x && *y != *z)
			H[0]*=sqrt(2);
	}
	sort(x,yn,z);
	gsl_matrix_set(M,index(x,yn,z,N),i,H[0]+gsl_matrix_get(M,index(x,yn,z,N),i));
}
void yhop_l2(int i, gsl_matrix *M, int N, double t, double H[2]){
	int a=0,b=0,c=0,d=0;
	int *x=&a, *y=&b, *z=&c, *yn=&d;
	findxyz(i,x,y,z,N);
	if(*y%2==0){
		if(*y==0) *yn=N-2;
		else if(*y==1) *yn=N-1;
		else *yn=*y-2;
		if(*x==*y && *y==*z)
			H[1]*=sqrt(2)/sqrt(6);
		else if(*x == *y || *y == *z){
			if(*yn != *x && *yn != *z)
				H[1]*=1/sqrt(2);
		}
		if(*yn==*x && *x==*z)
			H[1]*=sqrt(6)/sqrt(2);
		else if(*yn == *x || *yn == *z){
			if(*y != *x && *y != *z)
				H[1]*=sqrt(2);
		}
		sort(x,yn,z);
		gsl_matrix_set(M,index(x,yn,z,N),i,H[1]+gsl_matrix_get(M,index(x,yn,z,N),i));
	}
}
void yhop_r1(int i, gsl_matrix *M, int N, double t, double H[2]){
	int a=0,b=0,c=0,d=0;
	int *x=&a, *y=&b, *z=&c, *yn=&d;
	findxyz(i,x,y,z,N);
	if(N%2==0){
		if(*y==N-1) *yn=0;
		else *yn=*y+1;
	}

	else{
		if(*y==N-1) *yn=0;
		else *yn=*y+1;
	}
	if(*x==*y && *y==*z)
		H[0]*=sqrt(2)/sqrt(6);
	else if(*x == *y || *y == *z){
		if(*yn != *x && *yn != *z)
			H[0]*=1/sqrt(2);
	}
	if(*yn==*x && *x==*z)
		H[0]*=sqrt(6)/sqrt(2);
	else if(*yn == *x || *yn == *z){
		if(*y != *x && *y != *z)
			H[0]*=sqrt(2);
	}
	sort(x,yn,z);
	gsl_matrix_set(M,index(x,yn,z,N),i,H[0]+gsl_matrix_get(M,index(x,yn,z,N),i));
}
void yhop_r2(int i, gsl_matrix *M, int N, double t, double H[2]){
	int a=0,b=0,c=0,d=0;
	int *x=&a, *y=&b, *z=&c, *yn=&d;
	findxyz(i,x,y,z,N);
	if(*y%2==0){
		if(N%2==0){
			if(*y==N-1) *yn=1;
			else if(*y==N-2) *yn=0;
			else *yn=*y+2;
		}
		else{
			if(*y==N-1) *yn=1;
			else *yn=*y+2;
		}
		if(*x==*y && *y==*z)
			H[1]*=sqrt(2)/sqrt(6);
		else if(*x == *y || *y == *z){
			if(*yn != *x && *yn != *z)
				H[1]*=1/sqrt(2);
		}
		if(*yn==*x && *x==*z)
			H[1]*=sqrt(6)/sqrt(2);
		else if(*yn == *x || *yn == *z){
			if(*y != *x && *y != *z)
				H[1]*=sqrt(2);
		}
		sort(x,yn,z);
		gsl_matrix_set(M,index(x,yn,z,N),i,H[1]+gsl_matrix_get(M,index(x,yn,z,N),i));
	}
}

void zhop_l1(int i, gsl_matrix *M, int N, double t, double H[2]){
	int a=0,b=0,c=0,d=0;
	int *x=&a, *y=&b, *z=&c, *zn=&d;
	findxyz(i,x,y,z,N);
	if(*z==0) *zn=N-1;
	else *zn=*z-1;
	if(*x==*z && *y==*z)
		H[0]*=sqrt(2)/sqrt(6);
	else if(*z == *y || *x == *z){
		if(*zn != *x && *zn != *y)
			H[0]*=1/sqrt(2);
	}
	if(*zn==*x && *x==*y)
		H[0]*=sqrt(6)/sqrt(2);
	else if(*zn == *x || *zn == *y){
		if(*z != *x && *y != *z)
			H[0]*=sqrt(2);
	}
	sort(x,y,zn);
	gsl_matrix_set(M,index(x,y,zn,N),i,H[0]+gsl_matrix_get(M,index(x,y,zn,N),i));
}

void zhop_l2(int i, gsl_matrix *M, int N, double t, double H[2]){
	int a=0,b=0,c=0,d=0;
	int *x=&a, *y=&b, *z=&c, *zn=&d;
	findxyz(i,x,y,z,N);
	if(*z%2==0){
		if(*z==0) *zn=N-2;
		else if(*z==1) *zn=N-1;
		else *zn=*z-2;
		if(*x==*y && *y==*z)
			H[1]*=sqrt(2)/sqrt(6);
		else if(*x == *z || *y == *z){
			if(*zn != *x && *zn != *y)
				H[1]*=1/sqrt(2);
		}
		if(*zn==*x && *x==*y)
			H[1]*=sqrt(6)/sqrt(2);
		else if(*zn == *x || *zn == *y){
			if(*z != *x && *y != *z)
				H[1]*=sqrt(2);
		}
		sort(x,y,zn);
		gsl_matrix_set(M,index(x,y,zn,N),i,H[1]+gsl_matrix_get(M,index(x,y,zn,N),i));
	}
}

void zhop_r1(int i, gsl_matrix *M, int N, double t, double H[2]){
	int a=0,b=0,c=0,d=0;
	int *x=&a, *y=&b, *z=&c, *zn=&d;
	findxyz(i,x,y,z,N);
	if(N%2==0){
		if(*z==N-1) *zn=0;
		else *zn=*z+1;
	}

	else{
		if(*z==N-1) *zn=0;
		else *zn=*z+1;
	}
	if(*x==*y && *y==*z)
		H[0]*=sqrt(2)/sqrt(6);
	else if(*z == *x || *y == *z){
		if(*zn != *x && *zn != *y)
			H[0]*=1/sqrt(2);
	}
	if(*zn==*x && *x==*y)
		H[0]*=sqrt(6)/sqrt(2);
	else if(*zn == *x || *zn == *y){
		if(*z != *x && *y != *z)
			H[0]*=sqrt(2);
	}
	sort(x,y,zn);
	gsl_matrix_set(M,index(x,y,zn,N),i,H[0]+gsl_matrix_get(M,index(x,y,zn,N),i));
}

void zhop_r2(int i, gsl_matrix *M, int N, double t, double H[2]){
	int a=0,b=0,c=0,d=0;
	int *x=&a, *y=&b, *z=&c, *zn=&d;
	findxyz(i,x,y,z,N);
	if(*z%2==0){
		if(N%2==0){
			if(*z==N-1) *zn=1;
			else if(*z==N-2) *zn=0;
			else *zn=*z+2;
		}

		else{
			if(*z==N-1) *zn=1;
			else *zn=*z+2;
		}
		if(*x==*y && *y==*z)
			H[1]*=sqrt(2)/sqrt(6);
		else if(*x == *z || *y == *z){
			if(*zn != *x && *zn != *y)
				H[1]*=1/sqrt(2);
		}
		if(*zn==*x && *x==*y)
			H[1]*=sqrt(6)/sqrt(2);
		else if(*zn == *x || *zn == *y){
			if(*z != *x && *y != *z)
				H[1]*=sqrt(2);
		}
		sort(x,y,zn);
		gsl_matrix_set(M,index(x,y,zn,N),i,H[1]+gsl_matrix_get(M,index(x,y,zn,N),i));
	}
}

void onsite(int i, gsl_matrix *M, int N, double U){
	int a=0,b=0, c=0;
	int *x=&a, *y=&b, *z=&c;
	findxyz(i,x,y,z,N);
	if (*x==*y && *y==*z) gsl_matrix_set(M,index(x,y,z,N),i,9*U+gsl_matrix_get(M,index(x,y,z,N),i));
	else if(*x==*y||*x==*z||*y==*z) gsl_matrix_set(M,index(x,y,z,N),i,4*U+gsl_matrix_get(M,index(x,y,z,N),i));
}

int main (){


FILE *file;
char *m="w";
file=fopen("state",m);

int N=10;
int S=(int)(((double)1/12)*N*(N+1)*(2*N+1)+((double)1/4)*N*(N+1));
int Ar[(int)pow(S,2)][3];
double t=sqrt(2);
double H[]={t,1};


int e=1; /* e=0: print eigenvectors with eigenvalues. e=1: just eigenvalues. e=2: print groundstate energy as a function of U. */


if(e==2){
	for(double u=0.01;u<20;u+=0.01){



		gsl_eigen_symmv_workspace *W=gsl_eigen_symmv_alloc(S);
		gsl_vector *evals = gsl_vector_alloc(S);
		gsl_matrix *evecs = gsl_matrix_alloc(S,S);
		gsl_matrix *M=gsl_matrix_alloc(S,S);

		gsl_matrix_set_zero(M);

		double U=u;


		for(int i=0;i<S;i++){
			xhop_l1(i,M,N,t,H);
			H[0]=t;H[1]=1;
			xhop_l2(i,M,N,t,H);
			H[0]=t;H[1]=1;
			xhop_r1(i,M,N,t,H);
			H[0]=t;H[1]=1;
			xhop_r2(i,M,N,t,H);
			H[0]=t;H[1]=1;
			yhop_l1(i,M,N,t,H);
			H[0]=t;H[1]=1;
			yhop_l2(i,M,N,t,H);
			H[0]=t;H[1]=1;
			yhop_r1(i,M,N,t,H);
			H[0]=t;H[1]=1;
			yhop_r2(i,M,N,t,H);
			H[0]=t;H[1]=1;
			zhop_l1(i,M,N,t,H);
			H[0]=t;H[1]=1;
			zhop_l2(i,M,N,t,H);
			H[0]=t;H[1]=1;
			zhop_r1(i,M,N,t,H);
			H[0]=t;H[1]=1;
			zhop_r2(i,M,N,t,H);
			H[0]=t;H[1]=1;
			onsite(i,M,N,U);
		}

		gsl_eigen_symmv(M,evals,evecs,W);

		gsl_eigen_symmv_sort(evals,evecs,GSL_EIGEN_SORT_VAL_ASC);

		printf("%g %g\n",U, gsl_vector_get(evals,0));

		gsl_vector_free(evals);
		gsl_matrix_free(evecs);
		gsl_matrix_free(M);
		gsl_eigen_symmv_free(W);
	}
	
}

else{



gsl_eigen_symmv_workspace *W=gsl_eigen_symmv_alloc(S);
gsl_vector *evals = gsl_vector_alloc(S);
gsl_matrix *evecs = gsl_matrix_alloc(S,S);
gsl_matrix *M=gsl_matrix_alloc(S,S);



	double U=10000;


	for(int i=0;i<S;i++){
		xhop_l1(i,M,N,t,H);
		H[0]=t;H[1]=1;
		xhop_l2(i,M,N,t,H);
		H[0]=t;H[1]=1;
		xhop_r1(i,M,N,t,H);
		H[0]=t;H[1]=1;
		xhop_r2(i,M,N,t,H);
		H[0]=t;H[1]=1;
		yhop_l1(i,M,N,t,H);
		H[0]=t;H[1]=1;
		yhop_l2(i,M,N,t,H);
		H[0]=t;H[1]=1;
		yhop_r1(i,M,N,t,H);
		H[0]=t;H[1]=1;
		yhop_r2(i,M,N,t,H);
		H[0]=t;H[1]=1;
		zhop_l1(i,M,N,t,H);
		H[0]=t;H[1]=1;
		zhop_l2(i,M,N,t,H);
		H[0]=t;H[1]=1;
		zhop_r1(i,M,N,t,H);
		H[0]=t;H[1]=1;
		zhop_r2(i,M,N,t,H);
		H[0]=t;H[1]=1;
		onsite(i,M,N,U);
	}

	gsl_eigen_symmv(M,evals,evecs,W);

	gsl_eigen_symmv_sort(evals,evecs,GSL_EIGEN_SORT_VAL_ASC);

	if(e==0){

		for (int i=0;i<S;i++){
			for(int j=0;j<S;j++){
					int o=0,p=0,q=0;
					int *a=&o, *b=&p, *c=&q;
					findxyz(j,a,b,c,N);
					Ar[S*i+j][0]=*a;
					Ar[S*i+j][1]=*b;
					Ar[S*i+j][2]=*c; 
					}
			}

		for (int i=0; i<S; i++){
			printf("Eigenvalue %i = %g\n",i,gsl_vector_get(evals,i));
			printf("Eigenvector %i =\n",i);
			for(int j=0; j<S; j++){
					if(pow(gsl_matrix_get(evecs,j,i),2)>0.000001)
					printf("%i %i %i %g\n",Ar[S*i+j][0],Ar[S*i+j][1],Ar[S*i+j][2],gsl_matrix_get(evecs,j,i));
					if(Ar[S*i+j][2]==N) printf("\n");
			}
			printf("\n\n");
		}
	}

	if(e==1){

		for (int i=0;i<S;i++)
			fprintf(file,"%i %g\n",i,gsl_vector_get(evals,i));
	}


	gsl_vector_free(evals);
	gsl_matrix_free(evecs);
	gsl_matrix_free(M);
	gsl_eigen_symmv_free(W);
}
return 0;
}

