#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_complex.h>

void findxy(int i, int *x, int *y, int N){
	int q=1;
	while (i>q*(N-0.5*(q-1))-1) q++;
	*x= q-1;
	*y=i-(q-1)*(N-((q-1)-1)*0.5)+(q-1);
}

void xhop_l(int i, gsl_matrix *M, int N, double t){
	int a=0,b=0;
	int *x=&a, *y=&b;
	findxy(i,x,y,N);
	if(*x==0){
		if (*y==N-1) gsl_matrix_set(M,*y*(N-(*y-1)*0.5)+N-1-*y,i,t*sqrt(2)+gsl_matrix_get(M,*y*(N-(*y-1)*0.5)+N-1-*y,i));
		else if (*y==0) gsl_matrix_set(M,N-1,i,t/sqrt(2)+gsl_matrix_get(M,N-1,i));
		else gsl_matrix_set(M,*y*(N-(*y-1)*0.5)+N-1-*y,i,t+gsl_matrix_get(M,*y*(N-(*y-1)*0.5)+N-1-*y,i));
		if (*y==N-2) gsl_matrix_set(M,*y*(N-(*y-1)*0.5)+N-2-*y,i,sqrt(2)+gsl_matrix_get(M,*y*(N-(*y-1)*0.5)+N-2-*y,i));
		else if (*y==0) gsl_matrix_set(M,N-2,i,1/sqrt(2)+gsl_matrix_get(M,N-2,i));
		else gsl_matrix_set(M,*y*(N-(*y-1)*0.5)+N-2-*y,i,1+gsl_matrix_get(M,*y*(N-(*y-1)*0.5)+N-2-*y,i));
	}
	else if(*x==1){
		if(*y==1) gsl_matrix_set(M,*y,i,t/sqrt(2)+gsl_matrix_get(M,*y,i));
		else gsl_matrix_set(M,*y,i,t+gsl_matrix_get(M,*y,i));
	}

	else{
		if (*y==*x) gsl_matrix_set(M,(*x-1)*(N-(*x-2)*0.5)+*y+1-*x,i,t/sqrt(2)+gsl_matrix_get(M,(*x-1)*(N-(*x-2)*0.5)+*y+1-*x,i));
		else gsl_matrix_set(M,(*x-1)*(N-(*x-2)*0.5)+*y+1-*x,i,t+gsl_matrix_get(M,(*x-1)*(N-(*x-2)*0.5)+*y+1-*x,i));
		if(*x%2==0){
			if (*y==*x) gsl_matrix_set(M,(*x-2)*(N-(*x-3)*0.5)+*y+2-*x,i,1/sqrt(2)+gsl_matrix_get(M,(*x-2)*(N-(*x-3)*0.5)+*y+2-*x,i));
			else gsl_matrix_set(M,(*x-2)*(N-(*x-3)*0.5)+*y+2-*x,i,1+gsl_matrix_get(M,(*x-2)*(N-(*x-3)*0.5)+*y+2-*x,i));
		}
	}
}

void xhop_r(int i, gsl_matrix *M, int N, double t){
	int a=0,b=0;
	int *x=&a, *y=&b;
	findxy(i,x,y,N);
	if(*y>*x+1){
		gsl_matrix_set(M,(*x+1)*(N-*x*0.5)+*y-*x-1,i,t+gsl_matrix_get(M,(*x+1)*(N-*x*0.5)+*y-*x-1,i));
		if(*x%2==0){
			if(*y==*x+2) gsl_matrix_set(M,(*x+2)*(N-(*x+1)*0.5)+*y-*x-2,i,sqrt(2)+gsl_matrix_get(M,(*x+2)*(N-(*x+1)*0.5)+*y-*x-2,i));
			else gsl_matrix_set(M,(*x+2)*(N-(*x+1)*0.5)+*y-*x-2,i,1+gsl_matrix_get(M,(*x+2)*(N-(*x+1)*0.5)+*y-*x-2,i));
		}
	}
	if(*y==*x+1){
		if(*x==N-2){
			if (*y==*x+1) gsl_matrix_set(M,(*x+1)*(N-*x*0.5),i,t*sqrt(2)+gsl_matrix_get(M,(*x+1)*(N-*x*0.5),i));
			else gsl_matrix_set(M,(*x+1)*(N-*x*0.5),i,t+gsl_matrix_get(M,(*x+1)*(N-*x*0.5),i));
			if (*y==0) gsl_matrix_set(M,*y,i,sqrt(2)+gsl_matrix_get(M,*y,i));
			else gsl_matrix_set(M,*y,i,1+gsl_matrix_get(M,*y,i));
		}
		else{
			if(*y==*x+1) gsl_matrix_set(M,*y*(N-(*y-1)*0.5),i,t*sqrt(2)+gsl_matrix_get(M,*y*(N-(*y-1)*0.5),i));
			else gsl_matrix_set(M,*y*(N-(*y-1)*0.5),i,t+gsl_matrix_get(M,*y*(N-(*y-1)*0.5),i));
			if(*x%2==0){
				if(*y==*x+2) gsl_matrix_set(M,*y*(N-(*y-1)*0.5)+1,i,sqrt(2)+gsl_matrix_get(M,*y*(N-(*y-1)*0.5)+1,i));
				else gsl_matrix_set(M,*y*(N-(*y-1)*0.5)+1,i,1+gsl_matrix_get(M,*y*(N-(*y-1)*0.5)+1,i));
			}
		}
	}
	if(*y==*x){
		if(*x==N-1) gsl_matrix_set(M,N-1,i,t/sqrt(2)+gsl_matrix_get(M,N-1,i));
		else if(*x==N-2){
			gsl_matrix_set(M,(N-2)*(N-(N-3)*0.5)+1,i,t/sqrt(2)+gsl_matrix_get(M,(N-2)*(N-(N-3)*0.5)+1,i));
			gsl_matrix_set(M,N-2,i,1/sqrt(2)+gsl_matrix_get(M,N-2,i));
		}
		else{
			gsl_matrix_set(M,*y*(N-(*y-1)*0.5)+1,i,t/sqrt(2)+gsl_matrix_get(M,*y*(N-(*y-1)*0.5)+1,i));
			if(*y%2==0) gsl_matrix_set(M,*y*(N-(*y-1)*0.5)+2,i,1/sqrt(2)+gsl_matrix_get(M,*y*(N-(*y-1)*0.5)+2,i));
		}
	}
}
void yhop_l(int i, gsl_matrix *M, int N, double t){
	int a=0,b=0;
	int *x=&a, *y=&b;
	findxy(i,x,y,N);
	if(*y>*x+1){
		gsl_matrix_set(M,*x*(N-(*x-1)*0.5)+*y-*x-1,i,t+gsl_matrix_get(M,*x*(N-(*x-1)*0.5)+*y-*x-1,i));
		if(*y%2==0) {
			if(*y==*x+2) gsl_matrix_set(M,*x*(N-(*x-1)*0.5)+*y-*x-2,i,sqrt(2)+gsl_matrix_get(M,*x*(N-(*x-1)*0.5)+*y-*x-2,i));
			else gsl_matrix_set(M,*x*(N-(*x-1)*0.5)+*y-*x-2,i,1+gsl_matrix_get(M,*x*(N-(*x-1)*0.5)+*y-*x-2,i));
		}
	}
	if(*y==*x+1){
		gsl_matrix_set(M,*x*(N-(*x-1)*0.5)+*y-*x-1,i,t*sqrt(2)+gsl_matrix_get(M,*x*(N-(*x-1)*0.5)+*y-*x-1,i));
		if(*y%2==0) gsl_matrix_set(M,(*y-2)*(N-(*y-3)*0.5)+1,i,1+gsl_matrix_get(M,(*y-2)*(N-(*y-3)*0.5)+1,i));
	}
	if(*x==*y){
		if(*y==0){
			gsl_matrix_set(M,N-1,i,t/sqrt(2)+gsl_matrix_get(M,N-1,i));
			gsl_matrix_set(M,N-2,i,1/sqrt(2)+gsl_matrix_get(M,N-2,i));
		}
		else if(*y==1) gsl_matrix_set(M,1,i,t/sqrt(2)+gsl_matrix_get(M,1,i));
		else{
			gsl_matrix_set(M,(*y-1)*(N-(*y-2)*0.5)+1,i,t/sqrt(2)+gsl_matrix_get(M,(*y-1)*(N-(*y-2)*0.5)+1,i));
			if(*y%2==0) gsl_matrix_set(M,(*y-2)*(N-(*y-3)*0.5)+2,i,1/sqrt(2)+gsl_matrix_get(M,(*y-2)*(N-(*y-3)*0.5)+2,i));
		}
	}
}
void yhop_r(int i, gsl_matrix *M, int N, double t){
	int a=0,b=0;
	int *x=&a, *y=&b;
	findxy(i,x,y,N);
	if(*y==N-1){
		if(*x==0) gsl_matrix_set(M,*x,i,t*sqrt(2)+gsl_matrix_get(M,*x,i));
		else if(*x==*y) gsl_matrix_set(M,*x,i,t/sqrt(2)+gsl_matrix_get(M,*x,i));
		else gsl_matrix_set(M,*x,i,t+gsl_matrix_get(M,*x,i));
	}
	else if(*y==N-2){
		if (*x==*y) gsl_matrix_set(M,*x*(N-(*x-1)*0.5)+N-1-*x,i,t/sqrt(2)+gsl_matrix_get(M,*x*(N-(*x-1)*0.5)+N-1-*x,i));
		else gsl_matrix_set(M,*x*(N-(*x-1)*0.5)+N-1-*x,i,t+gsl_matrix_get(M,*x*(N-(*x-1)*0.5)+N-1-*x,i));
		if (*x==0) gsl_matrix_set(M,*x,i,sqrt(2)+gsl_matrix_get(M,*x,i));
		else gsl_matrix_set(M,*x,i,1+gsl_matrix_get(M,*x,i));
	}
	else{
		if(*x==*y) gsl_matrix_set(M,*x*(N-(*x-1)*0.5)+*y-*x+1,i,t/sqrt(2)+gsl_matrix_get(M,*x*(N-(*x-1)*0.5)+*y-*x+1,i));
		else gsl_matrix_set(M,*x*(N-(*x-1)*0.5)+*y-*x+1,i,t+gsl_matrix_get(M,*x*(N-(*x-1)*0.5)+*y-*x+1,i));
		if(*y%2==0) {
			if(*x==*y) gsl_matrix_set(M,*x*(N-(*x-1)*0.5)+*y-*x+2,i,1/sqrt(2)+gsl_matrix_get(M,*x*(N-(*x-1)*0.5)+*y-*x+2,i));
			else gsl_matrix_set(M,*x*(N-(*x-1)*0.5)+*y-*x+2,i,1+gsl_matrix_get(M,*x*(N-(*x-1)*0.5)+*y-*x+2,i));
		}
	}
}

void onsite(int i, gsl_matrix *M, int N, double U){
	int a=0,b=0;
	int *x=&a, *y=&b;
	findxy(i,x,y,N);
	if (*x==*y) gsl_matrix_set(M,*x*(N-(*x-1)*0.5),i,U+gsl_matrix_get(M,*x*(N-(*x-1)*0.5),i));
}

int main (){

/*for(int w=0;w<100;w++){*/

int N=6;
int S=(pow(N,2)+N)*0.5;
double t=sqrt(2);
double U=1000;
int Ar[(int)pow(S,2)][2];

gsl_eigen_symmv_workspace *W=gsl_eigen_symmv_alloc(S);
gsl_vector *evals = gsl_vector_alloc(S);
gsl_matrix *evecs = gsl_matrix_alloc(S,S);
gsl_matrix *M=gsl_matrix_alloc(S,S);;

for(int i=0;i<S;i++){
	xhop_l(i,M,N,t);
	xhop_r(i,M,N,t);
	yhop_l(i,M,N,t);
	yhop_r(i,M,N,t);
	onsite(i,M,N,U);
}


gsl_eigen_symmv(M,evals,evecs,W);

gsl_eigen_symmv_sort(evals,evecs,GSL_EIGEN_SORT_VAL_ASC);

for (int i=0;i<S;i++){
	for(int j=0;j<S;j++){
			int o=0,p=0;
			int *a=&o, *b=&p;
			findxy(j,a,b,N);
			Ar[S*i+j][0]=*a;
			Ar[S*i+j][1]=*b;
		}
	}


for (int i=0; i<S; i++){
	printf("Eigenvalue %i = %g\n",i,gsl_vector_get(evals,i));
	printf("Eigenvector %i =\n ",i);
	for(int j=0; j<S; j++){
			if(gsl_matrix_get(evecs,j,i)*gsl_matrix_get(evecs,j,i)>0.00000001)
			printf("%i %i %g\n",Ar[S*i+j][0],Ar[S*i+j][1],gsl_matrix_get(evecs,j,i));
			if(Ar[S*i+j][1]==N) printf("\n");
	}
	printf("\n\n");
}
for (int i=0;i<S;i++)
	printf("%i %g\n",i,gsl_vector_get(evals,i));



gsl_vector_free(evals);
gsl_matrix_free(evecs);
gsl_matrix_free(M);
gsl_eigen_symmv_free(W);

return 0;
}

