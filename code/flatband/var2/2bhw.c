#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_complex.h>


int index(int *x, int *y, int N){
	int n=0;
	for(int i=0; i<*x;i++)
		n+=(N-i);
	for(int i=0;i<*y-*x; i++)
		n++;
return n;
}


void sort(int *x, int *y){
	int xy[]={0,0};
	if(*x>*y){
		xy[0]=*y;
		xy[1]=*x;
	}
	else{
		xy[0]=*x;
		xy[1]=*y;
	}
	*x=xy[0];
	*y=xy[1];
}

void findxy(int i, int *x, int *y, int N){
	int q=1;
	while (i>q*(N-0.5*(q-1))-1) q++;
	*x= q-1;
	*y=i-(q-1)*(N-((q-1)-1)*0.5)+(q-1);
}


void xhop_l1(int i, gsl_matrix *M, int N,   double H[2]){
	int a=0,b=0,d=0;
	int *x=&a, *y=&b, *xn=&d;
	findxy(i,x,y,N);
	if(*x!=0){
		*xn=*x-1;
		if(*x == *y)
			H[0]*=1/sqrt(2);
		if(*xn == *y)
			H[0]*=sqrt(2);
		sort(xn,y);
		gsl_matrix_set(M,index(xn,y,N),i,H[0]+gsl_matrix_get(M,index(xn,y,N),i));
	}
}
void xhop_l2(int i, gsl_matrix *M, int N,   double H[2], int s){
	int a=0,b=0,d=0;
	int *x=&a, *y=&b, *xn=&d;
	findxy(i,x,y,N);
	if((*x%2==0 && s==0) ||(*x%2==1 && s==1)){
		if((*x!=0 && s==0) || (*x!=1 && s==1) ){
			*xn=*x-2;
			if(*x == *y)
				H[1]*=1/sqrt(2);
			if(*xn == *y)
				H[1]*=sqrt(2);
			sort(xn,y);
			gsl_matrix_set(M,index(xn,y,N),i,H[1]+gsl_matrix_get(M,index(xn,y,N),i));
		}
	}
}
void xhop_r1(int i, gsl_matrix *M, int N,   double H[2]){
	int a=0,b=0,d=0;
	int *x=&a, *y=&b, *xn=&d;
	findxy(i,x,y,N);
	if(*x!=N-1){
		*xn=*x+1;
		if(*x == *y)
			H[0]*=1/sqrt(2);
		if(*xn == *y)
			H[0]*=sqrt(2);
		sort(xn,y);
		gsl_matrix_set(M,index(xn,y,N),i,H[0]+gsl_matrix_get(M,index(xn,y,N),i));
	}
}
void xhop_r2(int i, gsl_matrix *M, int N,   double H[2], int s){
	int a=0,b=0,d=0;
	int *x=&a, *y=&b, *xn=&d;
	findxy(i,x,y,N);
	if((*x%2==0 && s==0) ||(*x%2==1 && s==1)){
		if((*x!=N-1 && s==0) || (*x!=N-2 && s==1) ){
			*xn=*x+2;
			if(*x == *y)
				H[1]*=1/sqrt(2);
			if(*xn == *y)
				H[1]*=sqrt(2);
			sort(xn,y);
			gsl_matrix_set(M,index(xn,y,N),i,H[1]+gsl_matrix_get(M,index(xn,y,N),i));
		}
	}
}

void yhop_l1(int i, gsl_matrix *M, int N,   double H[2]){
	int a=0,b=0,d=0;
	int *x=&a, *y=&b, *yn=&d;
	findxy(i,x,y,N);
	if(*y!=0){
		*yn=*y-1;
		if(*x == *y)
			H[0]*=1/sqrt(2);
		if(*yn == *x)
			H[0]*=sqrt(2);
		sort(x,yn);
		gsl_matrix_set(M,index(x,yn,N),i,H[0]+gsl_matrix_get(M,index(x,yn,N),i));
	}
}
void yhop_l2(int i, gsl_matrix *M, int N,   double H[2], int s){
	int a=0,b=0,d=0;
	int *x=&a, *y=&b, *yn=&d;
	findxy(i,x,y,N);
	if((*y%2==0 && s==0) ||(*y%2==1 && s==1)){
		if((*y!=0 && s==0) || (*y!=1 && s==1) ){
			*yn=*y-2;
			if(*x == *y)
				H[1]*=1/sqrt(2);
			if(*x == *yn)
				H[1]*=sqrt(2);
			sort(x,yn);
			gsl_matrix_set(M,index(x,yn,N),i,H[1]+gsl_matrix_get(M,index(x,yn,N),i));
		}
	}
}
void yhop_r1(int i, gsl_matrix *M, int N,   double H[2]){
	int a=0,b=0,d=0;
	int *x=&a, *y=&b, *yn=&d;
	findxy(i,x,y,N);
	if(*y!=N-1){
		*yn=*y+1;
		if(*x == *y)
			H[0]*=1/sqrt(2);
		if(*x == *yn)
			H[0]*=sqrt(2);
		sort(x,yn);
		gsl_matrix_set(M,index(x,yn,N),i,H[0]+gsl_matrix_get(M,index(x,yn,N),i));
	}
}
void yhop_r2(int i, gsl_matrix *M, int N,   double H[2], int s){
	int a=0,b=0,d=0;
	int *x=&a, *y=&b, *yn=&d;
	findxy(i,x,y,N);
	if((*y%2==0 && s==0) ||(*y%2==1 && s==1)){
		if((*y!=N-1 && s==0) || (*y!=N-2 && s==1) ){
			*yn=*y+2;
			if(*x == *y)
				H[1]*=1/sqrt(2);
			if(*x == *yn)
				H[1]*=sqrt(2);
			sort(x,yn);
			gsl_matrix_set(M,index(x,yn,N),i,H[1]+gsl_matrix_get(M,index(x,yn,N),i));
		}
	}
}

void onsite(int i, gsl_matrix *M, int N, double U){
	int a=0,b=0;
	int *x=&a, *y=&b;
	findxy(i,x,y,N);
	if(*x==*y) gsl_matrix_set(M,index(x,y,N),i,4*U+gsl_matrix_get(M,index(x,y,N),i));
}

void set(gsl_matrix *M, double H[2], int N, int S, double t, int s,double U){
	for(int i=0;i<S;i++){
		xhop_l1(i,M,N,H);
		H[0]=t;H[1]=1;
		xhop_l2(i,M,N,H,s);
		H[0]=t;H[1]=1;
		xhop_r1(i,M,N,H);
		H[0]=t;H[1]=1;
		xhop_r2(i,M,N,H,s);
		H[0]=t;H[1]=1;
		yhop_l1(i,M,N,H);
		H[0]=t;H[1]=1;
		yhop_l2(i,M,N,H,s);
		H[0]=t;H[1]=1;
		yhop_r1(i,M,N,H);
		H[0]=t;H[1]=1;
		yhop_r2(i,M,N,H,s);
		H[0]=t;H[1]=1;
		onsite(i,M,N,U);
	}
/*	for(int i=0;i<S;i++){
		for(int j=0;j<S;j++){
			printf("%10g",gsl_matrix_get(M,i,j));
		}	
		printf("\n");
	}
*/
}
void print(int S,int N, gsl_vector *evals, gsl_matrix *evecs, int e){


	FILE *fi;
	char *m = "w";
	fi = fopen("2state",m);

	int Ar[(int)pow(S,2)][2];

	if(e==0){

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
			fprintf(fi,"Eigenvalue %i = %g\n",i,gsl_vector_get(evals,i));
			fprintf(fi,"Eigenvector %i =\n",i);
			for(int j=0; j<S; j++){
					if(pow(gsl_matrix_get(evecs,j,i),2)>0.000000000000001)
					fprintf(fi,"%i %i %g\n",Ar[S*i+j][0],Ar[S*i+j][1],gsl_matrix_get(evecs,j,i));
					if(Ar[S*i+j][1]==N) fprintf(fi,"\n");
			}
			fprintf(fi,"\n\n");
		}
	}
	if(e==1){

		for (int i=0;i<S;i++)
			fprintf(fi,"%i %g\n",i,gsl_vector_get(evals,i));
	}

}
