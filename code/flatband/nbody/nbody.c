#include<stdio.h>


int main(){

int N=5;

FILE *t;
char *m="w";

t = fopen("test.c",m);


/*print counting functions*/
fprintf(t,"int fcount(");
for(int i=0; i<N; i++)
	fprintf(t,"int *x%i, ",i);
fprintf(t,"int *xc, int c){\n	");
fprintf(t, "int ar[%i];\n	",N);
for(int i=0; i<N; i++)
	fprintf(t, "ar[%i]=*x%i;\n	",i,i);
fprintf(t, "int n=0;\n	");
fprintf(t,"for(int i=0; i<%i;i++){\n	",N);
fprintf(t,"	if(*xc==ar[i] && i!=c ) n++;\n	");
fprintf(t,"}\nreturn n;\n}\n\n");


fprintf(t,"int tcount(");
for(int i=0; i<N; i++)
	fprintf(t,"int *x%i, ",i);
fprintf(t," *xn){\n	");
fprintf(t, "int ar[%i];\n	",N);
for(int i=0; i<N; i++)
	fprintf(t, "ar[%i]=*x%i;\n	",i,i);

fprintf(t, "int n=0;\n	");
fprintf(t,"for(int i=0; i<%i;i++){\n	",N);
fprintf(t,"	if(*xn==ar[i]) n++;\n	");
fprintf(t,"}\nreturn n;\n}\n\n");


/*print hopping functions*/
for(int i=0;i<N;i++){
	fprintf(t,"void %ihop_l1(int i, ",i);
	for(int j=0; j<N; j++)
		fprintf(t, "double *x%i, ",j);
	fprintf(t, "int N, double H[2]){");
	fprintf(t,"\n	int a0=0,");
	for(int j=1;j<N; j++)
		fprintf(t," a%i=0,",j);
	fprintf(t," a%i=0;\n	",N+1);
	fprintf(t,"double *x0=&a0, ");
	for(int j=1; j<N; j++)
		fprintf(t,"*x%i=&a%i, ",j,j);
	fprintf(t,"*x%in=&a%i;\n	",i,N+1);
	fprintf(t,"find(");
	for(int j=0; j<N-1; j++)
		fprintf(t,"x%i,",j);
	fprintf(t,"x%i);\n	",N-1);
	fprintf(t, "if(*x%i==0) *x%in=%i-1;\n	",i,i,N);
	fprintf(t, "else *x%in=*x%i-1;\n	",i,i);
	fprintf(t, "int o = tcount(");
	for(int j=0; j<N; j++)
		fprintf(t, "x%i, ",j);
	fprintf(t, "x%in);\n	",i);
	fprintf(t,"int f = fcount(");
	for(int j=0; j<N; j++)
		fprintf(t, "x%i, ",j);
	fprintf(t,"x%i, %i);\n	", i,i);
	fprintf(t,"H[0]*=sqrt((o+1)/f);\n	");
	fprintf(t,"sort(x0");
	for(int j=1; j<N; j++)
		fprintf(t,", x%i",j);
	fprintf(t,");\n	"); 
	fprintf(t,"gsl_matrix_set(M,index(");
	for(int j=0; j<N; j++){
		if(j!=i) fprintf(t,"x%i, ",j);
		else fprintf(t, "x%in, ",j);
	}
	fprintf(t," 0)");
	fprintf(t,", t, gsl_matrix_get(M,H[0]+index(");
	for(int j=0; j<N; j++){
		if(j!=i) fprintf(t,"x%i, ",j);
		else fprintf(t, "x%in, ",j);
	}
	fprintf(t,"0)));\n}");
	

fprintf(t,"\n");
}


return 0;
}
