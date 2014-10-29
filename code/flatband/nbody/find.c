#include<stdio.h>



int main(){

int N=4;
int n=0;
int L=6;
int i=10;


FILE *t;
char *m="w";

t = fopen("fockstate.c",m);


fprintf(t,"#include<stdio.h>\nint main(){\n");

fprintf(t,"int R=0, T=0, f=0;");
fprintf(t,"int Ar[%i];",L);
for(int j=0; j<L; j++){
	for(int i=0;i<N-2;i++)
		fprintf(t,"for(int i%i=%i; i%i<%i; i%i++){\n",i,n,i,L,i);
	fprintf(t,"R+=i%i;\n",N-3);
	fprintf(t,"} \n T+=R; R=0;");
	for(int i=0; i<N-3; i++)
		fprintf(t,"}");
	fprintf(t,"Ar[%i]=T;\n",j);
}

fprintf(t,"int I=Ar[0];\n while(i>I){f++;I+=Ar[f];}\n");



return 0;
}
