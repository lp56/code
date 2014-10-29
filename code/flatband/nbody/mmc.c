#include<stdio.h>

int main(){

int N=6;
int x=2;

FILE *t;
char *m="w";

t = fopen("mc.c",m);

fprintf(t,"#include<stdio.h>\nint main(){\n");
fprintf(t,"int x=%i;\nint N=%i;\nFILE *t;\nchar *m=\"w\";\nt=fopen(\"c.c\",m);\n\n",x,N);


for(int i=0; i<N-3; i++){
	fprintf(t,"for(int k0=0; k0<N; k0++){\n");
	for(int k=1; k<=i;k++)
		fprintf(t,"for(int k%i=k%i+1; k%i<N; k%i++){\n",k,k-1,k,k);
	fprintf(t,"if(k0 != %i){\n",x);
	fprintf(t,"fprintf(t,\"else if(\");\n");
	fprintf(t,"for(int l=0; l<N; l++){\n");
	fprintf(t,"if( l != %i && l != k0 ",x);
	for(int k=1; k<=i;k++)
		fprintf(t,"&& l !=k%i ",k);
	fprintf(t,") ");
	fprintf(t,"fprintf(t,\" x%%i==x%%i &&\",x,l);\n");
	fprintf(t,"}");
	fprintf(t,"fprintf(t,\" x%%i==x%%i){\\n \",x,x);\n");
	fprintf(t,"int ar[%i];\n",i+1);
	for(int k=0; k<=i; k++)
		fprintf(t,"ar[%i]=k%i;\n",k,k);
		
	for(int z=0; z<i; z++){
		if(z==0) fprintf(t,"fprintf(t,\"if(\");\n");
		else fprintf(t,"fprintf(t,\"else if(\");\n");
		fprintf(t,"for(int r0=0; r0<=%i; r0++){\n",i);
		for(int j=1; j<=z; j++)
			fprintf(t,"for(int r%i=r%i+1;r%i<%i;r%i++){\n",j,j-1,j,i,j);
		fprintf(t,"for(int p=0; p<=%i;p++){\n if(p!=r0 && p!=x",i);
		for(int q=1; q<=z; q++)
			fprintf(t,"&& p!=r%i ",q);
		fprintf(t,")\n");
		fprintf(t,"fprintf(t,\"x%%in == x%%i && \",x,ar[p]);");
		fprintf(t,"\n}");
		fprintf(t,"\nfprintf(t,\"0 == 0 ||\");");
		for(int j=0; j<=z;  j++)
			fprintf(t,"}");
		fprintf(t,"\nfprintf(t,\"0 == 1)\\n\");");
		fprintf(t,"\n");
	}

	for(int k=0; k<=i+1; k++)
		fprintf(t,"}");
	fprintf(t,"\n\n");	
}


fprintf(t,"\nreturn 0;}");
return 0;
}
