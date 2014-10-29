#include<stdio.h>
#include<math.h>

int main(){
int N=10;
int x=1;
int y=1;
int z=3;
	int n=0;
	for(int i=0; i<x;i++)
		n+=0.5*(pow((N-i),2)+(N-i));
	for(int i=0;i<y-x;i++)
		n+=(N-i-x);
	for(int i=0;i<z-y; i++)
		n++;
printf("%i\n",n);
return 0;
}




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

