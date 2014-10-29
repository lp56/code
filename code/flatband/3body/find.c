#include<stdio.h>
#include<math.h>


void findxyz(int i, int *x, int *y, int *z, int N){
	int q=0;
	int j=(pow(N,2)+N)*0.5;
	while (i>j-1){ q++; j+=(pow(N-q,2)+(N-q))*0.5;}
	*x=q;
	int k=N-q;
	while (i>j-(pow(N-*x,2)+(N-*x))*0.5+k-1){q++; k+=(N-q);}
	*y=q;
	*z=i-j+(pow(N-*x,2)+(N-*x))*0.5-k+(N-q)+*y;
}



int index(int *x, int *y, int *z, int N){
	int n=0;
	for(int i=0; i<*x;i++)
		n+=0.5*(pow((N-i),2)+(N-i));
	for(int i=0;i<*y-*x;i++)
		n+=(N-i-*x);
	for(int i=0;i<*z-*y; i++)
		n++;
return n;}

int main(){


int a=0,b=0,c=0;
int *x=&a, *y=&b, *z=&c;
int i;
int N=6;
while (i!=100){
scanf("%i",&i);

findxyz(i,x,y,z,N);

printf("%i\n",index(x,y,z,N));

}
return 0;
}



