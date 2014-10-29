#include<stdio.h>
#include<math.h>


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
			printf("b11\n");
		}
		else{
			xyz[0]=*y;
			xyz[1]=*z;
			xyz[2]=*x;
			printf("b12\n");
		}
	}
	else if(*y>=*x && *y>=*z){
		if(*z>*x){
			xyz[0]=*x;
			xyz[1]=*z;
			xyz[2]=*y;
			printf("b21\n");
		}
		else{
			xyz[0]=*z;
			xyz[1]=*x;
			xyz[2]=*y;
			printf("b22\n");
		}
	}
	else if(*z>=*y && *z>=*x){
		if(*x>*y){
			xyz[0]=*y;
			xyz[1]=*x;
			xyz[2]=*z;
			printf("b31\n");
		}
		else{
			xyz[0]=*x;
			xyz[1]=*y;
			xyz[2]=*z;
			printf("b32\n");
		}
	}
	printf("xyz[0]=%i xyz[1]=%i xyz[2]=%i\n",xyz[0],xyz[1],xyz[2]);
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

int main (){

int N=4;
int i=19;



int a=0, b=0, c=0, d=0;
int *x=&a, *y=&b, *z=&c, *xn=&d;


findxyz(i,x,y,z,N);


printf("x=%i y=%i z=%i\n",*x,*y,*z);

*xn=1;

sort(xn,y,z);

printf("xn=%i y=%i z=%i", *xn,*y,*z);

printf("index = %i\n",index(xn,y,z,N));

return 0;
}

 

