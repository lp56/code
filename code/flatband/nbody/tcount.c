#include<stdio.h>



int main(){

int x0=5, x1=5, x2=5, x3=6, x4=5, x5=6, x0n=6;

int ar[5];

ar[0]=x0;
ar[1]=x1;
ar[2]=x2;
ar[3]=x3;
ar[4]=x4;
ar[5]=x5;

int n=0;

for(int i=0; i<6; i++){
		if(x0n==ar[i]) n++;
}


printf("%i\n",n);

return 0;
}

