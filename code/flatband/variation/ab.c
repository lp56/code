#include<stdio.h>	
#include<math.h>
int main(){
double e1=0.5*(1-sqrt(17));
	double c=(2+e1)/(e1*e1-2), a=sqrt(2)/e1*(c+1), A=a*sqrt(2)+c, B=sqrt(2)*(1+c), C=1+sqrt(2)*a;

printf("a=%g,c=%g\n",a,c);

return 0;}
