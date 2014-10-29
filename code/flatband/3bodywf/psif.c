#include<stdio.h>
#include<math.h>
#include<complex.h>

#define PI 3.1415926535897832846

double complex Iexp(double x){
	return cos(x)+I*sin(x);
}

double complex p(double x, double y, double z, int i, int j, int k){
	if(i%2==1 && j%2==1 && k%2==1)
		return -2*sqrt(2)*cos(x)*cos(y)*cos(z)*Iexp(i*y+j*x+k*z);
	else if(j%2==1 && i%2==1)
		return 2*cos(x)*cos(y)*Iexp(i*y+j*x+k*z);
	else if(j%2==1 && k%2==1)
		return 2*cos(x)*cos(z)*Iexp(i*y+j*x+k*z);
	else if(i%2==1 && k%2==1)
		return 2*cos(y)*cos(z)*Iexp(i*y+j*x+k*z);
	else if(j%1==1)
		return -sqrt(2)*cos(x)*Iexp(i*y+j*x+k*z);
	else if(i%1==1)
		return -sqrt(2)*cos(y)*Iexp(i*y+j*x+k*z);
	else if(k%1==1)
		return -sqrt(2)*cos(z)*Iexp(i*y+j*x+k*z);
	else return Iexp(i*y+j*x+k*z);
}

double psif(double k[3],int x, int y, int z){
	int sp;
	if(x>y && y>z || y>z && z>x || z>x && x>y) sp=1;
	else sp=-1;
	return sp*(p(k[0],k[2],k[1],x,y,z)-p(k[0],k[1],k[2],x,y,z)+p(k[1],k[0],k[2],x,y,z)-
	p(k[1],k[2],k[0],x,y,z) +p(k[2],k[1],k[0],x,y,z)-p(k[2],k[0],k[1],x,y,z));
}


int main(){


double k[]={2*PI/6,4*PI/6,PI};

int x[]={1,2,3};


printf("%g\n",psif(k,x[0],x[1],x[2]));

return 0;}
