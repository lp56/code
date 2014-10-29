#include<stdio.h>
#include<math.h>
#include<complex.h>

#define PI 3.1415926535897832846

double complex Iexp(double x){
	return cos(x)+I*sin(x);
}

double complex p(double k[3], int x, int y, int z){
	if(x%2==1 && y%2==1 && z%2==1)
		return -2*sqrt(2)*cos(k[0])*cos(k[1])*cos(k[2])*Iexp(k[0]*x+k[1]*y+k[2]*z);
	else if(x%2==1 && y%2 ==1)
		return 2*cos(k[0])*cos(k[1])*Iexp(k[0]*x+k[1]*y+k[2]*z);
	else if(x%2==1 && z%2==1)
		return 2*cos(k[0])*cos(k[2])*Iexp(k[0]*x+k[1]*y+k[2]*z);
	else if(y%2==1 && z%2==1)
		return 2*cos(k[1])*cos(k[2])*Iexp(k[0]*x+k[1]*y+k[2]*z);
	else if(x%2==1)
		return sqrt(2)*cos(k[0])*Iexp(k[0]*x+k[1]*y+k[2]*z);
	else if(y%2==1)
		return sqrt(2)*cos(k[1])*Iexp(k[0]*x+k[1]*y+k[2]*z);
	else if(z%2==1)
		return sqrt(2)*cos(k[2])*Iexp(k[0]*x+k[1]*y+k[2]*z);
	else return Iexp(k[0]*x+k[1]*y+k[2]*z);
}

double complex psif(double k[3],int x, int y, int z){
	int sp;
	if((x>y && y>z) || (y>z && z>x) || (z>x && x>y)) sp=1;
	else sp=-1;
	return sp*(p(k,x,y,z)-p(k,x,z,y)+p(k,z,x,y)-p(k,z,y,x)+p(k,y,z,x)-p(k,y,x,z));
}



for 
