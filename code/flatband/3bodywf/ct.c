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
		return -sqrt(2)*cos(k[0])*Iexp(k[0]*x+k[1]*y+k[2]*z);
	else if(y%2==1)
		return -sqrt(2)*cos(k[1])*Iexp(k[0]*x+k[1]*y+k[2]*z);
	else if(z%2==1)
		return -sqrt(2)*cos(k[2])*Iexp(k[0]*x+k[1]*y+k[2]*z);
	else return Iexp(k[0]*x+k[1]*y+k[2]*z);
}



double complex psif(double k[3],int x, int y, int z){
	return p(k,x,y,z)-p(k,x,z,y)+p(k,z,x,y)-p(k,z,y,x)+p(k,y,z,x)-p(k,y,x,z);
}

double complex psifb(double k[3],int x, int y, int z){
	int sp;
	if((x>y && y>z) || (y>z && z>x) || (z>x && x>y)) sp=1;
	else sp=-1;
	return sp*psif(k,x,y,z);
}

double complex e(double k[3]){
	return 2*(psif(k,1,0,3)+psif(k,1,4,3));
}

double complex d(double k[3]){
	return 2*psif(k,3,4,7);
}

double complex g(double k[3]){
	return 2*psif(k,2,3,8);
}

int main(){

int M=3;

double k0[3],k1[3],k2[3],k3[3];

if(M==0){
	k0[0]=0;k0[1]=1;k0[2]=5;
	k1[0]=0;k1[1]=2;k1[2]=4;
	k2[0]=1;k2[1]=2;k2[2]=3;
	k3[0]=3;k3[1]=4;k3[2]=5;
}

if(M==3){
	k0[0]=0;k0[1]=4;k0[2]=5;
	k1[0]=1;k1[1]=3;k1[2]=5;
	k2[0]=2;k2[1]=3;k2[2]=4;
	k3[0]=0;k3[1]=1;k3[2]=2;
}


for(int i=0;i<3;i++){
	k0[i]*=2*PI/12;
	k1[i]*=2*PI/12;
	k2[i]*=2*PI/12;
	k3[i]*=2*PI/12;
}


double complex e0=e(k0), e1=e(k1), e2=e(k2), e3=e(k3), d0=d(k0), d1=d(k1), d2=d(k2), d3=d(k3),
g0=g(k0), g1=g(k1), g2=g(k2), g3=g(k3);

double complex l1=(g3*e0/e3-g0)/(e2*g3/e3-g2), l2=(g3*e1/e3-g1)/(g2-e2*g3/e3);

double complex B=(-d0+(d2-e2*d3/e3)*l1+e0*d3/e3)/(d1+(d2-e2*d3/e3)*l2-d3*e1/e3);

double complex C=(-g0-B*g1+(e0+B*e1)*g3/e3)/(g2-e2*g3/e3);

double complex D=(-g0-B*g1-C*g2)/g3;


double complex re= e0+e1*B+e2*C+e3*D;
double complex rd= d0+d1*B+d2*C+d3*D;
double complex rg= g0+g1*B+g2*C+g3*D;

printf("%g+i%g, %g+i%g, %g+i%g\n",creal(re),cimag(re),creal(rd),cimag(rd),creal(rg),cimag(rg));

return 0;}

