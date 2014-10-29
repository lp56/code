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
	return 2*psif(k,3,6,7);
}

double complex g(double k[3]){
	return 2*psif(k,3,4,8);
}

int main(){

FILE *t;
char *m="w";
t=fopen("state",m);

for(int M=0;M<2;M++){

fprintf(t,"\n\n\n eigenstate %i\n\n",M+1);

double k0[3],k1[3],k2[3],k3[3];

if(M==0){
	k0[0]=0;k0[1]=1;k0[2]=5;
	k1[0]=0;k1[1]=2;k1[2]=4;
	k2[0]=1;k2[1]=2;k2[2]=3;
	k3[0]=3;k3[1]=4;k3[2]=5;
}

if(M==1){
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

double complex PSI(int x1, int x2, int x3, double complex B, double complex C, double complex D){
	return psifb(k0,x1,x2,x3)-psifb(k1,x1,x2,x3)+psifb(k2,x1,x2,x3)+psifb(k3,x1,x2,x3);
}

int x[3], xm1[3], xm2[3];


double complex HPSI;

for(int i=0; i<12; i++){
	for(int j=i; j<12; j++){
		for(int k=j; k<12; k++){
			if(i!=j && j!=k && i!=k){

				x[0]=i;x[1]=j;x[2]=k;
				xm1[0]=i-1;xm1[1]=j-1;xm1[2]=k-1;
				xm2[0]=i-2;xm2[1]=j-2;xm2[2]=k-2;
				for(int r=0; r<3;r++){
					if(xm1[r]<0) xm1[r]+=12;
					if(xm2[r]<0) xm2[r]+=12;
				 HPSI=
	(1+pow(-1,x[0]))/2*(PSI((x[0]+2)%12,x[1],x[2],B,C,D)+PSI((xm2[0])%12,x[1],x[2],B,C,D))
	+(1+pow(-1,x[1]))/2*(PSI(x[0],(x[1]+2)%12,x[2],B,C,D)+PSI(x[0],(xm2[1])%12,x[2],B,C,D))
	+(1+pow(-1,x[2]))/2*(PSI(x[0],x[1],(x[2]+2)%12,B,C,D)+PSI(x[0],x[1],(xm2[2])%12,B,C,D))
	+sqrt(2)*(PSI((x[0]+1)%12,x[1],x[2],B,C,D)+PSI((xm1[0])%12,x[1],x[2],B,C,D)
	+PSI(x[0],(x[1]+1)%12,x[2],B,C,D)+PSI(x[0],(xm1[1])%12,x[2],B,C,D)
	+PSI(x[0],x[1],(x[2]+1)%12,B,C,D)+PSI(x[0],x[1],(xm1[2])%12,B,C,D));
	}
	
			double complex energy=HPSI/PSI(x[0],x[1],x[2],B,C,D);
	

	if (pow(cimag(PSI(x[0],x[1],x[2],B,C,D)),2)>pow(10,-20) || 
	pow(creal(PSI(x[0],x[1],x[2],B,C,D)),2)>pow(10,-20)){ fprintf(t,
	"%i,%i,%i     %g\n",x[0],x[1],x[2],cimag(PSI(x[0],x[1],x[2],B,C,D)));}

			}
		}
	}
}

}


return 0;}



