#include<stdio.h>
#include<complex.h>
#include<math.h>

int main(){

	FILE *fi=fopen("state","r");

	int N=250;
	double re,im,real[N],imag[N];
	double renorm, modover;
	complex double over;
	int i=0;

	while(fscanf(fi,"%lf %lf", &re,&im)!=EOF){ real[i]=re,imag[i]=im;i++;}

	renorm=sqrt((real[100]*real[100]+imag[100]*imag[100]+real[101]*real[101]+imag[101]*imag[101]));

	over=(real[100]+I*imag[100])/renorm/sqrt(2)+I*(real[101]+I*imag[101])/renorm/sqrt(2);

	modover=over*conj(over);
	modover=sqrt(modover);

	printf("%g\n",modover);

	return 0;
}
