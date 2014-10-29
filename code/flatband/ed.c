#include<stdio.h>

int main(){

FILE *f = fopen("dmrg_results.dat","r"), *g = fopen("edd.dat","w");

int i, in[15];
double d, da[15];
int j=0;
while(fscanf(f,"%i	%lf",&i,&d)!=EOF){in[j]=i;da[j]=d*i+i-2;j++;}

for(int k=0;k<15;k++){
	fprintf(g,"%i	%g\n",in[k],da[k]);
}

return 0;
}
