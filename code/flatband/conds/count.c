#include<stdio.h>


int main(){

int N=3;
int Ar[N];
int Max[N];
int count[2*N];
int x;
int y;

for(int i=0;i<2*N;i++)
	count[i]=0;

for(int i=0;i<N;i++){
	Ar[i]=i;
	Max[N-i-1]=2*N-1-i;
}

while(Ar[0]<Max[0]){
	Ar[N-1]++;
	for(int j=0;j<N;j++){
		if(Ar[N-j-1]>Max[N-j-1]){
			Ar[N-j-2]+=1;
			for(int k=0;k<j+1;k++)
				Ar[N-j-1+k]=Ar[N-j-2+k]+1;
		}
	}
	x=0;
	for(int i=0;i<N;i++)
		x+=Ar[i];
	y=x%(2*N);

	count[y]+=1;

if(y==3){
	printf("3. ");
	for(int i=0;i<N;i++)
	printf("%i ",Ar[i]);
	printf("\n");
}
else if(y==0){
	printf("0. ");
	for(int i=0;i<N;i++)
	printf("%i ",Ar[i]);
	printf("\n");

}

}

count[(int)(N*(N-1)/2)%(2*N)]+=1;

printf("\n");

for(int i=0; i<2*N; i++)
	printf("%3i  ",i);
printf("\n");
for(int i=0; i<2*N; i++)
	printf("%3i  ",count[i]);

printf("\n\n");

return 0;
}


