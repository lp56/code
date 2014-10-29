#include<stdio.h>
int main(){
int x=2;
int N=6;
FILE *t;
char *m="w";
t=fopen("c.c",m);

for(int k0=0; k0<N; k0++){
if(k0 != 2){
fprintf(t,"else if(");
for(int l=0; l<N; l++){
if( l != 2 && l != k0 ) fprintf(t," x%i==x%i &&",x,l);
}fprintf(t," x%i==x%i){\n ",x,x);
int ar[1];
ar[0]=k0;
}}

for(int k0=0; k0<N; k0++){
for(int k1=k0+1; k1<N; k1++){
if(k0 != 2){
fprintf(t,"else if(");
for(int l=0; l<N; l++){
if( l != 2 && l != k0 && l !=k1 ) fprintf(t," x%i==x%i &&",x,l);
}fprintf(t," x%i==x%i){\n ",x,x);
int ar[2];
ar[0]=k0;
ar[1]=k1;
fprintf(t,"if(");
for(int r0=0; r0<=1; r0++){
for(int p=0; p<=1;p++){
 if(p!=r0 && p!=x)
fprintf(t,"x%in == x%i && ",x,ar[p]);
}
fprintf(t,"0 == 0 ||");}
fprintf(t,"0 == 1)\n");
}}}

for(int k0=0; k0<N; k0++){
for(int k1=k0+1; k1<N; k1++){
for(int k2=k1+1; k2<N; k2++){
if(k0 != 2){
fprintf(t,"else if(");
for(int l=0; l<N; l++){
if( l != 2 && l != k0 && l !=k1 && l !=k2 ) fprintf(t," x%i==x%i &&",x,l);
}fprintf(t," x%i==x%i){\n ",x,x);
int ar[3];
ar[0]=k0;
ar[1]=k1;
ar[2]=k2;
fprintf(t,"if(");
for(int r0=0; r0<=2; r0++){
for(int p=0; p<=2;p++){
 if(p!=r0 && p!=x)
fprintf(t,"x%in == x%i && ",x,ar[p]);
}
fprintf(t,"0 == 0 ||");}
fprintf(t,"0 == 1)\n");
fprintf(t,"else if(");
for(int r0=0; r0<=2; r0++){
for(int r1=r0+1;r1<2;r1++){
for(int p=0; p<=2;p++){
 if(p!=r0 && p!=x&& p!=r1 )
fprintf(t,"x%in == x%i && ",x,ar[p]);
}
fprintf(t,"0 == 0 ||");}}
fprintf(t,"0 == 1)\n");
}}}}


return 0;}
