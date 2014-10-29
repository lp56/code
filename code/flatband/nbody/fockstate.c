#include<stdio.h>
int main(){
int R=0, T=0, f=0;int Ar[6];for(int i0=0; i0<6; i0++){
for(int i1=0; i1<6; i1++){
R+=i1;
} 
 T+=R; R=0;}Ar[0]=T;
for(int i0=0; i0<6; i0++){
for(int i1=0; i1<6; i1++){
R+=i1;
} 
 T+=R; R=0;}Ar[1]=T;
for(int i0=0; i0<6; i0++){
for(int i1=0; i1<6; i1++){
R+=i1;
} 
 T+=R; R=0;}Ar[2]=T;
for(int i0=0; i0<6; i0++){
for(int i1=0; i1<6; i1++){
R+=i1;
} 
 T+=R; R=0;}Ar[3]=T;
for(int i0=0; i0<6; i0++){
for(int i1=0; i1<6; i1++){
R+=i1;
} 
 T+=R; R=0;}Ar[4]=T;
for(int i0=0; i0<6; i0++){
for(int i1=0; i1<6; i1++){
R+=i1;
} 
 T+=R; R=0;}Ar[5]=T;
int I=Ar[0];
 while(i>I){f++;I+=Ar[f];}
