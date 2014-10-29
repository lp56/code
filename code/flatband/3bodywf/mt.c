#include<stdio.h>
#include<math.h>
int main(){

int x(int y){ return (1+pow(-1,y))/2;}
for(int i=0;i<12;i++){
printf("%i\n",x(i));
}
return 0;}
