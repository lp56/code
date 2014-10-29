#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_complex.h>
int index(int *x, int *y, int N);
void xhop_l1(int i, gsl_matrix *M, int N, double H[2]);
void xhop_l2(int i, gsl_matrix *M, int N, double H[2], int s);
void xhop_r1(int i, gsl_matrix *M, int N, double H[2]);
void xhop_r2(int i, gsl_matrix *M, int N, double H[2], int s);
void yhop_l1(int i, gsl_matrix *M, int N, double H[2]);
void yhop_l2(int i, gsl_matrix *M, int N, double H[2], int s);
void yhop_r1(int i, gsl_matrix *M, int N, double H[2]);
void yhop_r2(int i, gsl_matrix *M, int N, double H[2], int s);
void onsite(int i, gsl_matrix *M, int N, double U);
void set(gsl_matrix *M, double H[2], int N, int S,double t, int s, double U);
void alpha(double U, double H[2], double t, double alp[21], double *ec);

int delt(int i, int j){
	if(i == j) return 1;
	else return 0;
}

double vd(int i, int j){
	return sqrt(2)*delt(i,j)-delt(i,j+1)-delt(i,j-1);
}

double vdl(int i){
	return sqrt(2)*delt(i,0)-delt(i,1);
}

double vdr(int i){
	return sqrt(2)*delt(i,6)-delt(i,5);
}

int mod(int i){
	return sqrt(i*i);
}

double jev(double alp[21],int F, int h){

	int N=(F+1)/2;

	double s=0;
	int a=0, b=0, c=0, d=0;
	int *w=&a, *x=&b, *y=&c, *z=&d;
	int n;

	if(h==3||h==F-3){

		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int l=0; l<7;l++){
					for(int r=l+1; r<7; r++){
						*w=k; *x=m; *y=l; *z=r;
						s+=alp[index(w,x,7)]*alp[index(y,z,7)]*(vd(l,2)*vd(m,4)*vdr(r)*vdl(k)+vd(k,4)*vd(r,2)*delt(l+6,m));
					}
				}
			}
		}
	}

	else if(h==2||h==F-2){

		n=N-3;		
		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int l=0; l<7;l++){
					for(int r=l+1; r<7; r++){
						*w=k; *x=m; *y=l; *z=r;
						s+=alp[index(w,x,7)]*alp[index(y,z,7)]*vd(k,2)*vd(r,4)*delt(m,l+4);
					}
				}
			}
		}

		s*=pow(4,n);
	}

	else if(h==1||h==F-1){

		
		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int l=0; l<7;l++){
					for(int r=l+1; r<7; r++){
						*w=k; *x=m; *y=l; *z=r;
						s+=alp[index(w,x,7)]*alp[index(y,z,7)]*(vdr(r)*vdl(k)*delt(m,l+2)+delt(k,l+2)*delt(m,r+2));
					}
				}
			}
		}

	}

	else if(h==0){

		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int l=0; l<7;l++){
					for(int r=l+1; r<7; r++){
						*w=k; *x=m; *y=l; *z=r;
						s+=alp[index(w,x,7)]*alp[index(y,z,7)]*(delt(k,l)*delt(m,r));
					}
				}
			}
		}
		s*=pow(4,(N-2));
	}

	else{ 
		
		if(mod(h)%2==0) n=(2*F-2*mod(h)-6)/4;
		else n=(2*mod(h)-6)/4;
		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int l=0; l<7;l++){
					for(int r=l+1; r<7; r++){
						*w=k; *x=m; *y=l; *z=r;
						s+=alp[index(w,x,7)]*alp[index(y,z,7)]*vdl(k)*vd(m,4)*vdr(r)*vd(l,2);
					}
				}
			}
		}

		s*=pow(4,n);
	}
	return s;

}
double jevr(double alp[21], int F, int h){


	int N=(F+1)/2;

	double s=0;
	int  b=0, c=0, d=0, e=6;
	int  *x=&b, *y=&c, *z=&d, *f=&e;
	int n;



	if(F-4>h && h>1){

		n=h/2-3;
		if(h%2==0) return 0;
		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int p=0; p<6;p++){
						*x=k; *y=m; *z=p;
						s+=alp[index(x,y,7)]*alp[index(z,f,7)]*vdl(k)*vd(m,4)*vd(p,2)*(-sqrt(2));
				}
			}
		}

		s*=pow(4,n);
	}

	else if(h==1){

		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int p=0; p<6;p++){
						*x=k; *y=m; *z=p;
						s+=alp[index(x,y,7)]*alp[index(z,f,7)]*vdl(k)*delt(m,p+2)*(-sqrt(2));
				}
			}
		}
	}

	else if(h==F-1){

		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int p=0; p<6;p++){
						*x=k; *y=m; *z=p;
						s+=alp[index(x,y,7)]*alp[index(z,f,7)]*(delt(k+2,p)*(sqrt(2)*delt(m,5)+delt(m,6))+vdl(p)*delt(k,5)*delt(m,6));
				}
			}
		}
	}


	else if(h==F-2){

		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int p=0; p<6;p++){
						*x=k; *y=m; *z=p;
						s+=alp[index(x,y,7)]*alp[index(z,f,7)]*vd(p,2)*(delt(k,3)*delt(m,4)-sqrt(2)*delt(k,3)*delt(m,5)-delt(k,4)*delt(m,5));
				}
			}
		}

		s*=pow(4,N-3);
	}

	else if(h==F-3){

		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int p=0; p<6;p++){
						*x=k; *y=m; *z=p;
						s+=alp[index(x,y,7)]*alp[index(z,f,7)]*vd(p,2)*(delt(k,1)*delt(m,2)-sqrt(2)*delt(k,1)*delt(m,3)-delt(k,2)*delt(m,3));
				}
			}
		}
	}
	
	else if(h==F-4){


		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int p=0; p<6;p++){
						*x=k; *y=m; *z=p;
						s+=alp[index(x,y,7)]*alp[index(z,f,7)]*vd(p,2)*vd(m,4)*(delt(k,0)-sqrt(2)*vdl(k));
				}
			}
		}

		s*=pow(4,N-4);
	}

s *=sqrt(2);

return s;

}


double jevl(double alp[21], int F, int h){


	int N=(F+1)/2;

	double s=0;
	int  b=0, c=0, d=0, e=0;
	int  *x=&b, *y=&c, *z=&d, *f=&e;
	int n;



	if(F-3>h && h>4){

		n=(2*F-2*h-6)/4;
		if(h%2==1) return 0;
		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int p=1; p<7;p++){
						*x=k; *y=m; *z=p;
						s+=alp[index(x,y,7)]*alp[index(f,z,7)]*vdr(m)*vd(k,2)*vd(p,4)*(-sqrt(2));
				}
			}
		}

		s*=pow(4,n);
	}

	else if(h==1){

		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int p=1; p<7;p++){
						*x=k; *y=m; *z=p;
						s+=alp[index(x,y,7)]*alp[index(f,z,7)]*(vdr(p)*delt(m,1)*delt(k,0)+delt(m,p+2)*(sqrt(2)*delt(k,1)+delt(k,0)));
				}
			}
		}
	}

	else if(h==2){

		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int p=1; p<7;p++){
						*x=k; *y=m; *z=p;
						s+=alp[index(x,y,7)]*alp[index(f,z,7)]*vd(p,4)*(sqrt(2)*delt(m,3)*(sqrt(2)*delt(k,2)-delt(k,1))-delt(k,1)*delt(m,2)-delt(k,2)*delt(m,3));
				}
			}
		}
		s*=pow(4,N-3);
	}


	else if(h==3){

		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int p=1; p<7;p++){
						*x=k; *y=m; *z=p;
						s+=alp[index(x,y,7)]*alp[index(f,z,7)]*vd(p,2)*(sqrt(2)*delt(m,5)*(sqrt(2)*delt(k,4)-delt(k,3))-delt(k,3)*delt(m,4)-delt(k,4)*delt(m,5));
				}
			}
		}
	}

	else if(h==4){

		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int p=1; p<7;p++){
						*x=k; *y=m; *z=p;
						s+=alp[index(x,y,7)]*alp[index(f,z,7)]*vd(k,2)*vd(p,4)*(-delt(m,6)+sqrt(2)*delt(m,5));
				}
			}
		}
		s*=pow(4,N-4);
	}
	
	else if(h==F-1){


		for(int k=0;k<7;k++){
			for(int m=k+1;m<7;m++){
				for(int p=1; p<7;p++){
						*x=k; *y=m; *z=p;
						s+=alp[index(x,y,7)]*alp[index(f,z,7)]*vdr(m)*delt(p,k+2)*(-sqrt(2));
				}
			}
		}

	}

s *=sqrt(2);

return s;

}
