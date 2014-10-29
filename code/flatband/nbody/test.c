int fcount(int *x0, int *x1, int *x2, int *x3, int *x4, int *xc, int c){
	int ar[5];
	ar[0]=*x0;
	ar[1]=*x1;
	ar[2]=*x2;
	ar[3]=*x3;
	ar[4]=*x4;
	int n=0;
	for(int i=0; i<5;i++){
		if(*xc==ar[i] && i!=c ) n++;
	}
return n;
}

int tcount(int *x0, int *x1, int *x2, int *x3, int *x4,  *xn){
	int ar[5];
	ar[0]=*x0;
	ar[1]=*x1;
	ar[2]=*x2;
	ar[3]=*x3;
	ar[4]=*x4;
	int n=0;
	for(int i=0; i<5;i++){
		if(*xn==ar[i]) n++;
	}
return n;
}

void 0hop_l1(int i, double *x0, double *x1, double *x2, double *x3, double *x4, int N, double H[2]){
	int a0=0, a1=0, a2=0, a3=0, a4=0, a6=0;
	double *x0=&a0, *x1=&a1, *x2=&a2, *x3=&a3, *x4=&a4, *x0n=&a6;
	find(x0,x1,x2,x3,x4);
	if(*x0==0) *x0n=5-1;
	else *x0n=*x0-1;
	int o = tcount(x0, x1, x2, x3, x4, x0n);
	int f = fcount(x0, x1, x2, x3, x4, x0, 0);
	H[0]*=sqrt((o+1)/f);
	sort(x0, x1, x2, x3, x4);
	gsl_matrix_set(M,index(x0n, x1, x2, x3, x4,  0), t, gsl_matrix_get(M,H[0]+index(x0n, x1, x2, x3, x4, 0)));
}
void 1hop_l1(int i, double *x0, double *x1, double *x2, double *x3, double *x4, int N, double H[2]){
	int a0=0, a1=0, a2=0, a3=0, a4=0, a6=0;
	double *x0=&a0, *x1=&a1, *x2=&a2, *x3=&a3, *x4=&a4, *x1n=&a6;
	find(x0,x1,x2,x3,x4);
	if(*x1==0) *x1n=5-1;
	else *x1n=*x1-1;
	int o = tcount(x0, x1, x2, x3, x4, x1n);
	int f = fcount(x0, x1, x2, x3, x4, x1, 1);
	H[0]*=sqrt((o+1)/f);
	sort(x0, x1, x2, x3, x4);
	gsl_matrix_set(M,index(x0, x1n, x2, x3, x4,  0), t, gsl_matrix_get(M,H[0]+index(x0, x1n, x2, x3, x4, 0)));
}
void 2hop_l1(int i, double *x0, double *x1, double *x2, double *x3, double *x4, int N, double H[2]){
	int a0=0, a1=0, a2=0, a3=0, a4=0, a6=0;
	double *x0=&a0, *x1=&a1, *x2=&a2, *x3=&a3, *x4=&a4, *x2n=&a6;
	find(x0,x1,x2,x3,x4);
	if(*x2==0) *x2n=5-1;
	else *x2n=*x2-1;
	int o = tcount(x0, x1, x2, x3, x4, x2n);
	int f = fcount(x0, x1, x2, x3, x4, x2, 2);
	H[0]*=sqrt((o+1)/f);
	sort(x0, x1, x2, x3, x4);
	gsl_matrix_set(M,index(x0, x1, x2n, x3, x4,  0), t, gsl_matrix_get(M,H[0]+index(x0, x1, x2n, x3, x4, 0)));
}
void 3hop_l1(int i, double *x0, double *x1, double *x2, double *x3, double *x4, int N, double H[2]){
	int a0=0, a1=0, a2=0, a3=0, a4=0, a6=0;
	double *x0=&a0, *x1=&a1, *x2=&a2, *x3=&a3, *x4=&a4, *x3n=&a6;
	find(x0,x1,x2,x3,x4);
	if(*x3==0) *x3n=5-1;
	else *x3n=*x3-1;
	int o = tcount(x0, x1, x2, x3, x4, x3n);
	int f = fcount(x0, x1, x2, x3, x4, x3, 3);
	H[0]*=sqrt((o+1)/f);
	sort(x0, x1, x2, x3, x4);
	gsl_matrix_set(M,index(x0, x1, x2, x3n, x4,  0), t, gsl_matrix_get(M,H[0]+index(x0, x1, x2, x3n, x4, 0)));
}
void 4hop_l1(int i, double *x0, double *x1, double *x2, double *x3, double *x4, int N, double H[2]){
	int a0=0, a1=0, a2=0, a3=0, a4=0, a6=0;
	double *x0=&a0, *x1=&a1, *x2=&a2, *x3=&a3, *x4=&a4, *x4n=&a6;
	find(x0,x1,x2,x3,x4);
	if(*x4==0) *x4n=5-1;
	else *x4n=*x4-1;
	int o = tcount(x0, x1, x2, x3, x4, x4n);
	int f = fcount(x0, x1, x2, x3, x4, x4, 4);
	H[0]*=sqrt((o+1)/f);
	sort(x0, x1, x2, x3, x4);
	gsl_matrix_set(M,index(x0, x1, x2, x3, x4n,  0), t, gsl_matrix_get(M,H[0]+index(x0, x1, x2, x3, x4n, 0)));
}
