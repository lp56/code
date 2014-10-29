#include<math.h>
#include<gsl/gsl_vector.h>

void msort(gsl_vector *v)
{

int le;
int S= (int)v->size;
int n=1;

while(S>pow(2,n)){n++;}

int N=pow(2,n);

gsl_vector *l= gsl_vector_alloc(N);

for(int i=0; i<S;i++)
	gsl_vector_set(l,i,gsl_vector_get(v,i));

for(int i=S; i<N; i++)
	gsl_vector_set(l,i,INFINITY);

double * L = malloc(N*sizeof(double));

for(int m=0;m<n;m++)
	{
	
	le=pow(2,m);	

	if(m != 0)
		{
	for(int i=0; i<N;i++)
		gsl_vector_set(l,i,*(L+i));	
		}
	for(int k=0;k<N/2/le;k++)
		{
		int i=2*k*le, j=2*k*le;
		for(int n=2*k*le; n<(2+2*k)*le; n++)
			{
			if(gsl_vector_get(l,j)>gsl_vector_get(l,le+i)){
				*(L+n)=gsl_vector_get(l,le+i); i++;}
			else {*(L+n)=gsl_vector_get(l,j); j++;}
			if (i==(1+2*k)*le)
				{
				n++;
				for( ; n<(2+2*k)*le; n++)
					{*(L+n)=gsl_vector_get(l,j); j++;}
				break;
				}
			if (j==(1+2*k)*le)
				{
				n++;
				for( ; n<(2+2*k)*le; n++)
					{*(L+n)=gsl_vector_get(l,le+i); i++;}
				break;
				}
			}
		}
	}
for(int i=0; i<S; i++)
	gsl_vector_set(v,i,*(L+i));

free(L);
gsl_vector_free(l);
}

