#include <R.h>
#include <Rmath.h>
#include <stdlib.h>

/*******************************************************************************/

void local_alpha_update(double* warmalpha_ij, double* warmalpha_yi, double* my, double* yyi, double* xinner, double* lam, double* wweight, int* nnobs, double* nnobsdouble, int* kk, double* kkdouble, double* erci, double* ggamma, int* ytrain, double* outalpha_ij)
{

int i, j, ii, jj, q, iter, iii, jjj; 

int k = *kk;

double kdouble = *kkdouble;

int kminus = k-1;

int nobs = *nnobs;

double nobsdouble = *nnobsdouble;

double lambda = *lam;

double gamma = *ggamma;

double alpha_ij[nobs*k];

double fake[nobs*k];

double alpha_yi[nobs];

double change;

double temp, temp2, yici;

for (i=0;i<(nobs*k);i++)
{
alpha_ij[i] = warmalpha_ij[i];
}

for (i=0;i<nobs;i++)
{
alpha_yi[i] = warmalpha_yi[i];
}

/*--initiate fake--------------------------------------*/

for (i=0;i<nobs;i++)
	{
	    for (j=0;j<k;j++)
	    {
yici=0;

for (q=0;q<kminus;q++)
{
	for (ii=0;ii<nobs;ii++)
	{
yici += 2*alpha_yi[ii]*yyi[(ii+nobs*q)]*my[(j+k*q)]*xinner[(ii+i*nobs)];
	} /* ii */
	
	for (jj=0;jj<k;jj++)
	{
		for (ii=0;ii<nobs;ii++)
		{
yici -= alpha_ij[(ii+nobs*jj)]*my[(jj+k*q)]*my[(j+k*q)]*xinner[(ii+i*nobs)];
		} /* ii */
	} /* jj */

} /* q */

fake[(i+nobs*j)] = yici;				
	    }
	}



/*--update----------------------------------------------------*/


    for (iter=0;iter<10;iter++)
    {
    	for (i=0;i<nobs;i++)
	{
	    for (j=0;j<k;j++)
	    {
		if (ytrain[i]==(j+1))
		{	

yici = -(fake[(i+nobs*j)]/nobsdouble)/lambda - 2*erci[i]*alpha_ij[(i+nobs*j)] + kdouble-1;

temp = -(yici/2)/erci[i];

if (temp<0) {temp2=0;}
 	else {
		if (temp>(wweight[i]*gamma)) {temp2=(wweight[i]*gamma);}
			else {temp2=temp;}
		} 

change = temp2-alpha_ij[(i+nobs*j)];
		
alpha_yi[i]=temp2;
alpha_ij[(i+nobs*j)]=temp2;

/*--update fake[iii,jjj]-------------------------*/

for (iii=0;iii<nobs;iii++)
	{
	    for (jjj=0;jjj<k;jjj++)
	    {
for (q=0;q<kminus;q++)
{
fake[(iii+nobs*jjj)] += change*my[(j+k*q)]*my[(jjj+k*q)]*xinner[(i+iii*nobs)];
}
	    }
	}

/*--update fake-------------------------*/
		
		} /*if ==*/	
/*--------------if equal or not-----------------*/
		if (ytrain[i]!=(j+1))
		{	

yici = fake[(i+nobs*j)]/nobsdouble/lambda - 2*erci[i]*alpha_ij[(i+nobs*j)] + 1;

temp = -yici/2/erci[i];

if (temp<0) {temp2=0;}
 	else {
		if (temp>(wweight[i]*(1-gamma))) {temp2=(wweight[i]*(1-gamma));}
			else {temp2=temp;}
		} 

change = temp2-alpha_ij[(i+nobs*j)];
		
alpha_ij[(i+nobs*j)]=temp2;

/*--update fake[iii,jjj]-------------------------*/

for (iii=0;iii<nobs;iii++)
	{
	    for (jjj=0;jjj<k;jjj++)
	    {
for (q=0;q<kminus;q++)
{
fake[(iii+nobs*jjj)] -= change*my[(j+k*q)]*my[(jjj+k*q)]*xinner[(i+iii*nobs)];
}
	    }
	}

/*--update fake-------------------------*/

		} /*if !=*/

	    } /*for (j=0;j<k;j++)*/

	} /*for (i=0;i<nobs;i++)*/
    
/*check difference and compare to epsilon*/

    } /*for (iter=0;iter<2;iter++)*/

for (i=0;i<(nobs*k);i++)
{
outalpha_ij[i] = alpha_ij[i];
}

} /*void alpha_update*/
