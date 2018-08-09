#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "float.h"
#include "math.h"
#include "string.h"
#include "functions.h"
#include "randomlib.h"
#include "matmul2.h"
#include "pso.h"

/* This file contains the objective functions.
   Users can define their own functions if necessary.
   If new function is defined, it also should be included
   inside functionLookup() function

*/


double Sphere(double * x,int nd){
	double sum = 0;
	int i;
	for(i = 0; i < nd; i++){
	sum += pow((x)[i],2);
	}
	return(sum);
}

double Rosenbrock(double * x, int nd){
	double sum = 0;
	int i;
	for(i=0;i<(nd-1);i++){
		sum += 100*pow(x[i+1]- pow(x[i],2),2)+pow(x[i]-1,2);
	}
	return(sum);
}

double Rastrigin(double * x, int nd)
{

	double sum = 0, pi = 3.14159265359, temp;
	int i;
	for(i=0;i<nd;i++)
	{
		temp = (double)i/5.0;
		sum += pow((x)[i] - temp,2) - (10*cos(2*pi*((x)[i] - temp)));
	}
	return(10*nd + sum);
}

double Griewank(double * x, int nd){
	double sum = 0, prod = 1;
	int i;
	for(i=0;i<nd;i++)
	{
		sum += pow((x)[i]-1,2);
		prod *= cos(((x)[i]-1)/sqrt(i+1));
	}
	return(1 + sum/4000-prod);
}

double Ackley(double * x, int nd){
   double sum1 = 0, sum2 = 0, pi = 3.14159265359, temp;
   int i;
   for ( i = 0; i < nd; i++){
      temp = (double)i/5.0;
      sum1+= pow((x)[i]- temp,2);
      sum2+= cos(2*pi*((x)[i] - temp));
   }

   return -20*exp(-0.2*sqrt(sum1/nd)) - exp(sum2/nd) + 20 + exp(1);

}

double Levy(double * x, int nd)
{
	double sum = 0,pi = 3.14159265359;
	int i;
	for(i=0;i<(nd-1);i++)
	{
		sum += pow((x)[i]-1,2)*(1 + 10*pow(sin(pi * (x)[i] + 1),2));
	}
	double a = pow(sin(pi*(x)[0]),2);
	double b = (pow((x)[nd-1]-1,2)*(1 + pow(sin(2*pi*(x)[nd-1]),2)));
	return( a + b + sum);
}


double DropWave(double * x, int nd)
{

	double sum = 0.0;
	int i;

	for(i=0;i< nd ;i++)
	{
		sum += pow((x)[i],2);
	}

	return -(1 + cos(12*sqrt(sum)))/(2 + (0.5*(sum)));
}


double Testfun(double* x, int nd){
            double sum;
            sum =  pow(x[0] - 0.2, 2)* pow(x[1] - 0.3,2)*( 1 + 20*pow(x[0] + x[1] - 1.5, 2)*pow( x[0] + 2*(x)[1] - 2.3, 2));
            return sum;
}

double Schwef(double* x, int nd){

    int i;
    double sum;
    for ( i = 0; i < nd; i++){
        sum += x[i]*sin(sqrt(abs(x[i])));
    }

   sum = 418.9829*nd - sum;
  return sum;
}



fptr functionLookup(char* fName){
if (strcmp(fName,"Sphere")==0){return &Sphere;}
else if (strcmp(fName,"Rosenbrock")==0){return &Rosenbrock;}
else if (strcmp(fName,"Rastrigin")==0){return &Rastrigin;}
else if (strcmp(fName,"Griewank")==0){return &Griewank;}
else if (strcmp(fName,"Ackley")==0){return &Ackley;}
else if (strcmp(fName,"Levy")==0) {return &Levy;}
else if (strcmp(fName,"DropWave")==0){return &DropWave;}
else if (strcmp(fName,"Schwef")==0){return &Schwef;}
else if (strcmp(fName,"Testfun")==0){return &Testfun;}
else {printf("Function Name %s not recognized\n",fName);return NULL;}
}
