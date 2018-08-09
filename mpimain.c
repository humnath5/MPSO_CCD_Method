#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include <math.h>
#include <time.h>
#include "randomlib.h"
#include "matmul2.h"
#include "pso.h"
#include "assert.h"
#include "functions.h"

int main(int argc,char*argv[]){

	int i,j,k, seed1,seed2;
	srand((unsigned)time(NULL));
	seed1= rand()%30000; seed2= rand()%30000;

	char funName[80];
	strcpy(funName, "Rosenbrock");
	fptr f  = functionLookup(funName);

	int np = 500, nd = 5, ni = 500, maxRun = 50;
	double l  = -10, u = 10;
	/*
	np:  number of  PSO particles
	nd:  dimension of the problem to be optimized
	ni:  maximum iterations allowed per run
	maxRun: maximum runs (re-starts) allowed
	l: lower bound of each component of the domain
	u: upper bound of each component of the domain
	*/

	double*lb = (double*)malloc(sizeof(double)*nd);
	double*ub = (double*)malloc(sizeof(double)*nd);
	double *gbest =  (double*)malloc(sizeof(double)*nd);
	for(i = 0; i < nd; i++){
		lb[i] = l; ub[i] = u;
		gbest[i] = (ub[i]-lb[i])/3;
	 }

	/*==== Setting up MPI environment ====*/
	int myrank, comsize;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Barrier(MPI_COMM_WORLD);

	RandomInitialise(myrank + seed1, seed2);

	double start,finish, ellapse_time;
	double value = 9999999;

	/*=== Start optimization ===*/
	start = MPI_Wtime(NULL);

	pmpso_ccd(f, nd, np, ni, maxRun, lb, ub, &value, gbest,myrank);

	finish = MPI_Wtime(NULL);

	ellapse_time = difftime(finish, start);

	/* === Writing outputs ====*/
	if(myrank == 0){
	printf("\n ==== PSO Parameters Used ====\n");
	printf("Best value obtained: %lf\n",value);
	printf("Random seeds: %d, %d\n",seed1,seed2);
	printf("Number of processors used: %d\n",comsize);
	printf("Dimension of problem: %d\n",nd);
	printf("Number of PSO particles used: %d\n",np);
	printf("Maximum iterations per run allowed: %d\n",ni);
	printf("Maximum runs allowed : %d\n", maxRun);
	printf("Wall-clock time (Seconds): %lf\n",ellapse_time);

	printf("\n**========= Optimal Solution ======**\n");
	for(j=0;j<nd;j++){printf("%lf\n", gbest[j]);}
	}
	MPI_Finalize();
	free(lb);free(ub);free(gbest);
	return 0;
}
