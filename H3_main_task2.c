/*
 H3_main_task2.c

Main program for task 2 in HP3b. To compile this main-program, make sure to change the makefile. 

 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "func.h"
#define PI 3.141592653589

// Main program 
int main(){

	// Declaration of variables
	int i, j;
	double l;
	int max_grid_size, grid_size, grid_midpoint;
	double error, it_error;
	double h_sq;
	int gamma;
	int nbr_computations;

	// Initiation of variables
	error = 1.0; 
	l = 1.0;
	max_grid_size = 1281; // Maximal grid size used in the simulation
	grid_size = 41; // (Dynamic variable)
	grid_midpoint = (grid_size - 1)/2;
	h_sq = pow(l/(grid_size-1.0),2.0);
	gamma = 2;
	nbr_computations = 0;

	// Declaration of arrays
	double** u; 
	double** rho;

	u = (double**) malloc(max_grid_size * sizeof(double*));
	rho = (double**) malloc(max_grid_size * sizeof(double*));

	for(i = 0; i < max_grid_size; i++){
		u[i] = (double*) malloc(max_grid_size * sizeof(double));
		rho[i] = (double*) malloc(max_grid_size * sizeof(double));
	}

	// Initiation of arrays
	for(i = 0; i < max_grid_size; i++){
		for(j = 0; j < max_grid_size; j++){
			u[i][j] = 0.0;
			rho[i][j] = 0.0;			
		}
	}

	rho[grid_midpoint*4/5][grid_midpoint] = pow(h_sq,-1); //1/h^2
	rho[grid_midpoint*6/5][grid_midpoint] = pow(-h_sq,-1);


	// File to save data
	FILE *file;
	file = fopen("phi.data","w");

	// Call the multigrid function
	while (error >= pow(10,-5)){		
		error = multigrid(u, rho, grid_size, gamma, &nbr_computations);
		//printf("Error: %.10f \n", error);
	}
	
	printf("Computations: %i\n", nbr_computations);
	
	// Print the final solution to a file
	for(i = 0; i < grid_size; i++){
		for(j = 0; j < grid_size; j++){
			fprintf(file, "%f \t", u[i][j]);
		}
		
		fprintf(file, "\n");
		
	}

	// Close file
	fclose(file);

	// Free allocated memory
	for(i = 0; i < grid_size; i++){
			free(u[i]); 
			free(rho[i]); 
	}

	free(u); free(rho); 
	u = NULL; rho = NULL;
}

