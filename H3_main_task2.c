/*
 H3_main_task2.c

Main program for task 2 in HP3b. To compile this main-program, make sure to change the makefile. 

Note that the code only functions properly for 
gird_size = 10*2^x +1 due to the Rho used.

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
	int max_grid_size; // Maximum grid size during the whole simulation 
	int grid_size; // Start grid size (Dynamic variable during the simulation)
	int grid_midpoint;
	double abs_diff, it_abs_diff;
	double h_sq;
	int gamma;
	int nbr_computations; // Calculated how many Gauss-Seidel computations that are done

	// Initiation of variables
	abs_diff = 1.0; 
	l = 1.0;
	max_grid_size = 81;
	grid_size = 81;
	grid_midpoint = (grid_size - 1)/2;
	h_sq = pow(l/(grid_size-1.0),2.0);
	gamma = 1;
	nbr_computations = 0;

	// Declaration of arrays
	double** u; // The potential
	double** rho; // The charge distribution

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

	//printf("Grid size: %i\n", grid_size);

	// Call the multigrid function
	while (abs_diff >= pow(10,-5)){		
		abs_diff = multigrid(u, rho, grid_size, gamma, &nbr_computations);
	}
	

	// Print in terminal
	printf("Grid size: %i Computations: %i\n", grid_size, nbr_computations);

	
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

