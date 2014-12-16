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

	// Initiation of variables
	error = 1.0; 
	l = 1.0;
	max_grid_size = 81; // Maximal grid size used in the simulation
	grid_size = 21; // (Dynamic variable)
	grid_midpoint = (grid_size - 1)/2;
	h_sq = pow((grid_midpoint + 1)/l,2);
	gamma = 2;

	// Declaration of arrays
	double** u; 
	double** res_error; //Error given by LAP(res_error) = res
	double** residual;
	double** temp;
	double** rho;

	u = (double**) malloc(max_grid_size * sizeof(double*));
	res_error = (double**) malloc(max_grid_size * sizeof(double*));
	residual = (double**) malloc(max_grid_size * sizeof(double*));
	temp = (double**) malloc(max_grid_size * sizeof(double*));
	rho = (double**) malloc(max_grid_size * sizeof(double*));

	for(i = 0; i < max_grid_size; i++){
		u[i] = (double*) malloc(max_grid_size * sizeof(double));
		res_error[i] = (double*) malloc(max_grid_size * sizeof(double));
		residual[i] = (double*) malloc(max_grid_size * sizeof(double));
		temp[i] = (double*) malloc(max_grid_size * sizeof(double));
		rho[i] = (double*) malloc(max_grid_size * sizeof(double));
	}

	// Initiation of arrays
	for(i = 0; i < max_grid_size; i++){
		for(j = 0; j < max_grid_size; j++){
			u[i][j] = 0.0;
			rho[i][j] = 0.0;			
		}
	}

	rho[grid_midpoint*4/5][grid_midpoint] = 1;
	rho[grid_midpoint*6/5][grid_midpoint] = -1;

	// File to save data
	FILE *file;
	file = fopen("phi.data","w");

	// Call the multigrid function
	multigrid(u, rho, grid_size, gamma);

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
			free(res_error[i]); 
			free(residual[i]); 
			free(temp[i]);
	}

	free(u); free(res_error); free(residual); free(temp);

	u = NULL; res_error = NULL; residual = NULL; temp = NULL;

}

