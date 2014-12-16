/*
 H3_main_task1.c

Main program for task 1 in HP3b. To compile this main-program, make sure to change the makefile. 

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
	int i, j, m, n;
	double l;
	double d;
	double r;
	double phi;
	double x, y, dx;
	int max_grid_size, grid_size, grid_midpoint; // error and temp gridsize
	double error, it_error;
	double h_sq;

	// Initiation of variables
	error = 1.0; 
	l = 1;
	d = 0.2 * l;
	y = l / 2;
	max_grid_size = 21; // Maximal grid size used in the simulation
	grid_size = 21; // Smallest grid size: 11x11, next smallest grid size: 21x21 (Dynamic variable)
	grid_midpoint = (grid_size -1)/2;
	h_sq = pow((grid_midpoint+1)/l,2);

	// Declaration of arrays
	double** u; 
	double** res_error; //Error given by LAP(res_error) = res
	double** residual;
	double** temp;

	u = (double**) malloc(grid_size * sizeof(double*));
	res_error = (double**) malloc(grid_size * sizeof(double*));
	residual = (double**) malloc(grid_size * sizeof(double*));
	temp = (double**) malloc(grid_size * sizeof(double*));

	for(i = 0; i < grid_size; i++){
		u[i] = (double*) malloc(grid_size * sizeof(double));
		res_error[i] = (double*) malloc(grid_size * sizeof(double));
		residual[i] = (double*) malloc(grid_size * sizeof(double));
		temp[i] = (double*) malloc(grid_size * sizeof(double));
	}

	// Initiation of arrays
	for(i = 0; i < grid_size; i++){
		for(j = 0; j < grid_size; j++){
			u[i][j] = 0.0;
			
		}
	}



	/* TASK 1 */

	// File to save data
	FILE *file;
	file = fopen("phi.data","w");

	// Until error < 10^(-5)
	while(error >= pow(10,-5)){

		// Put the res and res_error to 0 at each iteration
		for(i = 0; i < grid_size; i++){
			for(j = 0; j < grid_size; j++){
				res_error[i][j] = 0.0;
				residual[i][j] = 0.0;
				temp[i][j] = 0.0;
			}
		}

		// Use Gauss-Seidel method to iterate three times
		for(i = 0; i < 3; i++){

			// Use Gauss-Seidel method, returns the error
			error = gauss_seidel(u, grid_size);

		}

		// Compute residual
		get_residual(u, residual, grid_size);

		// Restrict to coarser grid
		grid_size = decrease_grid(residual, grid_size);

		// Solve the residual equation exactly
		it_error = 1.0;		

		while(it_error >= 0.000001){
			it_error = get_error(residual, res_error, grid_size);
		}

		// Get fine grid
		grid_size = increase_grid(res_error, grid_size);
		
		// Interpolate
		for(i = 0; i < grid_size; i++){
			for(j = 0; j < grid_size; j++){
				u[i][j] += res_error[i][j];
			}
		}

		// Again, use Gauss-Seidel method to iterate three times
		for(i = 0; i < 3; i++){

			// Use Gauss-Seidel method, returns the error
			error = gauss_seidel(u, grid_size);

		}
	
	}

	// Print the final solution to a file
	for(i = 0; i < grid_size; i++){
		for(j = 0; j < grid_size; j++){
			fprintf(file, "%f \t", u[i][j]);
		}
		
		fprintf(file, "\n");
		
	}

	// Close file
	fclose(file);

	// Free allocated memory DOES NOT WORK
	for(i = 0; i < grid_size; i++){
			free(u[i]); 
			free(res_error[i]); 
			free(residual[i]); 
			free(temp[i]);
	}

	free(u); free(res_error); free(residual); free(temp);

	u = NULL; res_error = NULL; residual = NULL; temp = NULL;

}

