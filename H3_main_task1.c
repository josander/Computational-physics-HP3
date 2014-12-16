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
	int m_max, n_max;
	double lambda;
	double epsilon0;
	double l, l_inv;
	double d;
	double r, r_c, r_plus, r_minus;
	double phi;
	double x, y, dx;
	int max_grid_size, grid_size, grid_midpoint; // error and temp gridsize
	double error, itError;
	double h_sq;

	// Initiation of variables
	m_max = 10; // 10, 50, 100
	n_max = 10;
	error = 1.0; 
	lambda = 1;
	epsilon0 = 1;
	l = 1;
	l_inv = 1/l;
	d = 0.2 * l;
	r_c = l/2;
	r_plus = r_c + d / 2.0;
	r_minus = r_c - d / 2.0;
	y = l / 2;
	max_grid_size = 21; // Maximal grid size used in the simulation
	grid_size = 21; // Smallest grid size: 11x11, next smallest grid size: 21x21 (Dynamic variable)
	grid_midpoint = (grid_size -1)/2;
	h_sq = pow((grid_midpoint+1)/l,2);

	// Declaration of arrays
	double** u; 
	double** rError; //Error given by LAP(rError) = res
	double** residual;
	double** temp;

	u = (double**) malloc(grid_size * sizeof(double*));
	rError = (double**) malloc(grid_size * sizeof(double*));
	residual = (double**) malloc(grid_size * sizeof(double*));
	temp = (double**) malloc(grid_size * sizeof(double*));

	for(i = 0; i < grid_size; i++){
		u[i] = (double*) malloc(grid_size * sizeof(double));
		rError[i] = (double*) malloc(grid_size * sizeof(double));
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

		// Put the res and rError to 0 at each iteration
		for(i = 0; i < grid_size; i++){
			for(j = 0; j < grid_size; j++){
				rError[i][j] = 0.0;
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
		itError = 1.0;		

		
		while(itError > 0.00000001){
			itError = get_error(residual, rError, grid_size);
				//printf("%.10f \n", itError );

	
		grid_size = increase_grid(rError, grid_size);
		
		// Interpolate
		for(i = 0; i < grid_size; i++){
			for(j = 0; j < grid_size; j++){
				u[i][j] += rError[i][j];
				
			}
		}

		// Again, use Gauss-Seidel method to iterate three times
		for(i = 0; i < 3; i++){

			// Use Gauss-Seidel method, returns the error
			error = gauss_seidel(u, grid_size);



		}
		//printf("%f \n", error );
	
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
	/*for(i = 0; i < grid_size; i++){
			free(u1[i]); 
			free(u2[i]); 
			free(temp[i]);
	}

	free(u1); free(u2); free(temp);*/

	u = NULL; rError = NULL; residual = NULL; 
}}

