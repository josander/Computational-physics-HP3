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
	int max_grid_size, grid_size, grid_midpoint;
	double error;
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
	dx = 0.0001;
	max_grid_size = 21; // Maximal grid size used in the simulation
	grid_size = 21; // Smallest grid size: 11x11, next smallest grid size: 21x21 (Dynamic variable)
	grid_midpoint = (grid_size -1)/2;
	h_sq = pow((grid_midpoint+1)/l,2);

	// Declaration of arrays
	double** u; 
	double** rError; //Error given by LAP(rError) = res
	double** temp;

	u = (double**) malloc(grid_size * sizeof(double*));
	rError = (double**) malloc(grid_size * sizeof(double*));
	temp = (double**) malloc(grid_size * sizeof(double*));

	for(i = 0; i < grid_size; i++){
		u[i] = (double*) malloc(grid_size * sizeof(double));
		rError[i] = (double*) malloc(grid_size * sizeof(double));
		temp[i] = (double*) malloc(grid_size * sizeof(double));
	}

	// Initiation of arrays
	for(i = 0; i < grid_size; i++){
		for(j = 0; j < grid_size; j++){
			u[i][j] = 0.0;
			rError[i][j] = 0.0;

			temp[i][j] = 0.0;
		}
	}



	/* TASK 1 */

	// File to save data
	FILE *file;
	file = fopen("phi.data","w");

	// Until error < 10^(-5)
	//while(error >= pow(10,-5)){

		// Use Gauss-Seidel method to iterate three times
		for(i = 0; i < 3; i++){

			// Use Gauss-Seidel method, returns the error
			error = gauss_seidel(u, grid_size);


			// Change pointers

		}

		// Print the final solution to a file
		for(i = 0; i < grid_size; i++){
			for(j = 0; j < grid_size; j++){
				fprintf(file, "%f \t", u[i][j]);
			}
		
			fprintf(file, "\n");
		
		}


		// Restrict to coarser grid
		grid_size = decrease_grid(u, grid_size);
		grid_size = increase_grid(u, grid_size); // Only for tests


		// Solve the residual equation exactly

		// Interpolate

		// Again, use Gauss-Seidel method to iterate three times
		for(i = 0; i < 3; i++){

			// Use Gauss-Seidel method, returns the error
			error = gauss_seidel(u, grid_size);

			// Change pointers

		}

	//}

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

	u = NULL; rError = NULL; temp = NULL;

}

