/*
 H3_main_task1.c

Main program for task 1 in HP3b. To compile this main-program, make sure to change the makefile. 

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
	int nbr_computations; // Calculated how many Gauss-Seidel computations that are done

	// Initiation of variables
	abs_diff = 1.0; 
	l = 1.0;
	max_grid_size = 1281; // Maximal grid size used in the simulation
	grid_size = 161; // Smallest grid size: 11x11, next smallest grid size: 21x21 (Dynamic variable)
	grid_midpoint = (grid_size -1.0)/2.0;
	h_sq = pow(l/(grid_size-1.0),2.0);
	nbr_computations = 0;

	// Declaration of arrays
	double** u; // Potential
	double** res_error; //Error given by LAP(res_error) = res
	double** residual;
	double** rho; // Charge distribution

	u = (double**) malloc(max_grid_size * sizeof(double*));
	res_error = (double**) malloc(max_grid_size * sizeof(double*));
	residual = (double**) malloc(max_grid_size * sizeof(double*));
	rho = (double**) malloc(max_grid_size * sizeof(double*));

	for(i = 0; i < max_grid_size; i++){
		u[i] = (double*) malloc(max_grid_size * sizeof(double));
		res_error[i] = (double*) malloc(max_grid_size * sizeof(double));
		residual[i] = (double*) malloc(max_grid_size * sizeof(double));
		rho[i] = (double*) malloc(max_grid_size * sizeof(double));
	}
	// Initiation of arrays
	for(i = 0; i < max_grid_size; i++){
		for(j = 0; j < max_grid_size; j++){
			u[i][j] = 0.0;
			rho[i][j] = 0.0;
		
		}
	}
	
	//Initiate the dipole
	rho[grid_midpoint*4/5][grid_midpoint] = pow(h_sq,-1); //1/h^2
	rho[grid_midpoint*6/5][grid_midpoint] = pow(-h_sq,-1);


	// File to save data
	FILE *file;
	file = fopen("phi.data","w");

	// Until abs_diff < 10^(-5)
	while(abs_diff >= pow(10,-6)){

		// Put the res and res_error to zero at each iteration
		for(i = 0; i < grid_size; i++){
			for(j = 0; j < grid_size; j++){
				residual[i][j] = 0.0;
				res_error[i][j] = 0.0;	
			}
		}

		// Iterate Gauss-Seidel relaxation three times
		for(i = 0; i < 3; i++){

			// Use Gauss-Seidel method, returns the abs_diff
			abs_diff = gauss_seidel(u, rho, grid_size, &nbr_computations);

		}
		
		// Compute residual
		get_residual(u, rho, residual, grid_size);


		// Reduce the gridsize of the residual
		grid_size = decrease_grid(residual, grid_size);

		// Solve the residual equation to 10^-5
		it_abs_diff = 1.0;		
		while(it_abs_diff >= pow(10,-5)){
			it_abs_diff = gauss_seidel(residual, res_error, grid_size, &nbr_computations);

		
		}

		// Increase the gridsize of thesolution to the residual-equation
		grid_size = increase_grid(res_error, grid_size);

		
		// Interpolate the solution to the residual-equation
		for(i = 0; i < grid_size; i++){
			for(j = 0; j < grid_size; j++){
				u[i][j] += res_error[i][j];		
			}
		}
		
		// Again, use Gauss-Seidel method to iterate three times
		for(i = 0; i < 3; i++){

			// Use Gauss-Seidel method, returns the abs_diff
			abs_diff = gauss_seidel(u, rho, grid_size, &nbr_computations);
		}

	}
	printf("Gridsize: %i Computations: %i\n",grid_size, nbr_computations);
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
	}
	free(u); free(res_error); free(residual); 
	u = NULL; res_error = NULL; residual = NULL; 
}

