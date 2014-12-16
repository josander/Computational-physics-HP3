/*
 H3_main.c
Main program
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
	int grid_size, grid_midpoint;
	double error;
	double h_sq;

	// Initiation of variables
	m_max = 10; // 10, 50, 100
	n_max = 10;
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
	grid_size = 11; //Smallest grid size: 11x11, next smallest grid size: 21x21
	grid_midpoint = (grid_size -1)/2;
	error = 1.0;
	h_sq = pow((grid_midpoint+1)/l,2);

	// Declaration of arrays
	double** u1; 
	double** u2;
	double** temp;
	u1 = (double**) malloc(grid_size * sizeof(double*));
	u2 = (double**) malloc(grid_size * sizeof(double*));
	temp = (double**) malloc(grid_size * sizeof(double*));

	for(i = 0; i < grid_size; i++){
		u1[i] = (double*) malloc(grid_size * sizeof(double));
		u2[i] = (double*) malloc(grid_size * sizeof(double));
		temp[i] = (double*) malloc(grid_size * sizeof(double));
	}

	// Initiation of arrays
	for(i = 0; i < grid_size; i++){
		for(j = 0; j < grid_size; j++){
			u1[i][j] = 0.0;
			u2[i][j] = 0.0;
			temp[i][j] = 0.0;
		}
	}

	// Initiate arrays with a dipole
	u1[grid_midpoint][grid_midpoint*4/5] = -h_sq;
	u1[grid_midpoint][grid_midpoint*6/5] = h_sq;	

	/* TASK 1 */

	// File to save data
	FILE *file;
	file = fopen("phi.data","w");

	// Use Gauss-Seidel method to iterate three times
	for(i = 0; i < 3; i++){

		// Use Gauss-Seidel method, returns the error
		error = gauss_seidel(u1, u2, grid_size, error, h_sq);

		// Change pointers
		temp = u1; 
		u1 = u2;
		u2 = temp;
	}

	// Print the final solution to a file
	for(i = 0; i < grid_size; i++){
		for(j = 0; j < grid_size; j++){
			fprintf(file, "%f \t", u2[i][j]);
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

	u1 = NULL; u2 = NULL; temp = NULL;

}

