/*
func.c
Contains functions for homeproblem 3/b

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.141592653589

// Function for the Gauss-Seidel method. Returns the maximal error
double gauss_seidel(double **u_new, double **u_old, int grid_size, double error, double h_sq){

	
	int i, j;
	int grid_midpoint = (grid_size -1)/2;	
	error = 0.0;
	// Gauss-Seidel
	for(i = 1; i < grid_size - 1; i++){
		for(j = 1; j < grid_size - 1; j++){
			
			if(j == grid_midpoint){
				if(i == grid_midpoint*4/5){
					u_new[i][j] = -h_sq;
				}
				else if(i == grid_midpoint*6/5){
					u_new[i][j] = h_sq;
				}
				else{
					u_new[i][j] = 0.25 * (u_old[i+1][j] + u_new[i-1][j] + u_old[i][j+1] + u_new[i][j-1]);
			
				}
			}else{ 

				u_new[i][j] = 0.25 * (u_old[i+1][j] + u_new[i-1][j] + u_old[i][j+1] + u_new[i][j-1]);
			}


			// Calc maximal error
			if(fabs(u_new[i][j] - u_old[i][j]) > error){
				error = fabs(u_new[i][j] - u_old[i][j]);
			}
		}
	}
	return(error);

}

