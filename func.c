/*
func.c
Contains functions for homeproblem 3/b

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.141592653589

// Function for the Gauss-Seidel method. Returns the maximal error and calculates the next iteration. 
double gauss_seidel(double **u_new, double **u_old, int grid_size){

	
	int i, j;
	int grid_midpoint = (grid_size -1)/2;	
	double error = 0.0;
	double h_sq = pow(1/(grid_size - 1),2);

	// Gauss-Seidel
	for(i = 1; i < grid_size - 1; i++){
		for(j = 1; j < grid_size - 1; j++){
			

			//If the point is on the dipole
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
			// Else, iterate GS-algorithm
			}else{ 

				u_new[i][j] = 0.25 * (u_old[i+1][j] + u_new[i-1][j] + u_old[i][j+1] + u_new[i][j-1]);
			}


			// Calculate maximal error
			if(fabs(u_new[i][j] - u_old[i][j]) > error){
				error = fabs(u_new[i][j] - u_old[i][j]);
			}
		}
	}

	return(error);
}

void get_residual(double **grid, double **res, int grid){

	int i,j;
	double h_sq = pow(1/(grid_size - 1),2);

	for(i = 1; i < grid_size - 1; i++){
		for(j = 1; j < grid_size - 1; j++){
		// -LAP(PHI)
		res[i][j] = 4*grid[i][j] - grid[i + 1][j] - grid[i - 1][j] -grid[i][j + 1] - grid[i][j - 1];
		res[i][j] *= h_sq;
		}
	}
	//Add the charge
	res[(grid_size -1)/2)*4/5][(grid_size -1)/2] += -h_sq;  
	res[(grid_size -1)/2)*6/5][(grid_size -1)/2] += h_sq;


}

