/*
func.c
Contains functions for homeproblem 3/b

Note that the code only functions properly for 
gird_size = N*10 + 1 

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.141592653589
#define MINGRID 11

// Function for the Gauss-Seidel method for LAP(A) = B. Returns the maximal error. 
double gauss_seidel(double **A, double **B, int grid_size){

	
	int i, j;

	int grid_midpoint = (grid_size - 1)/2;	
	
	double h_sq = pow(1.0/(grid_size - 1),2);
	double temp;
	double it_error = 0.0;

	// Gauss-Seidel
	for(i = 1; i < grid_size - 1; i++){
		for(j = 1; j < grid_size - 1; j++){
			
			temp = A[i][j];
			A[i][j] = 0.25 * (A[i+1][j] + A[i-1][j] + A[i][j+1] + A[i][j-1] - B[i][j]*h_sq);
			
			// Calculate maximal error
			if(fabs(A[i][j] - temp) > it_error){
				it_error = fabs(A[i][j] - temp);
			
				
			}
		}
	}

	return(it_error);
}
//Function that calculates the residual of LAP(A) = B
void get_residual(double **A, double **B , double **res, int grid_size){

	int i,j;
	double h_sq = pow(1.0/(grid_size - 1.0),2);
	
	

	for(i = 1; i < grid_size - 1; i++){
		for(j = 1; j < grid_size - 1; j++){
			res[i][j] = B[i][j] - (-4.0*A[i][j] + A[i-1][j] + A[i+1][j] + A[i][j-1] +A[i][j+1])/h_sq;
		}
	}		
}


// Function that increase the grid size, into a finer grid. Returns the new grid size.
int increase_grid(double **A, int grid_size){

	int i,j;
	int new_grid_size = 2 * (grid_size - 1) + 1;

	double** temp;
	temp = (double**) malloc(new_grid_size * sizeof(double*));

	for(i = 0; i < new_grid_size; i++){
		temp[i] = (double*) malloc(new_grid_size * sizeof(double));
	}



	// For all the outer points in the temp array, initialize with zeros
	for(i = 0; i < new_grid_size; i++){
		temp[0][i] = 0.0;
		temp[new_grid_size - 1][i] = 0.0;
		temp[i][0] = 0.0;
		temp[i][new_grid_size - 1] = 0.0;
	}

	// For all grid points that equals the old grid
	for(i = 0; i < grid_size - 1; i++){
		for(j = 0; j < grid_size - 1; j++){
			temp[2*i][2*j] = A[i][j];
			temp[2*i][2*j+1] = 0.5 * (A[i][j] + A[i][j+1]);
			temp[2*i+1][2*j] = 0.5 * (A[i+1][j] + A[i][j]); 
			temp[2*i+1][2*j+1] = 0.25 * (A[i][j] + A[i+1][j+1] + A[i][j+1] + A[i+1][j]);
		}
	}

	// Write the temp-array to the original array
	for(i = 0; i < new_grid_size; i++){
		for(j = 0; j < new_grid_size; j++){
			A[i][j] = temp[i][j];
		}
	}

	for(i = 0; i < new_grid_size; i++){
		free(temp[i]); 
	}

	free(temp);
	temp = NULL;	

	return new_grid_size;

	
}

// Function that decreases the grid size, into a coarser grid. Returns the new grid size.
int decrease_grid(double **A, int grid_size){

	int i,j;
	int new_grid_size = (int) (0.5 * (grid_size - 1) + 1);
	
	double** temp;
	temp = (double**) malloc(new_grid_size * sizeof(double*));

	for(i = 0; i < new_grid_size; i++){
		temp[i] = (double*) malloc(new_grid_size * sizeof(double));
	}
	

	// For all the outer points in the temp array, initialize with zeros
	for(i = 0; i < new_grid_size; i++){
		temp[0][i] = 0.0;
		temp[new_grid_size - 1][i] = 0.0;
		temp[i][0] = 0.0;
		temp[i][new_grid_size - 1] = 0.0;
	}

	// For all the inner grid points, with eight neighbours
	for(i = 1; i < new_grid_size - 1; i++){
		for(j = 1; j < new_grid_size - 1; j++){
			temp[i][j] = 0.25 * A[2*i+1][2*j+1] + (A[2*i][2*j] + A[2*i][2*j+2] + A[2*i+2][2*j] + A[2*i+2][2*j+2])/16.0 + (A[2*i+1][2*j+2] + A[2*i+1][2*j+2] + A[2*i][2*j+1] + A[2*i+2][2*j+1])/8.0;
		}
	}

	// Write the temp-array to the original array
	for(i = 0; i < new_grid_size; i++){
		for(j = 0; j < new_grid_size; j++){
			A[i][j] = temp[i][j];
		}
	}

	for(i = 0; i < new_grid_size; i++){
		free(temp[i]); 
	}

	free(temp);
	temp = NULL;	

	return new_grid_size;

}

// Solves the Poisson equation LAP(A) = B through the multigrid method

double multigrid(double **A, double **B, int grid_size, int gamma){
	
	double error = 1.0;
	int i,j;
	int n_smooth = 5;
	

	// If the most coarse grid, solve the equation using GS 
	if (grid_size <= MINGRID){
		while (error >= pow(10,-5)){		
			error = gauss_seidel(A, B, grid_size);
			return(error);
			
		
		}
				
	}else{ // If a finer grid than the most coarse run the MG


		// Declaration of arrays

		double** res;
		double** res_error;
		res = (double**) malloc(grid_size * sizeof(double*));
		res_error = (double**) malloc(grid_size * sizeof(double*));
		
		for(i = 0; i < grid_size; i++){
			res[i] = (double*) malloc(grid_size * sizeof(double));
			res_error[i] = (double*) malloc(grid_size * sizeof(double));
		}

		// Initiate arrays
		for(i = 0; i < grid_size; i++){
			for(j = 0; j < grid_size; j++){
				res[i][j] = 0.0;
				res_error[i][j] = 0.0;
			}
		}


		// Presmooth A
		for(i = 0; i < n_smooth; i++){
			error = gauss_seidel(A, B, grid_size);	
		}
		
		// Calculate the residual
		get_residual(A, B, res, grid_size);


			
		// Decrease the grid size of res
		grid_size = decrease_grid(res, grid_size);
	
		// Recursive solution to the residual equation
		for(i = 0; i < gamma; i++){		
			error = multigrid(res_error, res, grid_size, gamma);
			
		}
		
		// Increas res_error to original size of A 
		grid_size = increase_grid(res_error, grid_size);


		// Update A with the error
		for(i = 0; i < grid_size; i++){
			for(j = 0; j < grid_size; j++){
				A[i][j] += res_error[i][j];		
			}
		}
		
		// Postsmooth A
		for(i = 0; i < n_smooth; i++){
			error = gauss_seidel(A, B, grid_size);	
		}
		

		// Free allocated memory
		for(i = 0; i < grid_size; i++){
			free(res[i]); 
			free(res_error[i]); 
		}

		free(res); free(res_error);
		res = NULL; res_error = NULL;

		return(error);
	}


}

