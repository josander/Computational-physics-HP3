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

// Function for the Gauss-Seidel method. Returns the maximal error and calculates the next iteration. 
double gauss_seidel(double **A, double **B, int grid_size){

	
	int i, j;
	int grid_midpoint = (grid_size -1)/2;	
	
	double h_inv_sq = pow(1/(grid_size - 1),2);
	double temp;
	double itError = 0.0;

	// Gauss-Seidel
	for(i = 1; i < grid_size - 1; i++){
		for(j = 1; j < grid_size - 1; j++){
			
			temp = u[i][j];
			A[i][j] = 0.25 * (A[i+1][j] + A[i-1][j] + A[i][j+1] + A[i][j-1]) - h_inv_sq*B[i][j];
			


			// Calculate maximal error
			if(fabs(u[i][j] - temp) > itError){
				itError = fabs(u[i][j] - temp);
			}
		}
	}

	return(itError);
}
//Function that calculates the residual
void get_residual(double **A, double **B , double **res,int grid_size){

	int i,j;
	double h_inv_sq = pow(1.0/(grid_size - 1),2);
	for(i = 1; i < grid_size - 1; i++){
		for(j = 1; j < grid_size - 1; j++){
		// -LAP(PHI)
			res[i][j] = 4*A[i][j] - A[i + 1][j] - A[i - 1][j] -A[i][j + 1] - A[i][j - 1];
				
			res[i][j] *= h_inv_sq;
			res[i][j] += B[i][j]

		}

	}
	

}
// Function that iterates GS to solve the R-eq  LAP(E) = R 
double get_error(double **res, double **error, int grid_size){

	int i,j;
	double h_sq = pow(1.0/(grid_size - 1),2);	
	double h_inv_sq = pow((grid_size - 1),2);
	double temp;
	double itError = 0.0; // NOTE: iteration error, not the solution to the R-eq
	//GS
	for(i = 1; i < grid_size - 1; i++){
		for(j = 1; j < grid_size - 1; j++){
			temp = error[i][j];
			error[i][j] = 0.25 * (error[i+1][j] + error[i-1][j] + error[i][j+1] + error[i][j-1] - h_sq*res[i][j] );

			// Calculate maximal error
			if(fabs(error[i][j] - temp) > itError){
				itError = fabs(error[i][j] - temp);
		
			}		
		}

	}

	
return(itError);


}

// Function that increase the grid size, into a finer grid. Returns the new grid size.
int increase_grid(double **A, int grid_size){

	int i,j;
	int new_grid_size = 2 * (grid_size - 1) + 1;
	double temp[new_grid_size][new_grid_size];

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

	return new_grid_size;
}

// Function that decreases the grid size, into a coarser grid. Returns the new grid size.
int decrease_grid(double **A, int grid_size){

	int i,j;
	int new_grid_size = (int) (0.5 * (grid_size - 1) + 1);
	double temp[new_grid_size][new_grid_size];
	

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
			temp[i][j] = 0.25 * A[2*i+1][2*j+1] + (A[2*i][2*j] + A[2*i][2*j+2] + A[2*i+2][2*j] + A[2*i+2][2*j+2])/16 + (A[2*i+1][2*j+2] + A[2*i+1][2*j+2] + A[2*i][2*j+1] + A[2*i+2][2*j+1])/8;
		}
	}

	// Write the temp-array to the original array
	for(i = 0; i < new_grid_size; i++){
		for(j = 0; j < new_grid_size; j++){
			A[i][j] = temp[i][j];
		}
	}

	return new_grid_size;

}

// solves LAP(A) = B 
void multigrid(double **A, double **B, int grid_size, int gamma){
	
	double error = 1.0;
	int i,j;
	int n_smooth = 3;

	// Initiate residual and the soultion to the r-eq
	double** res;
	double** res_error;

	res = (double**) malloc(grid_size * sizeof(double*));
	res_error = (double**) malloc(grid_size * sizeof(double*));
	for(i = 0; i < grid_size; i++){
		res[i] = (double*) malloc(grid_size * sizeof(double));
		res_error[i] = (double*) malloc(grid_size * sizeof(double));

	}
		
	for(i = 0; i < grid_size; i++){
		for(j = 0; j < grid_size; j++){
			res[i][j] = 0.0;
			res_error[i][j] = 0.0;
			
		}
	}


	
	if (grid_size == MINGRID){
		while (error >= pow(10,-5)){		
			error = gauss_seidel(A, B, grid_size);
		}
	}else{
		// Presmooth A
		for(i = 0; i < n_smooth; i++){
			error = gauss_seidel(A, B, grid_size);	
		}
		// Calculate the residual
		get_residual(A, res ,B ,grid_size);
		// Decrase the grid_size of res
		grid_size = decrease_grid(res, grid_size);
		
		// Recursive solution to the residual equation
		for(i = 0; i < gamma; i++){
			multigrid(res_error, res, grid_size, gamma);
		}

		//Increas res_error to original size of A 
		grid_size = increase_grid(res_error, grid_size);

		// Update A
		for(i = 0; i < grid_size; i++){
			for(j = 0; j < grid_size; j++){
				A[i][j] += res_error[i][j];
			}
		}
		// Postsmooth of A
		for(i = 0; i < n_smooth; i++){
			error = gauss_seidel(A, B, grid_size);	
		}

		

	}

	for(i = 0; i < grid_size; i++){
		free(res[i]); 
		free(res_error[i]); 
	}

	free(res); free(res_error);
	res = NULL; res_error = NULL;
}

