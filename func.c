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

// Function for the Gauss-Seidel method. Returns the maximal error and calculates the next iteration. 
double gauss_seidel(double **u, int grid_size){

	
	int i, j;
	int grid_midpoint = (grid_size -1)/2;	
	
	double h_inv_sq = pow(1/(grid_size - 1),2);
	double temp;
	double itError = 0.0;

	// Gauss-Seidel
	for(i = 1; i < grid_size - 1; i++){
		for(j = 1; j < grid_size - 1; j++){
			
			temp = u[i][j];
			//If the point is on the dipole
			if(j == grid_midpoint){
				if(i == grid_midpoint*4/5){
					u[i][j] = 0.25 * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) + 1; //h_inv_sq * h_sq = 1
				}
				else if(i == grid_midpoint*6/5){
					u[i][j] = 0.25 * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) - 1;
				}
				else{
					u[i][j] = 0.25 * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]);
			
				}
			//  GS-algorithm
			}else{ 

				u[i][j] = 0.25 * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]);
			}


			// Calculate maximal error
			if(fabs(u[i][j] - temp) > itError){
				itError = fabs(u[i][j] - temp);
			}
		}
	}

	return(itError);
}
//Function that calculates the residual
void get_residual(double **grid, double **res, int grid_size){

	int i,j;
	double h_inv_sq = pow(1.0/(grid_size - 1),2);
	for(i = 1; i < grid_size - 1; i++){
		for(j = 1; j < grid_size - 1; j++){
		// -LAP(PHI)
			res[i][j] = 4*grid[i][j] - grid[i + 1][j] - grid[i - 1][j] -grid[i][j + 1] - grid[i][j - 1];
				
			res[i][j] *= h_inv_sq;

		}

	}
	
	//Add the chargedist 
	res[((grid_size -1)/2)*4/5][(grid_size -1)/2] += -h_inv_sq;  
	res[((grid_size -1)/2)*6/5][(grid_size -1)/2] += h_inv_sq;


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
		//printf("\n");
	}

			//printf("%f\t", temp[15][1]);

	// Write the temp-array to the original array
	for(i = 0; i < new_grid_size; i++){
		for(j = 0; j < new_grid_size; j++){
		//	printf("%f ", temp[i][j]);

			A[i][j] = temp[i][j];


		}
		//printf("\n");
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

