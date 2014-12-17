/*
func.h for HP3b
*/

#ifndef _func_h
#define _func_h

extern double gauss_seidel(double **, double **,int);
extern void get_residual(double **, double **, double **,int);
extern int increase_grid(double **, int);
extern int decrease_grid(double **, int);
extern void multigrid(double **, double **, int , int);


#endif
