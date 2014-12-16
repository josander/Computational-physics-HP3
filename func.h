/*
func.h for HP3b
*/

#ifndef _func_h
#define _func_h

extern double gauss_seidel(double **, double **, int, double);
extern double get_residual();
extern int increase_grid(double **, int);
extern int decrease_grid(double **, int);

#endif
