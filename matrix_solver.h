#ifndef MATRIX_SOLVER_H_
#define MATRIX_SOLVER_H_

#define EPS 1e-14

double f_abs(double a);

void reverse_stroke(double *A, double *B, double *X, int n, int THREAD_NUMBER);

//solving matrix a with values b and writing results to x
int solve(double *A, double *B, double *X, int n, int THREAD_NUMBER);
int parallel_solve(double *A, double *B, double *X, int n, int THREAD_NUMBER);

#endif
