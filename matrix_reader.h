#ifndef MATRIX_READER_H_
#define MATRIX_READER_H_

//function to create element of matrix 
double f(int i, int j, int n, int mode);

//reading matrix from a file
int inint_A_file(double *A, int n, char *file);

void inint_A_form(double *A, int n, int k);

#endif
