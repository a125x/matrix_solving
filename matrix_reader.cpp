#include <stdio.h>
#include "matrix_reader.h"

//function to create element of matrix 
double f(int i, int j, int n, int mode)
{
    switch (mode)
    {
      case 1: 
          return n - (i > j ? i : j);
      case 2: 
          return (i > j ? i : j) + 1;
      case 3: 
          return (i > j ? i - j : j - i);
      case 4:
          return 1. / (i + j + 1);
      default:
          return 0.;
    }

    return 0.;
}

//initializing matrix A
void inint_A_form(double *A, int n, int k)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
          {
            A[i * n + j] = f(i, j, n, k);
          }
    }
}

//reading matrix from a file
int inint_A_file(double *A, int n, char *file)
{
    FILE *in = fopen(file, "r");

    if (!in)
      return 1;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (fscanf(in, "%lf", A + i * n + j) != 1)
            {
                fclose(in);
                return 2;
            }
        }
    }

    fclose (in);
    return 0;
}

