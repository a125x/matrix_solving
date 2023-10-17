#include "matrix_solver.h"

double f_abs(double a)
{
    return (a > 0 ? a : -a);
}

void reverse_stroke(double *A, double *B, double *X, int n)
{
    double src = 0.;

    for (int i = n - 1; i >= 0; i--)
    {
        src = B[i];

        for (int j = n - 1; j > i; j--)
            src -= A[i * n + j] * X[j];

        X[i] = src;
    }
}

//solving matrix a with values b and writing results to x
int solve(double *A, double *B, double *X, int n)
{
    int ind_max = 0;
    double mnoz = 0., maxx = 0.;

    for (int steps = 0; steps < n; steps++)
    {
        ind_max = steps;
        maxx = f_abs(A[steps * n + steps]);

        for (int i = steps; i < n; i++)
            if (f_abs(A[steps + i * n]) > maxx)
            {
                maxx = f_abs(A[steps + i * n]);
                ind_max = i;
            }

        if (maxx < EPS)
            return -1;

        if (ind_max != steps)
        {
            double tmp = B[steps];
            B[steps] = B[ind_max];
            B[ind_max] = tmp;
            for (int j = steps; j < n; j++)
            {
                tmp = A[ind_max * n + j];
                A[ind_max * n + j] = A[steps * n + j];
                A[steps * n + j] = tmp;
            }
        }
     
        mnoz = 1. / A[steps * n + steps];
        B[steps] *= mnoz;
        
        for (int j = steps; j < n; j++)
            A[steps * n + j] *= mnoz;
        
        for (int i = steps + 1; i < n; i++)
        {
            mnoz = A[i * n + steps];
            for (int j = steps; j < n; j++)
                A[i * n + j] -= mnoz * A[steps * n + j];
            
            B[i] -= mnoz * B[steps];
        }
    }
    
    reverse_stroke (A, B, X, n);
    
    return 0;
}
