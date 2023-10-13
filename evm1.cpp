#include <stdio.h>
#include <time.h>
#include "matrix_solver.h"
#include "matrix_reader.h"

//function to create matrix 
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

void init_B(double *A, double *B, int n)
{
    double temp = 0;

    for (int i = 0; i < n; i++)
    {
        temp = 0;

        for (int j = 0; j < n; j += 2)
        {
            temp += A[i * n + j];
        }

        B[i] = temp;
    }
}

void print_mat(double *A, int n, int m, int p)
{
    for (int i = 0; i < n && i < p; i++)
    {
        for (int j = 0; j < m && j < p; j++)
        {
            printf("%10.3e ", A[i * m + j]);
        }

        printf("\n");
    }
}

double res_ctr(double *A, double *B, double *X, int n)
{
    double sum = 0.;
    double sum_str = 0.;
    double norm_b = 0.;
    
    for (int i = 0; i < n; i++)
    {
        sum_str = 0.;
        for (int j = 0; j < n; j++)
        {
            sum_str += A[i * n + j] * X[j];
        }
        
        sum += f_abs (sum_str - B[i]);
        norm_b += f_abs (B[i]);
    }
    
    return (norm_b > EPS ? sum / norm_b : 0);
}


int main(int argc, char *argv[])
{
    //A is a matrix, B is a right part 
    double *A = nullptr, *B = nullptr, *X = nullptr;

    //n is matrix dimentions, m is amount of outputs of matrix,
    //k is a formula number
    int n = 0, k = 0, m = 0, task = 11;

    //checking if input data are correct
    if (!((argc == 4 || argc == 5) 
        && (sscanf (argv[1], "%d", &n) == 1) && (n > 0)
        && (sscanf (argv[2], "%d", &m) == 1)
        && (sscanf (argv[3], "%d", &k) == 1)
        && ((k == 0 && argc == 5) || (k > 0 && k <= 4 && argc == 4))))
    {
        printf ("Usage: %s n m k [file]\n", argv[0]);
        return 0;
    }
    
    //creating matrix a, vector b and vector of solutions x
    A = new double[n * n];
    B = new double[n];
    X = new double[n];

    //memory allocation error
    if (!A || !B || !X)
    {
        if (A) delete [] A;
        if (B) delete [] B;
        if (X) delete [] X;
        printf("ERROR: Memory\n");
      
        return 0;
    }
    
    double t1 = 0., t2 = 0.;
    int res = 0;
    char *file = nullptr;
    
    //obtaining matrix from a file
    if (k)
        inint_A_form(A, n, k);
    else
    {
        file = argv[4];
        res = inint_A_file(A, n, file);
        if (res)
        {
            delete [] A;
            delete [] B;
            delete [] X;
            printf("ERROR: Wrong file\n");
            return 0;
        }
    }

    //initializing vector b
    init_B(A, B, n);

    //printing matrix a
    for (int i = 0; i < m; i++)
    {
        printf("===========");
    }

    printf("\nA = \n");
    print_mat (A, n, n, m);
    
    //printing vector b
    for (int i = 0; i < m; i++)
    {
        printf("===========");
    }    

    printf("\nB = \n");
    print_mat (B, 1, n, m);

    //solving matrix
    t1 = clock();
    res = solve(A, B, X, n);
    t1 = (clock() - t1) / CLOCKS_PER_SEC;
    
    //printing error if we can't find the solution
    if (res)
    {
        delete [] A;
        delete [] B;
        delete [] X;
        printf ("Impossible to find the solution\n");

        return 0;
    }

    //printing solution vector x
    for (int i = 0; i < m; i++)
    {
        printf("===========");
    }     

    printf("\nX = \n");
    print_mat(X, 1, n, m);
  
    if (k)
        inint_A_form (A, n, k);
    else
    {
        file = argv[4];
        res = inint_A_file(A, n, file);
        if (res)
        {
            delete [] A;
            delete [] B;
            delete [] X;
            printf("ERROR: Wrong file\n");
            return 0;
        }
    }    
   
    init_B(A, B, n);

    double resid = 0.;
    t2 = clock();
    resid = res_ctr(A, B, X, n);
    t2 = (clock() - t2) / CLOCKS_PER_SEC;
 
    printf ("%s : Task = %d Res = %e T1 = %.2f T2 = %.2f K = %d N = %d M = %d\n", argv[0], task, resid, t1, t2, k, n, m);

    delete [] A;
    delete [] B;
    delete [] X;

    return 0;
}
