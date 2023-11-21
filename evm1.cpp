#include <stdio.h>
#include <time.h>
#include <math.h>
#include <chrono>
#include "matrix_solver.h"
#include "matrix_reader.h"

//initializing vector B
void init_B(double *A, double *B, int n)
{
    double temp = 0;

    for (int i = 0; i < n; i++)
    {
        temp = 0;

        for (int j = 0; j < n; j += 2)
            temp += A[i * n + j];

        B[i] = temp;
    }
}

//printing any matrix
void print_mat(double *A, int n, int m, int p)
{
    for (int i = 0; i < n && i < p; i++)
    {
        for (int j = 0; j < m && j < p; j++)
            printf("%10.3e ", A[i * m + j]);

        printf("\n");
    }
}

double *count_norm_err(double *A, int n)
{
    double *X = nullptr;
    X = new double[n];

    double *vec = nullptr;
    vec = new double[n];

    double *Y = nullptr;
    Y = new double[n];

    for (int i = 0; i < n; i++)
        if (i % 2 == 0)
            X[i] = 1;
        else
            X[i] = 0;

    for (int i = 0; i < n; i++)
    {
        Y[i] = 0.;

        for (int j = 0; j < n; j++)
            Y[i] += A[i * n + j] * X[j];
    }
    
    int res = 0;
    res = parallel_solve(A, Y, vec, n);

    if (res)
    {
        delete [] X;
        delete [] Y;
        printf ("Impossible to find the solution\n");

        return nullptr;
    }

    for (int i = 0; i < n; i++)
        vec[i] = vec[i] - X[i];

    delete [] X;
    delete [] Y;

    return vec;
}

double norm1(double *vec, int n) 
{
    double res = 0;

    for(int i = 0; i < n; i++) 
        res += f_abs(vec[i]);
    
    return res;
}

double res_ctr1(double *A, double *B, double *X, int n)
{
    return (norm1(B, n) > EPS ? norm1(X, n) / norm1(B, n) : 0);
}

double norm2(double *vec, int n) 
{
    double res = 0;

    for(int i = 0; i < n; i++) 
        res += vec[i] * vec[i];
    
    return sqrt(res);
}

double res_ctr2(double *A, double *B, double *X, int n)
{
    return (norm2(B, n) > EPS ? norm2(X, n) / norm2(B, n) : 0);
}

double norminf(double *vec, int n) 
{
    double res = 0;

    for(int i = 0; i < n; i++) 
        if(f_abs(vec[i]) > res) 
            res = f_abs(vec[i]);

    return res;
}

double res_ctrinf(double *A, double *B, double *X, int n)
{
    return (norminf(B, n) > EPS ? norminf(X, n) / norminf(B, n) : 0);
}

int main(int argc, char *argv[])
{
    //A is a matrix, B is a right part 
    double *A = nullptr, *B = nullptr, *X = nullptr;

    //n is matrix dimentions, m is amount of outputs of matrix,
    //k is a formula number
    int n = 0, k = 0, m = 0;

    //checking if input data are correct
    if (!((argc == 4 || argc == 5) 
        && (sscanf(argv[1], "%d", &n) == 1) && (n > 0)
        && (sscanf(argv[2], "%d", &m) == 1)
        && (sscanf(argv[3], "%d", &k) == 1)
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
        printf("-----------");

    printf("\nA = \n");
    print_mat(A, n, n, m);
    
    //printing vector b
    for (int i = 0; i < m; i++)
        printf("-----------");

    printf("\nB = \n");
    print_mat(B, 1, n, m);

    //solving matrix
    auto start = std::chrono::high_resolution_clock::now();
    res = parallel_solve(A, B, X, n);
    auto end = std::chrono::high_resolution_clock::now();
    auto t1 = std::chrono::duration<double, std::milli>(end - start).count();
    
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
        printf("-----------");

    printf("\nX = \n");
    print_mat(X, 1, n, m);
  
    //initializing matrix again to find the error
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

    //counting all the resids and norms
    double resid1, resid2, residinf, nrmer1, nrmer2, nrmerinf = 0.;

    resid1 = res_ctr1(A, B, X, n);
    resid2 = res_ctr2(A, B, X, n);
    residinf = res_ctrinf(A, B, X, n); 

    double *vec = nullptr;
    vec = new double[n];
    vec = count_norm_err(A, n);

    nrmer1 = norm1(vec, n);
    nrmer2 = norm2(vec, n);
    nrmerinf = norminf(vec, n);

    //finding time to solve for single-threaded version
    //initializing matrix and b

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
    
    start = std::chrono::high_resolution_clock::now();
    res = solve(A, B, X, n);
    end = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::duration<double, std::milli>(end - start).count();
    
    //printing all the necessary stuff
    for (int i = 0; i < m; i++)
        printf("-----------");
    printf ("\n%s : \nNormErr1 = %e, NormErr2 = %e, NormErrinf = %e, \nRes1 = %e, Res2 = %e, Resinf = %e, \nTime parallel = %.5f, Time single = %.5f, \nK = %d, N = %d, M = %d\n", argv[0], nrmer1, nrmer2, nrmerinf, resid1, resid2, residinf, t1, t2, k, n, m);

    //cleaning up the memory
    delete [] A;
    delete [] B;
    delete [] X;

    return 0;
}
