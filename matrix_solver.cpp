#include "matrix_solver.h"
#include <thread>
#include <mutex>
using namespace std;

double f_abs(double a)
{
    return (a > 0 ? a : -a);
}

//finding the solution of the given triangular matrix
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

//preparing matrix to solve by making it triangular
int solve(double *A, double *B, double *X, int n)
{
    //maxx is a maximal element in the column
    //and ind_max is its index
    int ind_max = 0;
    double mnoz = 0., maxx = 0.;

    for (int steps = 0; steps < n; steps++)
    {
        //choosing maximal element and finding its index
        ind_max = steps;
        maxx = f_abs(A[steps * n + steps]);

        for (int i = steps; i < n; i++)
            //if elem > max it becomes max and changing 
            //index correspondingly
            if (f_abs(A[steps + i * n]) > maxx)
            {
                maxx = f_abs(A[steps + i * n]);
                ind_max = i;
            }

        //if maximal element is lower than constant return error
        if (maxx < EPS)
            return -1;

        //if choosen max element is not already what we're 
        //working with, changing its position to the current one
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
    
        //mnoz is equal to the element on the diagonal on the
        //steps row to the power of -1 and multiplying matrix
        //row and the answer to it creating nice 1-equal 
        //diagonal
        mnoz = 1. / A[steps * n + steps];
        B[steps] *= mnoz;
        
        for (int j = steps; j < n; j++)
            A[steps * n + j] *= mnoz;
       
        //removing all the column with this Xsteps below by
        //subtracting steps row from all the rows below 
        //mnoz times (mnoz becomes equal to Xsteps i for all
        //the is below Xsteps i'm so tired rn fuck
        for (int i = steps + 1; i < n; i++)
        {
            mnoz = A[i * n + steps];
            for (int j = steps; j < n; j++)
                A[i * n + j] -= mnoz * A[steps * n + j];
            
            B[i] -= mnoz * B[steps];
        }
    }
    
    //solving prepared matrix
    reverse_stroke (A, B, X, n);
    
    return 0;
}

void par_sub(double mnoz, double *A, double *B, int n, int i, int steps)
{ 
    mnoz = A[i * n + steps];
    for (int j = steps; j < n; j++)
        A[i * n + j] -= mnoz * A[steps * n + j];
            
    B[i] -= mnoz * B[steps];
}

void par_print()
{
    printf("parprint\n");
}

//preparing matrix to solve by making it triangular
int parallel_solve(double *A, double *B, double *X, int n)
{
    //maxx is a maximal element in the column
    //and ind_max is its index
    int ind_max = 0;
    double mnoz = 0., maxx = 0.;

    for (int steps = 0; steps < n; steps++)
    {
        //choosing maximal element and finding its index
        ind_max = steps;
        maxx = f_abs(A[steps * n + steps]);

        for (int i = steps; i < n; i++)
            //if elem > max it becomes max and changing 
            //index correspondingly
            if (f_abs(A[steps + i * n]) > maxx)
            {
                maxx = f_abs(A[steps + i * n]);
                ind_max = i;
            }

        //if maximal element is lower than constant return error
        if (maxx < EPS)
            return -1;

        //if choosen max element is not already what we're 
        //working with, changing its position to the current one
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
    
        //mnoz is equal to the element on the diagonal on the
        //steps row to the power of -1 and multiplying matrix
        //row and the answer to it creating nice 1-equal 
        //diagonal
        mnoz = 1. / A[steps * n + steps];
        B[steps] *= mnoz;
        
        for (int j = steps; j < n; j++)
            A[steps * n + j] *= mnoz;
       
        //removing all the column with this Xsteps below by
        //subtracting steps row from all the rows below 
        //mnoz times (mnoz becomes equal to Xsteps i for all
        //the is below Xsteps
        vector<thread> threads;
        for (int i = steps + 1; i < n; i++)
            threads.push_back(thread(par_sub, mnoz, A, B, n, i, steps)); 
            //threads.push_back(thread(par_print));

        for (int i = 0; i < threads.size(); ++i)
            threads[i].join();
    }
    
    //solving prepared matrix
    reverse_stroke (A, B, X, n);
    return 0;
}


