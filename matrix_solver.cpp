#include "matrix_solver.h"
#include <thread>
#include <mutex>
#include <vector>
#include <iostream>
using namespace std;

//const int THREAD_NUMBER = 8;

double f_abs(double a)
{
    return (a > 0 ? a : -a);
}

void par_back(double *A, double *B, double *X, int n, int start, int end)
{
    double src = 0.;

    for (int i = start; i >= end; i--)
    {
        src = B[i];
        
        for (int j = n - 1; j > i; j--)
            src -= A[i * n + j] * X[j];
        
        X[i] = src;
    }
}

//finding the solution of the given triangular matrix
void reverse_stroke(double *A, double *B, double *X, int n, int THREAD_NUMBER)
{
    const int GAP = n / THREAD_NUMBER;
    int start = n-1;
    int end = start;
    if (end - GAP < n)
            end -= GAP;
        else
            end = 0;

    vector<thread> threads;
    
    for (int i = 0; i < THREAD_NUMBER; i++)
    {
        threads.push_back(thread(par_back, A, B, X, n, start, end));
        start -= GAP;
        if (end - GAP < n)
            end -= GAP;
        else
            end = 0;
    }

    for (int i = 0; i < threads.size(); ++i)
        threads[i].join();

}

//preparing matrix to solve by making it triangular
int solve(double *A, double *B, double *X, int n, int THREAD_NUMBER)
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
        for (int i = steps + 1; i < n; i++)
        {
            mnoz = A[i * n + steps];
            for (int j = steps; j < n; j++)
                A[i * n + j] -= mnoz * A[steps * n + j];
            
            B[i] -= mnoz * B[steps];
        }
    }
    
    //solving prepared matrix
    reverse_stroke(A, B, X, n, THREAD_NUMBER);
    
    return 0;
}

void par_sub(double mnoz, double *A, double *B, int start, int cap, int n, int steps)
{   
    for (int i = start; i < cap; i++)
    {
        mnoz = A[i * n + steps];
        for (int j = steps; j < n; j++)
            A[i * n + j] -= mnoz * A[steps * n + j];
            
        B[i] -= mnoz * B[steps];
    }
}

void par_max(double tmp, int ind_max, double *A, double *B, int start, int cap, int n, int steps)
{
    for (int j = start; j < cap; j++)
    {
        tmp = A[ind_max * n + j];
        A[ind_max * n + j] = A[steps * n + j];
        A[steps * n + j] = tmp;
    }
}

//preparing matrix to solve by making it triangular
int parallel_solve(double *A, double *B, double *X, int n, int THREAD_NUMBER)
{
    //maxx is a maximal element in the column
    //and ind_max is its index
    int ind_max = 0;
    double maxx = 0.;

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

        int GAP = n / THREAD_NUMBER;
        int start = 0;
        int end = start;
          
        if (end + GAP < n)
                end += GAP;
            else
                end = n;

        vector<thread> threads_max;
                
        //if choosen max element is not already what we're 
        //working with, changing its position to the current one
        if (ind_max != steps)
        {
            double tmp = B[steps];
            B[steps] = B[ind_max];
            B[ind_max] = tmp;
            
            for (int i = 0; i < THREAD_NUMBER; i++)
            {
                threads_max.push_back(thread(par_max, tmp, ind_max, A, B, start, end, n, steps));
                start += GAP;
                if (end + GAP < n)
                    end += GAP;
                else
                    end = n;
            }

            for (int i = 0; i < threads_max.size(); ++i)
            threads_max[i].join();
        }
        
        //mult is equal to the element on the diagonal on the
        //steps row to the power of -1 and multiplying matrix
        //row and the answer to it creating nice 1-equal 
        //diagonal
        double mult = 0;

        mult = 1. / A[steps * n + steps];
        B[steps] *= mult;
        
        for (int j = steps; j < n; j++)
            A[steps * n + j] *= mult;
       
        //removing all the column with this Xsteps below by
        //subtracting steps row from all the rows below 
        //mult times (mult becomes equal to Xsteps i for all
        //the is below Xsteps
        
        GAP = (n - steps) / THREAD_NUMBER;

        start = steps+1;
        end = start;
        if (end + GAP < n)
                    end += GAP;
                else
                    end = n;

        vector<thread> threads;
        for (int i = 0; i < THREAD_NUMBER; i++)
        {
            threads.push_back(thread(par_sub, mult, A, B, start, end, n, steps));
            start += GAP;
            if (end + GAP < n)
                end += GAP;
            else
                end = n;
        }

        for (int i = 0; i < threads.size(); ++i)
            threads[i].join();
    }
    
    //solving prepared matrix
    reverse_stroke (A, B, X, n, THREAD_NUMBER);
    return 0;
}
