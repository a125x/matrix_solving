#include <stdio.h>
#include "matrix_reader.h"

//reading matrix from a file
int inint_A_file(double *A, int n, char *file)
{
    FILE *in = fopen(file, "r");
    int i = 0, j = 0;

    if (!in)
      return 1;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
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

