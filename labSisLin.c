#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main ()
{
    int n = 3;

    double *A = malloc(sizeof(double) * n*n);
    A[0] = 3;
    A[1] = 2;
    A[2] = -1;
    A[3] = 2;
    A[4] = -4;
    A[5] = 2;
    A[6] = -1;
    A[7] = 1;
    A[8] = 5;

    double *B = malloc(sizeof(double) * n);
    B[0] = 8;
    B[1] = -4;
    B[2] = 3;

    double *X = malloc(sizeof(double) * n);
    X[0] = 0;
    X[1] = 2;
    X[2] = -1;

    double erro, tIter, tTotal;
    tIter = 0;
    tTotal = 0;
    erro = jacobi(A, B, X, n, &tIter, &tTotal);

    printf("\nO sistema foi resolvido e as raízes são {%f, %f, %f}. O tempo médio das iterações foi de %.14f e o tempo total foi de %.14f\n", X[0], X[1], X[2], tIter, tTotal);

    tIter = 0;
    tTotal = 0;

    X[0] = 0;
    X[1] = 2;
    X[2] = -1;

    erro = gaussSeidel(A, B, X, n, &tIter, &tTotal);
    printf("\nO sistema foi resolvido e as raízes são {%f, %f, %f}. O tempo médio das iterações foi de %.14f e o tempo total foi de %.14f\n", X[0], X[1], X[2], tIter, tTotal);

    free(A);
    free(B);
    free(X);
}

