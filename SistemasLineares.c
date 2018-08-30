#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

void printX(double *X, int n) {
    for (int i = 0; i < n; i++) {
        printf("%f, ", X[i]);
    }  
    printf("\n");
}

double averageTime(double* times, int n) {
    double sum = 0;
    int i;

    for (i = 0; i < n; i++) {
        sum += times[i];
    }

    return (sum/n);
}

double calcErroAbsAprox(double *X, double *Xold, int n) {
    double max = 0;
    double error = 0;
    int i = 0;

    for (i = 0; i < n; i++) {
        error = fabs(X[i] - Xold[i]);
        /*printf("Erro %d:  %f\n", i, error);*/
        if (error > max) { 
            max = error;
        } 
    }

    return max;
}

// Método de Jacobi
double  jacobi (double  *A, double *B,  double *X,  int n,
		double *tIteracao, double *tTotal)
{
    /* Alocar Xold */
    double *Xold = malloc(sizeof(double) * n);

    /* Aloca vetor com os tempos das iterações */
    double *itTimes = malloc(sizeof(double) * 500);

    int itCount = 0;
    int i, j;
    double sum = 0;
    double error;  


    double tInicio = timestamp();

    while (itCount < MAXIT) {
        double tIterInicio = timestamp();

        memcpy(Xold, X, sizeof(double) * n);
        
        /* Rodar pro X passado (xold) */
        for (i = 0; i < n; i++) {
            sum = 0;

            /* Soma toda a linha i e multiplica por cada elemento do X */
            for (j = 0; j < n; j++) {
                sum += (A[i*n + j] * Xold[j]);

                /* if (i != j) { 
                    sum += (A[i*n + j] * X[i]);
                }*/
                /*printf("A[i*n + j] = %f, Xold[j] = %f\n", A[i*n + j], Xold[j]);*/
            }
            /*printf("B[%d] = %f, sum = %f\n", i, B[i], sum);*/

            /* Ao invés de executar um if para verificar se i = j, subtrai o elemento somado "a mais"  //todo ver se vale a pena */
            sum -= (A[i*n + i] * Xold[i]);

            /* Faz o passo da divisão pra cálculo do novo X */
            /*printf("B[%d] = %f, sum = %f\n", i, B[i], sum);*/
            X[i] = (B[i] - sum) / A[i*n + i];
        }

        error = calcErroAbsAprox(X, Xold, n);

        printf("Jacobi. Iteração %d: erro de %.14f\n", (itCount+1), error);

        /* Controle do tempo da iteração */
        itTimes[itCount] = timestamp() - tIterInicio;

        /*printf("x old: \n");
        printX(Xold, n);
        printf("x new: \n");
        printX(X, n);*/

        if (error <= EPS) {
            *tTotal = timestamp() - tInicio;
            *tIteracao = averageTime(itTimes, itCount);
            
            free(Xold);
            free(itTimes);
            /* Calcular tempo da iteração e colocar no tIteracao */
            return error;
        }
        itCount += 1;        
    }

    /* Caso estoure o número de iterações */

    *tTotal = timestamp() - tInicio;
    *tIteracao = averageTime(itTimes, itCount);

    free(Xold);
    free(itTimes);
    return -1;
}

// Método de Gauss-Seidel
double  gaussSeidel (double  *A, double *B,  double *X,  int n,
		     double *tIteracao, double *tTotal)
{
        /* Alocar Xold */
    double *Xold = malloc(sizeof(double) * n);

    /* Aloca vetor com os tempos das iterações */
    double *itTimes = malloc(sizeof(double) * 500);

    int itCount = 0;
    int i, j;
    double sum = 0;
    double error;  


    double tInicio = timestamp();

    while (itCount < MAXIT) {
        double tIterInicio = timestamp();

        memcpy(Xold, X, sizeof(double) * n);
        
        /* Rodar pro X passado (xold) */
        for (i = 0; i < n; i++) {
            sum = 0;

            /* Soma toda a linha i e multiplica por cada elemento do X */
            for (j = 0; j < n; j++) {
                sum += (A[i*n + j] * X[j]);

                /* if (i != j) { 
                    sum += (A[i*n + j] * X[i]);
                }*/
                /*printf("A[i*n + j] = %f, Xold[j] = %f\n", A[i*n + j], Xold[j]);*/
            }
            /*printf("B[%d] = %f, sum = %f\n", i, B[i], sum);*/

            /* Ao invés de executar um if para verificar se i = j, subtrai o elemento somado "a mais"  //todo ver se vale a pena */
            sum -= (A[i*n + i] * X[i]);

            /* Faz o passo da divisão pra cálculo do novo X */
            /*printf("B[%d] = %f, sum = %f\n", i, B[i], sum);*/
            X[i] = (B[i] - sum) / A[i*n + i];
        }

        error = calcErroAbsAprox(X, Xold, n);

        printf("Jacobi. Iteração %d: erro de %.14f\n", (itCount+1), error);

        /* Controle do tempo da iteração */
        itTimes[itCount] = timestamp() - tIterInicio;

        /*printf("x old: \n");
        printX(Xold, n);
        printf("x new: \n");
        printX(X, n);*/

        if (error <= EPS) {
            *tTotal = timestamp() - tInicio;
            *tIteracao = averageTime(itTimes, itCount);
            
            free(Xold);
            free(itTimes);
            /* Calcular tempo da iteração e colocar no tIteracao */
            return error;
        }
        itCount += 1;        
    }

    /* Caso estoure o número de iterações */

    *tTotal = timestamp() - tInicio;
    *tIteracao = averageTime(itTimes, itCount);

    free(Xold);
    free(itTimes);
    return -1;
}


