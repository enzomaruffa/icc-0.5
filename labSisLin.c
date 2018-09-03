#include <stdio.h>
#include <math.h>
#include <string.h>

#include "utils.h"
#include "SistemasLineares.h"

int main ()
{
    /* Lê o n */
    int n = 0;
    printf("Digite o número de incógnitas \n");
    scanf("%d", &n);

    int i = 0;
    double *A = malloc(sizeof(double) * n*n);
    for (i = 0; i < n*n; i++) {
        printf("Digite o termo A[%d][%d]:\n", (int)(i/n), i%n);
        scanf("%lf", &(A[i]));
    }

    double *B = malloc(sizeof(double) * n);
    for (i = 0; i < n; i++) {
        printf("Digite o termo B[%d]:\n", i);
        scanf("%lf", &(B[i]));
    }
    
    /* Imprime a matriz lida */
    printf("\nCONFERINDO!");
    printf("\nMatriz de dimensão: %d\n", n);
    printf("Coeficientes da matriz: ");
    for (i = 0; i < n*n; i++) {
        if (i%n == 0) {
            printf("\n");
        }
        printf("%f ", A[i]);
    }
    printf("\nTermos de B: ");
    for (i = 0; i < n; i++) {
        printf("%f ", B[i]);
    }

    printf("\n");
    double *X = malloc(sizeof(double) * n);
    for (i = 0; i < n; i++) {
        printf("Digite o chute X[%d]:\n", i);
        scanf("%lf", &(X[i]));
    }

    printf("\nChutes: ");
    for (i = 0; i < n; i++) {
        printf("%f ", X[i]);
    }


    double erro, tIter, tTotal;

    /* Executando o método de Jacobi */
    printf("\n---------------------------------\n");
    printf("Executando o método de Jacobi\n");
    double *XJacobi = malloc(sizeof(double) * n);
    memcpy(XJacobi, X, sizeof(double) * n);
    
    tIter = 0;
    tTotal = 0;

    erro = jacobi(A, B, XJacobi, n, &tIter, &tTotal);

    printf("\nO sistema foi resolvido. O tempo médio das iterações foi de %.14f e o tempo total foi de %.14f\n", tIter, tTotal);
    printf("Raízes encontradas: { ");
    for (i = 0; i < n; i++) {
        printf("%f ", XJacobi[i]);
    }
    printf("} \n");


    /* Executando o método de Gauss */
    printf("\n---------------------------------\n");
    printf("Executando o método de Gauss-Seidel\n");
    double *XGauss = malloc(sizeof(double) * n);
    memcpy(XGauss, X, sizeof(double) * n);

    tIter = 0;
    tTotal = 0;

    erro = gaussSeidel(A, B, XGauss, n, &tIter, &tTotal);
    printf("\nO sistema foi resolvido. O tempo médio das iterações foi de %.14f e o tempo total foi de %.14f\n", tIter, tTotal);
    printf("Raízes encontradas: { ");
    for (i = 0; i < n; i++) {
        printf("%f ", XGauss[i]);
    }
    printf("} \n");


    /* Fim */

    free(XJacobi);
    free(XGauss);

    free(A);
    free(B);
    free(X);
}

