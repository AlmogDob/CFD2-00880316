#include <stdio.h>
#define MATRIX2D_IMPLEMENTATION
#include "Matrix2D.h"

int offset2d(int i, int j, int ni, int nj) 
{
    assert(i < ni);
    assert(j < nj);

    return j * ni + i;
}

void print_array(double *a, int len)
{
    printf("--------------------\n");
    for (int i = 0; i < len; i++) {
        printf("%2d: %g \n", i, a[i]);
    }
    printf("\n");
    printf("--------------------\n");
}

#define ni 3
#define nj 3
int main(void)
{
    Mat2D m1;

    mat2D_alloc(&m1, nj, ni);
    MAT2D_AT(m1, 0, 0) = 1;
    MAT2D_AT(m1, 0, 1) = 2;
    MAT2D_AT(m1, 0, 2) = 3;
    MAT2D_AT(m1, 1, 0) = 6;
    MAT2D_AT(m1, 1, 1) = 7;
    MAT2D_AT(m1, 1, 2) = 8;

    // MAT2D_PRINT(m1);
    // print_array(m1.elements, ni*nj);
     
    /* ---------------------------------------- */
    double *a1 = (double *)calloc(ni*nj, sizeof(double));
    print_array(a1, ni*nj);

    a1[offset2d(0, 0, ni, nj)] = 1;
    a1[offset2d(1, 0, ni, nj)] = 2;
    a1[offset2d(2, 0, ni, nj)] = 3;
    a1[offset2d(0, 1, ni, nj)] = 4;
    a1[offset2d(1, 1, ni, nj)] = 5;
    a1[offset2d(2, 1, ni, nj)] = 6;
    a1[offset2d(0, 2, ni, nj)] = 7;
    a1[offset2d(1, 2, ni, nj)] = 8;
    a1[offset2d(2, 2, ni, nj)] = 9;
    print_array(a1, ni*nj);

    return 0;
}
