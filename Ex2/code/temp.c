#include <stdio.h>
#define MATRIX2D_IMPLEMENTATION
#include "Matrix2D.h"

int main(void)
{
    Mat2D m1, m2;

    mat2D_alloc(&m1, 3, 5);
    mat2D_alloc(&m2, 3, 1);
    
    mat2D_fill(m1, 2);
    MAT2D_AT(m1, 0, 0) = 7;
    MAT2D_AT(m1, 1, 1) = 7;
    MAT2D_AT(m1, 2, 2) = 7;

    mat2D_get_col(m2, m1, 0);

    MAT2D_PRINT(m1);
    MAT2D_PRINT(m2);

    return 0;
}
