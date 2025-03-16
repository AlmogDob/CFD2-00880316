/* This one-file library is heavily inspired by Tsoding's nn.h implementation of matrix
creation and operation. you can find the source code in:
https://github.com/tsoding/nn.h .
featured in this video of his:
https://youtu.be/L1TbWe8bVOc?list=PLpM-Dvs8t0VZPZKggcql-MmjaBdZKeDMw .*/


#ifndef MATRIX2D_H_
#define MATRIX2D_H_

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef MATRIX2D_MALLOC
#define MATRIX2D_MALLOC malloc
#endif //MATRIX2D_MALLOC

#ifndef MATRIX2D_ASSERT
#include <assert.h>
#define MATRIX2D_ASSERT assert
#endif //MATRIX2D_ASSERT

typedef struct {
    size_t rows;
    size_t cols;
    size_t stride;
    double *elements;
    
} Mat2D;

#define MAT2D_AT(m, i, j) (m).elements[mat2D_offset2d((m), (i), (j))]
#define MAT2D_PRINT(m) mat2D_print(m, #m, 0)

double rand_double(void);
void mat2D_alloc(Mat2D *m, size_t rows, size_t cols);
void mat2D_free(Mat2D m);
size_t mat2D_offset2d(Mat2D m, size_t i, size_t j);
void mat2D_fill(Mat2D m, double x);
void mat2D_rand(Mat2D m, double low, double high);
void mat2D_dot(Mat2D dst, Mat2D a, Mat2D b);
void mat2D_sum(Mat2D dst, Mat2D a);
void mat2D_mult(Mat2D m, double factor);
void mat2D_print(Mat2D m, const char *name, size_t padding);
void mat2D_identity_mat(Mat2D m);
void mat2D_copy(Mat2D res, Mat2D src);

#endif // MATRIX2D_H_

#ifdef MATRIX2D_IMPLEMENTATION

double rand_double(void)
{
    return (double) rand() / (double) RAND_MAX;
}

void mat2D_alloc(Mat2D *m, size_t rows, size_t cols)
{
    m->rows = rows;
    m->cols = cols;
    m->stride = cols;
    m->elements = (double*)MATRIX2D_MALLOC(sizeof(*(m->elements))*rows*cols);
    MATRIX2D_ASSERT(m->elements != NULL);
}

void mat2D_free(Mat2D m)
{
    free(m.elements);
}

size_t mat2D_offset2d(Mat2D m, size_t i, size_t j)
{
    assert(i < m.rows && j < m.cols);
    return i * m.cols + j;
}

void mat2D_fill(Mat2D m, double x)
{
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            MAT2D_AT(m, i, j) = x;
        }
    }
}

void mat2D_rand(Mat2D m, double low, double high)
{
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            MAT2D_AT(m, i, j) = rand_double()*(high - low) + low;
        }
    }
}

void mat2D_dot(Mat2D dst, Mat2D a, Mat2D b)
{
    MATRIX2D_ASSERT(a.cols == b.rows);
    size_t n = a.cols;
    MATRIX2D_ASSERT(a.rows == dst.rows);
    MATRIX2D_ASSERT(b.cols == dst.cols);

    for (size_t i = 0; i < dst.rows; i++) {
        for (size_t j = 0; j < dst.cols; j++) {
            for (size_t k = 0; k < n; k++) {
                MAT2D_AT(dst, i, j) += MAT2D_AT(a, i, k)*MAT2D_AT(b, k, j);
            }
        }
    }

}

void mat2D_sum(Mat2D dst, Mat2D a)
{
    MATRIX2D_ASSERT(dst.rows == a.rows);
    MATRIX2D_ASSERT(dst.cols == a.cols);
    for (size_t i = 0; i < dst.rows; ++i) {
        for (size_t j = 0; j < dst.cols; ++j) {
            MAT2D_AT(dst, i, j) += MAT2D_AT(a, i, j);
        }
    }
}

void mat2D_mult(Mat2D m, double factor)
{
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            MAT2D_AT(m, i, j) *= factor;
        }
    }
}

void mat2D_print(Mat2D m, const char *name, size_t padding)
{
    printf("%*s%s = [\n", (int) padding, "", name);
    for (size_t i = 0; i < m.rows; ++i) {
        printf("%*s    ", (int) padding, "");
        for (size_t j = 0; j < m.cols; ++j) {
            printf("%f ", MAT2D_AT(m, i, j));
        }
        printf("\n");
    }
    printf("%*s]\n", (int) padding, "");
}

void mat2D_identity_mat(Mat2D m)
{
    MATRIX2D_ASSERT(m.cols == m.rows);
    for (size_t i = 0; i < m.rows; ++i) {
        for (size_t j = 0; j < m.cols; ++j) {
            MAT2D_AT(m, i, j) = i == j ? 1 : 0;
            // if (i == j) {
            //     MAT2D_AT(m, i, j) = 1;
            // }
            // else {
            //     MAT2D_AT(m, i, j) = 0;
            // }
        }
    }
}

void mat2D_copy(Mat2D res, Mat2D src)
{
    MATRIX2D_ASSERT(res.cols == src.cols);
    MATRIX2D_ASSERT(res.rows == src.rows);

    for (size_t i = 0; i < res.rows; ++i) {
        for (size_t j = 0; j < res.cols; ++j) {
            MAT2D_AT(res, i, j) = MAT2D_AT(src, i, j);
        }
    }
    
}

#endif // MATRIX2D_IMPLEMENTATION