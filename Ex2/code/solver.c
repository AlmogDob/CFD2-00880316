#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#define MATRIX2D_IMPLEMENTATION
#include "Matrix2D.h"

#define MAXDIR 1000
#define MAXWORD MAXDIR
#define PI M_PI
#define dprintSTRING(expr) do{printf("%s:%d: ", __FILE__, __LINE__); printf(#expr " = %s\n", expr);} while(0)   /* macro for easy debuging*/
#define dprintINT(expr) do{printf("%s:%d: ", __FILE__, __LINE__); printf(#expr " = %d\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintFLOUT(expr) do{printf("%s:%d: ", __FILE__, __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintDOUBLE(expr) do{printf("%s:%d: ", __FILE__, __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintSCHEME(expr) do{ printf("%s:%d: scheme = ", __FILE__, __LINE__); if (flags & EXPLICIT_SW) { printf("EXPLICIT_SW\n"); } else if (flags & IMPLICIT_SW) { printf("IMPLICIT_SW\n"); } else if (flags & EXPLICIT_ROE) { printf("EXPLICIT_ROE\n"); } } while(0)     /* macro for easy debuging*/

typedef enum {
    EXPLICIT_SW  = (1 << 0),
    IMPLICIT_SW  = (1 << 1),
    EXPLICIT_ROE = (1 << 2),
} t_flag;

#ifdef __linux__
    #include <sys/stat.h>
    #include <dirent.h>
    #define ON_LINUX 1
    int create_empty_dir(char *parent_directory);
#endif

void read_input(char *input_file, int *N, double *x_min, double *x_max, double* L, double *delta_x, double *Re_inf, double *M_inf, double *CFL, double *gamma, double *T_inf, double *Pr_inf, int *max_iterations, double *final_time, t_flag *flags);
double calc_norm_mu(double norm_T, double gamma, double T_inf);
double calc_norm_kappa(double norm_T, double gamma, double T_inf);
double calc_norm_energy(double norm_p, double norm_rho, double norm_a, double gamma);
void calc_T_matrix(Mat2D T_matrix, double gamma, double norm_rho, double norm_a, double norm_u);

int main(int argc, char const *argv[])
{
/* declarations */
    char input_file[MAXDIR], output_dir[MAXDIR], temp_word[MAXWORD];
    int N, max_iterations;
    double x_min, x_max, L, delta_x, Re_inf, M_inf, CFL, gamma, T_inf, Pr_inf, final_time;
    t_flag flags = 0;

/* getting the input file and output file */
    if (--argc != 2 && argc != 4) {
        fprintf(stderr, "%s:%d: [Error] not right usage... Usage: main 'input file' 'output directory'\n", __FILE__, __LINE__);
        return -1;
    }

    strncpy(input_file, (*(++argv)), MAXDIR);

    if (input_file[MAXDIR-1] != '\0') {
        fprintf(stderr, "%s:%d: [Error] input too long\n", __FILE__, __LINE__);
        return -1;
    }

    strncpy(output_dir, (*(++argv)), MAXDIR);

    if (output_dir[MAXDIR-1] != '\0') {
        fprintf(stderr, "%s:%d: [Error] input too long\n", __FILE__, __LINE__);
        return -1;
    }

/*------------------------------------------------------------*/
/* reading the input */
    read_input(input_file, &N, &x_min, &x_max, &L, &delta_x, &Re_inf, &M_inf, &CFL, &gamma, &T_inf, &Pr_inf, &max_iterations, &final_time, &flags);

/* Checking the input */
    dprintINT(N);
    dprintDOUBLE(x_min);
    dprintDOUBLE(x_max);
    dprintDOUBLE(L);
    dprintDOUBLE(delta_x);
    dprintDOUBLE(Re_inf);
    dprintDOUBLE(M_inf);
    dprintDOUBLE(CFL);
    dprintDOUBLE(gamma);
    dprintDOUBLE(T_inf);
    dprintDOUBLE(Pr_inf);
    dprintDOUBLE(final_time);
    dprintINT(max_iterations);
    dprintSCHEME(flags);
    printf("--------------------\n");

    Mat2D mat = mat2D_alloc(3, 3);
    mat2D_fill(mat, 3);
    MAT2D_PRINT(mat);

/*------------------------------------------------------------*/
/* creating output directory */
    if (ON_LINUX) {
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
            return -1;
        }
    }

/*------------------------------------------------------------*/
/* allocating the matrices */

/*------------------------------------------------------------*/
/* initialization */

/*------------------------------------------------------------*/
/* the loop */
    
/*------------------------------------------------------------*/
/* output */

/*------------------------------------------------------------*/
/* freeing the memory */

    return 0;
}

#if ON_LINUX
/* create empty dir at 'parent directory'.
// if allready exists, delete all the files inside
returns 0 on success
this function handles the errors so on fail just quit
argument list:
parent_directory - char pointer to the directory name */
int create_empty_dir(char *parent_directory)
{
    char path_to_remove[BUFSIZ];

    if (mkdir(parent_directory, 0777) == -1) {
        if (errno == EEXIST) {
            DIR *dir = opendir(parent_directory);
            if (dir == NULL) {
                fprintf(stderr, "%s:%d: [Error] problem opening '%s': %s\n", __FILE__, __LINE__, parent_directory, strerror(errno));
                return 1;
            }
            struct dirent* entity;
            entity = readdir(dir);
            while (entity != NULL) {   /* directory is not empty */
                strncpy(path_to_remove, parent_directory, BUFSIZ);
                strncat(path_to_remove, "/", BUFSIZ/2);
                strncat(path_to_remove, entity->d_name, BUFSIZ);
                if (entity->d_type == DT_REG) {
                    if (remove(path_to_remove) != 0) {
                        fprintf(stderr, "%s:%d: [Error] problem removing '%s': %s\n", __FILE__, __LINE__, path_to_remove, strerror(errno));
                        return 1;
                    }
                }
                entity = readdir(dir);
            }
            closedir(dir);
            return 0;
        }

        fprintf(stderr, "%s:%d: [Error] problem making '%s': %s\n", __FILE__, __LINE__, parent_directory, strerror(errno));
        return 1;
    }
    return 0;
}
#endif

/* read input parameters from input file
argument list:
input-file     - file pointer to input file
N              - int pointer
x_min          - double pointer
x_max          - double pointer
Re_inf         - double pointer
M_inf          - double pointer
CFL            - double pointer
gamma          - double pointer
T_inf          - double pointer
Pr_inf         - double pointer
max_iterations - int pointer max number of desired iterations 
final_time     - double pointer to the final time of the program
flags          - bit flags */
void read_input(char *input_file, int *N, double *x_min, double *x_max, double *L, double *delta_x, double *Re_inf, double *M_inf, double *CFL, double *gamma, double *T_inf, double *Pr_inf, int *max_iterations, double *final_time, t_flag *flags)
{
    char current_word[MAXWORD];
    float temp_f;
    FILE *fp = fopen(input_file, "rt");
    if (fp == NULL) {
        fprintf(stderr, "%s:%d:[Error] problem opening input file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }

    while(fscanf(fp, "%s", current_word) != EOF) {
        if (!strcmp(current_word, "N")) {
            fscanf(fp, "%d", N);
        } else if (!strcmp(current_word, "x_max")) {
            fscanf(fp, "%g", &temp_f);
            *x_max = (double)temp_f;
        } else if (!strcmp(current_word, "x_min")) {
            fscanf(fp, "%g", &temp_f);
            *x_min = (double)temp_f;
        } else if (!strcmp(current_word, "Re_inf")) {
            fscanf(fp, "%g", &temp_f);
            *Re_inf = (double)temp_f;
        } else if (!strcmp(current_word, "M_inf")) {
            fscanf(fp, "%g", &temp_f);
            *M_inf = (double)temp_f;
        } else if (!strcmp(current_word, "CFL")) {
            fscanf(fp, "%g", &temp_f);
            *CFL = (double)temp_f;
        } else if (!strcmp(current_word, "gamma")) {
            fscanf(fp, "%g", &temp_f);
            *gamma = (double)temp_f;
        } else if (!strcmp(current_word, "T_inf")) {
            fscanf(fp, "%g", &temp_f);
            *T_inf = (double)temp_f;
        } else if (!strcmp(current_word, "Pr_inf")) {
            fscanf(fp, "%g", &temp_f);
            *Pr_inf = (double)temp_f;
        } else if (!strcmp(current_word, "iterations")) {
            fscanf(fp, "%g", &temp_f);
            *max_iterations = (int)temp_f;
        } else if (!strcmp(current_word, "final_time")) {
            fscanf(fp, "%g", &temp_f);
            *final_time = (double)temp_f;
        } else if (!strcmp(current_word, "scheme")) {
            fscanf(fp, "%s", current_word);
            if (!strcmp(current_word, "EXPLICIT_SW")) {
                *flags |= EXPLICIT_SW;
            } else if (!strcmp(current_word, "IMPLICIT_SW")) {
                *flags |= IMPLICIT_SW;
            } else if (!strcmp(current_word, "EXPLICIT_ROE")) {
                *flags |= EXPLICIT_ROE;
            }
        }
    }
    assert(*x_max > *x_min);
    *L = *x_max - * x_min;
    *delta_x = (*L) / (*N);

    fclose(fp);
}

double calc_norm_mu(double norm_T, double gamma, double T_inf)
{
    float constant = 110.4;

    return (sqrt(gamma * norm_T) * norm_T * (T_inf + constant)) / (gamma * T_inf * norm_T + constant);
}

double calc_norm_kappa(double norm_T, double gamma, double T_inf)
{
    float constant = 194;

    return (sqrt(gamma * norm_T) * norm_T * (T_inf + constant)) / (gamma * T_inf * norm_T + constant);
}

double calc_norm_energy(double norm_p, double norm_rho, double norm_a, double gamma)
{
    return (norm_p) / (gamma - 1) + 0.5 * norm_rho * norm_a * norm_a;
}

void calc_T_matrix(Mat2D T_matrix, double gamma, double norm_rho, double norm_a, double norm_u)
{
    MAT2D_AT(T_matrix, 1, 1) = 1;
    MAT2D_AT(T_matrix, 1, 2) = norm_rho / 2 / norm_a;
    MAT2D_AT(T_matrix, 1, 3) = - norm_rho / 2 / norm_a;
    MAT2D_AT(T_matrix, 2, 1) = norm_u;
    MAT2D_AT(T_matrix, 2, 2) = norm_rho / 2 / norm_a * (norm_u + norm_a);
    MAT2D_AT(T_matrix, 2, 3) = - norm_rho / 2 / norm_a * (norm_u - norm_a);
    MAT2D_AT(T_matrix, 3, 1) = norm_u * norm_u / 2;
    MAT2D_AT(T_matrix, 3, 2) = norm_rho / 2 / norm_a * (norm_u * norm_u / 2 + norm_u * norm_a + norm_a * norm_a / (gamma - 1));
    MAT2D_AT(T_matrix, 3, 3) = - norm_rho / 2 / norm_a * (norm_u * norm_u / 2 - norm_u * norm_a + norm_a * norm_a / (gamma - 1));
}

void calc_T_inverse_matrix(Mat2D T_inverse_matrix, double gamma, double norm_rho, double norm_a, double norm_u)
{
    MAT2D_AT(T_inverse_matrix, 1, 1) = 1;
    MAT2D_AT(T_inverse_matrix, 1, 2) = 1;
    MAT2D_AT(T_inverse_matrix, 1, 3) = 1;
    MAT2D_AT(T_inverse_matrix, 2, 1) = 1;
    MAT2D_AT(T_inverse_matrix, 2, 2) = 1;
    MAT2D_AT(T_inverse_matrix, 2, 3) = 1;
    MAT2D_AT(T_inverse_matrix, 3, 1) = 1;
    MAT2D_AT(T_inverse_matrix, 3, 2) = 1;
    MAT2D_AT(T_inverse_matrix, 3, 3) = 1;
}

