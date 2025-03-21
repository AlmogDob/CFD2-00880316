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
#define dprintSIZE_T(expr) do{printf("%s:%d: ", __FILE__, __LINE__); printf(#expr " = %zu\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintSCHEME(expr) do{ printf("%s:%d: scheme = ", __FILE__, __LINE__); if (flags & EXPLICIT_SW) { printf("EXPLICIT_SW -> "); } else if (flags & IMPLICIT_SW) { printf("IMPLICIT_SW -> "); } else if (flags & EXPLICIT_ROE) { printf("EXPLICIT_ROE -> "); } if (flags & CONST_DELTA_TIME) { printf("CONST DELTA TIME\n"); } else { printf("VARIABLE DELTA TIME\n"); }} while(0)     /* macro for easy debuging*/

typedef enum {
    EXPLICIT_SW      = (1 << 0),
    IMPLICIT_SW      = (1 << 1),
    EXPLICIT_ROE     = (1 << 2),
    CONST_DELTA_TIME = (1 << 3),
} t_flag;

#ifdef __linux__
    #include <sys/stat.h>
    #define __USE_MISC
    #include <dirent.h>
    #define ON_LINUX 1
    int create_empty_dir(char *parent_directory);
#endif

void read_input(char *input_file, int *N, double *x_min, double *x_max, double *L, double *norm_delta_x, double *Re_inf, double *M_inf, double *CFL, double *gamma, double *Pr_inf, double *R_specific, double *c_v, double *T_inf, int *max_iterations, double *const_norm_delta_time, double *final_time, t_flag *flags);
void output_metadata(char *output_dir, int N, double x_min, double x_max, double Re_inf, double M_inf, double CFL, double gamma, double Pr_inf, double R_specific, double T_inf, int max_iterations, double final_time, t_flag flags);
void print_mat2D_to_file(FILE *fp, Mat2D m);
int offset3d(int i, int j, int k, int cols, int rows, int layers);
double calc_norm_mu(double norm_T, double gamma, double T_inf);
double calc_norm_kappa(double norm_T, double gamma, double T_inf);
double calc_norm_energy(double gamma, double norm_rho, double norm_u, double norm_p);
double calc_norm_p(double gamma, double norm_rho, double norm_u, double norm_e);
double calc_norm_total_enthalpy_at_i(Mat2D Q, double gamma, int i);
double calc_norm_T(Mat2D Q, double gamma, int i);
double calc_norm_delta_t(Mat2D Q, double gamma, double R_specific, double T_inf, double norm_delta_x, double const_norm_delta_time, double CFL, int N, t_flag flags);
void calc_T_matrix_at_i(Mat2D T_matrix, Mat2D Q, double gamma, int i);
void calc_T_inverse_matrix_at_i(Mat2D T_inverse_matrix, Mat2D Q, double gamma, int i);
void calc_lambda_plus_matrix_at_i(Mat2D lambda_plus_matrix, Mat2D Q, double gamma, int i, double epsilon);
void calc_lambda_minus_matrix_at_i(Mat2D lambda_minus_matrix, Mat2D Q, double gamma, int i, double epsilon);
void calc_A_plus_or_minus_at_i(Mat2D Q, Mat2D norm_A, Mat2D work_3_3_mat1, Mat2D work_3_3_mat2, Mat2D work_3_3_mat3, Mat2D work_3_3_mat5, double gamma, double epsilon, char sign, int i);
void calc_norm_P_minus_Rx_matrix_at_i(Mat2D norm_P_minus_Rx, Mat2D Q, double norm_delta_x, double gamma, double T_inf, int i);
void calc_norm_R_matrix_at_i(Mat2D norm_R, Mat2D Q, double gamma, double T_inf, double Pr_inf, double c_v, int i);
double initialize_Q(char *init_conditions_file, Mat2D Q, double gamma, int N);
void apply_BC(Mat2D Q, int N);
void calc_vector_of_E(Mat2D E, Mat2D Q, double gamma, int N);
void calc_vector_of_tilde_norm_E_at_half_sw(Mat2D tilde_norm_E, Mat2D Q, Mat2D work_3_3_mat1, Mat2D work_3_3_mat2, Mat2D work_3_3_mat3, Mat2D work_3_3_mat4, Mat2D work_3_3_mat5, Mat2D work_3_1_mat1, Mat2D work_3_1_mat2, Mat2D work_3_1_mat3, double gamma, double epsilon, int N);
void calc_vector_of_norm_V1_at_half(Mat2D V1, Mat2D Q, double gamma, double Pr_inf, double T_inf, double norm_delta_x, int N);
double calc_delta_Q_implicit_steger_warming(Mat2D delta_Q, Mat2D Q, Mat2D work_3_Np2_mat1, Mat2D work_3_Np1_mat1, Mat2D work_3_Np1_mat2, Mat2D work_3_3_mat1, Mat2D work_3_3_mat2, Mat2D work_3_3_mat3, Mat2D work_3_3_mat4, Mat2D work_3_3_mat5, Mat2D work_3_1_mat1, Mat2D work_3_1_mat2, Mat2D work_3_1_mat3, double *work_Np2_3_3_array1, double *work_Np2_3_3_array2, double *work_Np2_3_3_array3, double *work_3_Np2_array1, double gamma, double epsilon, double M_inf, double Re_inf, double Pr_inf, double c_v, double T_inf, double norm_delta_t, double  norm_delta_x, int N);
int btri3s(double *a, double *b, double *c, double *f, int kd, int ks, int ke);
double calc_delta_Q(Mat2D delta_Q, Mat2D Q, Mat2D work_3_Np2_mat1, Mat2D work_3_Np1_mat1, Mat2D work_3_Np1_mat2, Mat2D work_3_3_mat1, Mat2D work_3_3_mat2, Mat2D work_3_3_mat3, Mat2D work_3_3_mat4, Mat2D work_3_3_mat5, Mat2D work_3_1_mat1, Mat2D work_3_1_mat2, Mat2D work_3_1_mat3, double *work_Np2_3_3_array1, double *work_Np2_3_3_array2, double *work_Np2_3_3_array3, double *work_3_Np2_array1, double gamma, double epsilon, double M_inf, double Re_inf, double Pr_inf, double c_v, double T_inf, double norm_delta_t, double  norm_delta_x, int N, t_flag flags);

int main(int argc, char const *argv[])
{
/* declarations */
    char input_file[MAXDIR], init_conditions_file[MAXDIR], output_dir[MAXDIR], temp_word[MAXWORD];
    int N, max_iterations;
    double x_min, x_max, L, norm_delta_x, Re_inf, M_inf, CFL, gamma, Pr_inf, R_specific, c_v, T_inf, const_norm_delta_time, final_time, epsilon = 1e-2, norm_delta_t, current_norma, elapsed_norm_time = 0, *work_Np2_3_3_array1, *work_Np2_3_3_array2, *work_Np2_3_3_array3, *work_3_Np2_array1;
    t_flag flags = 0;
    Mat2D init_Q, current_Q, next_Q, delta_Q, work_3_Np2_mat1, work_3_Np2_mat2, work_3_Np1_mat1, work_3_Np1_mat2, work_3_3_mat1, work_3_3_mat2, work_3_3_mat3, work_3_3_mat4, work_3_3_mat5, work_3_1_mat1, work_3_1_mat2, work_3_1_mat3;
    FILE *output_Q_file, *output_iter_data_file;

/* getting the input file and output file */
    if (--argc != 3) {
        fprintf(stderr, "%s:%d: [Error] not right usage... Usage: main 'input file' 'initial conditions file' 'output directory'\n", __FILE__, __LINE__);
        return -1;
    }

    strncpy(input_file, (*(++argv)), MAXDIR);

    if (input_file[MAXDIR-1] != '\0') {
        fprintf(stderr, "%s:%d: [Error] input too long\n", __FILE__, __LINE__);
        return -1;
    }

    strncpy(init_conditions_file, (*(++argv)), MAXDIR);

    if (init_conditions_file[MAXDIR-1] != '\0') {
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
    read_input(input_file, &N, &x_min, &x_max, &L, &norm_delta_x, &Re_inf, &M_inf, &CFL, &gamma, &Pr_inf, &R_specific, &c_v, &T_inf, &max_iterations, &const_norm_delta_time, &final_time, &flags);

/* Checking the input */
    dprintINT(N);
    dprintDOUBLE(x_min);
    dprintDOUBLE(x_max);
    dprintDOUBLE(L);
    dprintDOUBLE(norm_delta_x);
    dprintDOUBLE(Re_inf);
    dprintDOUBLE(M_inf);
    dprintDOUBLE(CFL);
    dprintDOUBLE(gamma);
    dprintDOUBLE(Pr_inf);
    dprintDOUBLE(R_specific);
    dprintDOUBLE(c_v);
    dprintDOUBLE(T_inf);
    dprintDOUBLE(final_time);
    dprintDOUBLE(const_norm_delta_time);
    dprintINT(max_iterations);
    dprintSCHEME(flags);
    printf("--------------------\n");

/*------------------------------------------------------------*/
/* creating output directory */
    if (ON_LINUX) {
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
            return -1;
        }
        sprintf(temp_word, "/N%d", N);
        strcat(output_dir, temp_word);
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
            return -1;
        }
        if (flags & EXPLICIT_SW) {
            strcat(output_dir, "/EXPLICIT_SW");
            if (create_empty_dir(output_dir) != 0) {
                fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
                return -1;
            }
        } else if (flags & IMPLICIT_SW) {
            strcat(output_dir, "/IMPLICIT_SW");
            if (create_empty_dir(output_dir) != 0) {
                fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
                return -1;
            }
        }
        sprintf(temp_word, "/M%g", M_inf);
        strcat(output_dir, temp_word);
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
            return -1;
        }
        sprintf(temp_word, "/CFL%g", CFL);
        strcat(output_dir, temp_word);
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
            return -1;
        }
    }

    strcpy(temp_word, output_dir);
    strcat(temp_word, "/output_Q.txt");
    output_Q_file = fopen(temp_word, "wt");

    strcpy(temp_word, output_dir);
    strcat(temp_word, "/output_iter_data.txt");
    output_iter_data_file = fopen(temp_word, "wt");
    fprintf(output_iter_data_file, "%s, %s, %s, %s\n", "No", "norma", "norm_delta_time", "elapsed_norm_time");

/*------------------------------------------------------------*/
/* allocating the matrices */
    mat2D_alloc(&init_Q         , 3, N+2);
    mat2D_alloc(&current_Q      , 3, N+2);
    mat2D_alloc(&next_Q         , 3, N+2);
    mat2D_alloc(&delta_Q        , 3, N+2);
    mat2D_alloc(&work_3_Np2_mat1, 3, N+2);
    mat2D_alloc(&work_3_Np2_mat2, 3, N+2);

    mat2D_alloc(&work_3_Np1_mat1, 3, N+1);
    mat2D_alloc(&work_3_Np1_mat2, 3, N+1);

    mat2D_alloc(&work_3_3_mat1, 3, 3);
    mat2D_alloc(&work_3_3_mat2, 3, 3);
    mat2D_alloc(&work_3_3_mat3, 3, 3);
    mat2D_alloc(&work_3_3_mat4, 3, 3);
    mat2D_alloc(&work_3_3_mat5, 3, 3);

    mat2D_alloc(&work_3_1_mat1, 3, 1);
    mat2D_alloc(&work_3_1_mat2, 3, 1);
    mat2D_alloc(&work_3_1_mat3, 3, 1);

    work_Np2_3_3_array1 = (double *)calloc((N + 2) * 3 * 3, sizeof(double));
    work_Np2_3_3_array2 = (double *)calloc((N + 2) * 3 * 3, sizeof(double));
    work_Np2_3_3_array3 = (double *)calloc((N + 2) * 3 * 3, sizeof(double));
    work_3_Np2_array1   = (double *)calloc((N + 2) * 3, sizeof(double));

/*------------------------------------------------------------*/
/* initialization */
    current_norma = initialize_Q(init_conditions_file, init_Q, gamma, N);
    apply_BC(init_Q, N);

    mat2D_copy(current_Q, init_Q);
    mat2D_copy(next_Q, init_Q);

    norm_delta_t = calc_norm_delta_t(current_Q, gamma, R_specific, T_inf, norm_delta_x, const_norm_delta_time, CFL, N, flags);

    print_mat2D_to_file(output_Q_file, current_Q);
    fprintf(output_iter_data_file, "%d, %g, %g, %g\n", 0, current_norma, norm_delta_t, elapsed_norm_time);

/*------------------------------------------------------------*/
/* the loop */
    for (int i = 0; i < 10000; i++) {
        apply_BC(current_Q, N);

        norm_delta_t = calc_norm_delta_t(current_Q, gamma, R_specific, T_inf, norm_delta_x, const_norm_delta_time, CFL, N, flags);

        current_norma = calc_delta_Q(delta_Q, current_Q, work_3_Np2_mat1, work_3_Np1_mat1, work_3_Np1_mat2, work_3_3_mat1, work_3_3_mat2, work_3_3_mat3, work_3_3_mat4, work_3_3_mat5, work_3_1_mat1, work_3_1_mat2, work_3_1_mat3, work_Np2_3_3_array1, work_Np2_3_3_array2, work_Np2_3_3_array3, work_3_Np2_array1, gamma, epsilon, M_inf, Re_inf, Pr_inf, c_v, T_inf, norm_delta_t, norm_delta_x, N, flags);

        mat2D_add(next_Q, delta_Q);
        mat2D_copy(current_Q, next_Q);

        elapsed_norm_time += norm_delta_t;

        if (i % 100 == 0) {
            printf("%d\n",i);
        }
        if (0 == i) {
            fprintf(output_iter_data_file, "%d, %g, %g, %g\n", i+1, current_norma, norm_delta_t, elapsed_norm_time);
            print_mat2D_to_file(output_Q_file, current_Q);
        } else if (i % 1 == 0) {
            fprintf(output_iter_data_file, "%d, %g, %g, %g\n", i+1, current_norma, norm_delta_t, elapsed_norm_time);
            print_mat2D_to_file(output_Q_file, current_Q);
        }

        if (1.5 < elapsed_norm_time) {
            break;
        }
    }

/*------------------------------------------------------------*/
/* output */
    output_metadata(output_dir, N, x_min, x_max, Re_inf, M_inf, CFL, gamma, Pr_inf, R_specific, T_inf, max_iterations, final_time, flags);

/*------------------------------------------------------------*/
/* freeing the memory */
    mat2D_free(init_Q         );
    mat2D_free(current_Q      );
    mat2D_free(next_Q         );
    mat2D_free(delta_Q        );
    mat2D_free(work_3_Np2_mat1);
    mat2D_free(work_3_Np2_mat2);

    mat2D_free(work_3_Np1_mat1);
    mat2D_free(work_3_Np1_mat2);

    mat2D_free(work_3_3_mat1);
    mat2D_free(work_3_3_mat2);
    mat2D_free(work_3_3_mat3);
    mat2D_free(work_3_3_mat4);
    mat2D_free(work_3_3_mat5);

    mat2D_free(work_3_1_mat1);
    mat2D_free(work_3_1_mat2);
    mat2D_free(work_3_1_mat3);

    free(work_Np2_3_3_array1);
    free(work_Np2_3_3_array2);
    free(work_Np2_3_3_array3);

    free(work_3_Np2_array1);

    fclose(output_Q_file);
    fclose(output_iter_data_file);
    
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
void read_input(char *input_file, int *N, double *x_min, double *x_max, double *L, double *norm_delta_x, double *Re_inf, double *M_inf, double *CFL, double *gamma, double *Pr_inf, double *R_specific, double *c_v, double *T_inf, int *max_iterations, double *const_norm_delta_time, double *final_time, t_flag *flags)
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
        } else if (!strcmp(current_word, "Pr_inf")) {
            fscanf(fp, "%g", &temp_f);
            *Pr_inf = (double)temp_f;
        } else if (!strcmp(current_word, "R_specific")) {
            fscanf(fp, "%g", &temp_f);
            *R_specific = (double)temp_f;
        } else if (!strcmp(current_word, "c_v")) {
            fscanf(fp, "%g", &temp_f);
            *c_v = (double)temp_f;
        } else if (!strcmp(current_word, "T_inf")) {
            fscanf(fp, "%g", &temp_f);
            *T_inf = (double)temp_f;
        } else if (!strcmp(current_word, "iterations")) {
            fscanf(fp, "%g", &temp_f);
            *max_iterations = (int)temp_f;
        } else if (!strcmp(current_word, "norm_delta_time")) {
            fscanf(fp, "%g", &temp_f);
            *const_norm_delta_time = (double)temp_f;
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
    *norm_delta_x = (*L) / (*N) / (*L);

    if (*const_norm_delta_time) {
        *flags |= CONST_DELTA_TIME;
    }

    fclose(fp);
}

void output_metadata(char *output_dir, int N, double x_min, double x_max, double Re_inf, double M_inf, double CFL, double gamma, double Pr_inf, double R_specific, double T_inf, int max_iterations, double final_time, t_flag flags)
{
    char temp_word[MAXWORD], scheme[MAXWORD];
    FILE *metadata_file;

    strcpy(temp_word, output_dir);
    strcat(temp_word, "/metadata.txt");
    metadata_file = fopen(temp_word, "wt");

    if (flags & EXPLICIT_SW) {
        strcpy(scheme, "EXPLICIT SW");
    } else if(flags & IMPLICIT_SW) {
        strcpy(scheme, "IMPLICIT SW");
    } else if(flags & EXPLICIT_ROE) {
        strcpy(scheme, "EXPLICIT ROE");
    }

    fprintf(metadata_file, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n" , "N", "x_min", "x_max", "Re_inf", "M_inf", "CFL", "scheme", "gamma", "Pr_inf", "R_specific", "T_inf", "iterations", "final_time");
    fprintf(metadata_file, "%d, %g, %g, %g, %g, %g, %s, %g, %g, %g, %g, %d, %g" , N, x_min, x_max, Re_inf, M_inf, CFL, scheme, gamma, Pr_inf, R_specific, T_inf, max_iterations, final_time);


    fclose(metadata_file);
}

void print_mat2D_to_file(FILE *fp, Mat2D m)
{
    fprintf(fp, "\n");
    for (size_t i = 0; i < m.rows; i++) {
        for (size_t j = 0; j < m.cols; j++) {
            fprintf(fp, "%g ", MAT2D_AT(m, i, j));
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

/* converts a 2D index into 1D index
argument list:
i    - first direction
j    - second direction
cols - first direction size
rows - second direction size */
int offset2d(int i, int j, int cols, int rows)
{
    assert(i < cols);
    assert(j < rows);

    return j * cols + i;
}

/* converts a 3D index into 1D index
argument list:
i     - first direction
j     - second direction
k     - third direction
cols  - first direction size
rows  - second direction size
layer - third direction size */
int offset3d(int i, int j, int k, int cols, int rows, int layers)
{
    assert(i < cols);
    assert(j < rows);
    assert(k < layers);

    return (k * rows + j) * cols + i;
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

double calc_norm_energy(double gamma, double norm_rho, double norm_u, double norm_p)
{
    return (norm_p) / (gamma - 1) + 0.5 * norm_rho * norm_u * norm_u;
}

double calc_norm_p(double gamma, double norm_rho, double norm_u, double norm_e)
{
    return (gamma - 1) * (norm_e - 0.5 * norm_rho * norm_u * norm_u);
}

double calc_norm_total_enthalpy_at_i(Mat2D Q, double gamma, int i)
{
    double norm_rho = MAT2D_AT(Q, 0, i);
    double norm_u = MAT2D_AT(Q, 1, i) / norm_rho;
    double norm_e = MAT2D_AT(Q, 2, i);
    double norm_p = calc_norm_p(gamma, norm_rho, norm_u, norm_e);

    return (norm_e + norm_p) / norm_rho;
}

double calc_norm_T(Mat2D Q, double gamma, int i)
{
    double norm_rho = MAT2D_AT(Q, 0, i);
    double norm_u = MAT2D_AT(Q, 1, i) / norm_rho;
    double norm_e = MAT2D_AT(Q, 2, i);
    double norm_p = calc_norm_p(gamma, norm_rho, norm_u, norm_e);

    return norm_p / norm_rho;
}

double calc_norm_delta_t(Mat2D Q, double gamma, double R_specific, double T_inf, double norm_delta_x, double const_norm_delta_time, double CFL, int N, t_flag flags)
{
    double delta_t_CFL, delta_t_r, norm_rho, norm_u, norm_e, norm_p, norm_T, norm_a, max_allowed_norm_delta_time; 
    double maximum, norm_u_max = 0, norm_T_max = 0;

    for (int i = 1; i <= N; i++) {
        norm_rho = MAT2D_AT(Q, 0, i);
        norm_u = MAT2D_AT(Q, 1, i) / norm_rho;
        norm_e = MAT2D_AT(Q, 2, i);
        norm_p = calc_norm_p(gamma, norm_rho, norm_u, norm_e);
        norm_T = norm_p / (norm_rho * R_specific);
        norm_a = sqrt((gamma * norm_p) / norm_rho);

        maximum = fmax(fabs(norm_u), fmax(fabs(norm_u + norm_a), fabs(norm_u - norm_a)));

        if (maximum > norm_u_max) {
            norm_u_max = maximum;
        }
        if (norm_T > norm_T_max) {
            norm_T_max = norm_T;
        }
    }

    delta_t_CFL = (CFL * norm_delta_x) / norm_u_max;
    delta_t_r = (norm_delta_x * norm_delta_x) / (2 * calc_norm_mu(norm_T_max, gamma, T_inf));

    max_allowed_norm_delta_time = fmin(delta_t_CFL, delta_t_r);

    if (flags & CONST_DELTA_TIME) {
        if ((const_norm_delta_time > max_allowed_norm_delta_time)) {
            if (flags & IMPLICIT_SW) {
                if ((const_norm_delta_time > delta_t_r)) {
                    printf("[Warning] dt might be too big\n");
                    printf("          const dt         -> %g\n", const_norm_delta_time);
                    printf("          max allowed dt r -> %g\n", delta_t_r);
                }
            } else {
                printf("[Warning] dt might be too big\n");
                printf("          const dt           -> %g\n", const_norm_delta_time);
                printf("          max allowed dt r   -> %g\n", delta_t_r);
                printf("          max allowed dt CFL -> %g\n", delta_t_CFL);
            }
        }
        return const_norm_delta_time;
    } else {
        return max_allowed_norm_delta_time;
    }
}

void calc_T_matrix_at_i(Mat2D T_matrix, Mat2D Q, double gamma, int i)
{
    double norm_rho = MAT2D_AT(Q, 0, i);
    double norm_u = MAT2D_AT(Q, 1, i) / norm_rho;
    double norm_e = MAT2D_AT(Q, 2, i);
    double norm_p = calc_norm_p(gamma, norm_rho, norm_u, norm_e);
    double norm_a = sqrt(gamma * norm_p / norm_rho);

    MAT2D_AT(T_matrix, 0, 0) = 1;
    MAT2D_AT(T_matrix, 0, 1) = norm_rho / 2 / norm_a;
    MAT2D_AT(T_matrix, 0, 2) = - norm_rho / 2 / norm_a;
    MAT2D_AT(T_matrix, 1, 0) = norm_u;
    MAT2D_AT(T_matrix, 1, 1) = norm_rho / 2 / norm_a * (norm_u + norm_a);
    MAT2D_AT(T_matrix, 1, 2) = - norm_rho / 2 / norm_a * (norm_u - norm_a);
    MAT2D_AT(T_matrix, 2, 0) = norm_u * norm_u / 2;
    MAT2D_AT(T_matrix, 2, 1) = norm_rho / 2 / norm_a * (norm_u * norm_u / 2 + norm_u * norm_a + norm_a * norm_a / (gamma - 1));
    MAT2D_AT(T_matrix, 2, 2) = - norm_rho / 2 / norm_a * (norm_u * norm_u / 2 - norm_u * norm_a + norm_a * norm_a / (gamma - 1));
}

void calc_T_inverse_matrix_at_i(Mat2D T_inverse_matrix, Mat2D Q, double gamma, int i)
{
    double norm_rho = MAT2D_AT(Q, 0, i);
    double norm_u = MAT2D_AT(Q, 1, i) / norm_rho;
    double norm_e = MAT2D_AT(Q, 2, i);
    double norm_p = calc_norm_p(gamma, norm_rho, norm_u, norm_e);
    double norm_a = sqrt(gamma * norm_p / norm_rho);

    MAT2D_AT(T_inverse_matrix, 0, 0) = 1 - (gamma - 1) / 2 * norm_u * norm_u / norm_a / norm_a;
    MAT2D_AT(T_inverse_matrix, 0, 1) = (gamma - 1) * norm_u * norm_u / norm_a / norm_a;
    MAT2D_AT(T_inverse_matrix, 0, 2) = - (gamma - 1) / norm_a / norm_a;
    MAT2D_AT(T_inverse_matrix, 1, 0) = 1 / norm_rho / norm_a * ((gamma - 1) * norm_u * norm_u - norm_u * norm_a);
    MAT2D_AT(T_inverse_matrix, 1, 1) = 1 / norm_rho / norm_a * (norm_a - (gamma - 1) * norm_u);
    MAT2D_AT(T_inverse_matrix, 1, 2) = (gamma - 1) / norm_rho / norm_a;
    MAT2D_AT(T_inverse_matrix, 2, 0) = - 1 / norm_rho / norm_a * ((gamma - 1) * norm_u * norm_u + norm_u * norm_a);
    MAT2D_AT(T_inverse_matrix, 2, 1) = 1 / norm_rho / norm_a * (norm_a + (gamma - 1) * norm_u);
    MAT2D_AT(T_inverse_matrix, 2, 2) = - (gamma - 1) / norm_rho / norm_a;
}

void calc_lambda_plus_matrix_at_i(Mat2D lambda_plus_matrix, Mat2D Q, double gamma, int i, double epsilon)
{
    double norm_rho = MAT2D_AT(Q, 0, i);
    double norm_u = MAT2D_AT(Q, 1, i) / norm_rho;
    double norm_e = MAT2D_AT(Q, 2, i);
    double norm_p = calc_norm_p(gamma, norm_rho, norm_u, norm_e);
    double norm_a = sqrt(gamma * norm_p / norm_rho);
    
    mat2D_identity_mat(lambda_plus_matrix);
    mat2D_mult(lambda_plus_matrix, norm_u);
    MAT2D_AT(lambda_plus_matrix, 1, 1) += norm_a;
    MAT2D_AT(lambda_plus_matrix, 2, 2) -= norm_a;

    for (int i = 0; i < 3; i++) {
        MAT2D_AT(lambda_plus_matrix, i, i) = (MAT2D_AT(lambda_plus_matrix, i, i) + sqrt(MAT2D_AT(lambda_plus_matrix, i, i) * MAT2D_AT(lambda_plus_matrix, i, i) + epsilon * epsilon)) / 2;
    }
}

void calc_lambda_minus_matrix_at_i(Mat2D lambda_minus_matrix, Mat2D Q, double gamma, int i, double epsilon)
{
    double norm_rho = MAT2D_AT(Q, 0, i);
    double norm_u = MAT2D_AT(Q, 1, i) / norm_rho;
    double norm_e = MAT2D_AT(Q, 2, i);
    double norm_p = calc_norm_p(gamma, norm_rho, norm_u, norm_e);
    double norm_a = sqrt(gamma * norm_p / norm_rho);

    mat2D_identity_mat(lambda_minus_matrix);
    mat2D_mult(lambda_minus_matrix, norm_u);
    MAT2D_AT(lambda_minus_matrix, 1, 1) += norm_a;
    MAT2D_AT(lambda_minus_matrix, 2, 2) -= norm_a;

    for (int i = 0; i < 3; i++) {
        MAT2D_AT(lambda_minus_matrix, i, i) = (MAT2D_AT(lambda_minus_matrix, i, i) - sqrt(MAT2D_AT(lambda_minus_matrix, i, i) * MAT2D_AT(lambda_minus_matrix, i, i) + epsilon * epsilon)) / 2;
    }
}

void calc_A_plus_or_minus_at_i(Mat2D Q, Mat2D norm_A, Mat2D work_3_3_mat1, Mat2D work_3_3_mat2, Mat2D work_3_3_mat3, Mat2D work_3_3_mat4, double gamma, double epsilon, char sign, int i)
{
    assert('p' == sign || 'm' == sign);

    Mat2D norm_T, norm_inverse_T, norm_lambda, m_temp;

    norm_T         = work_3_3_mat1;
    norm_inverse_T = work_3_3_mat2;
    norm_lambda    = work_3_3_mat3;
    m_temp         = work_3_3_mat4;

    calc_T_matrix_at_i(norm_T, Q, gamma, i);
    calc_T_inverse_matrix_at_i(norm_inverse_T, Q, gamma, i);

    mat2D_fill(m_temp, 0);
    mat2D_fill(norm_A, 0);
    if ('p' == sign) {
        calc_lambda_plus_matrix_at_i(norm_lambda, Q, gamma, i, epsilon);
    } else if ('m' == sign) {
        calc_lambda_minus_matrix_at_i(norm_lambda, Q, gamma, i, epsilon);
    }
    mat2D_dot(m_temp, norm_lambda, norm_inverse_T);
    mat2D_dot(norm_A, norm_T, m_temp);
}

void calc_lambda_ave_pm_matrix_at_i(Mat2D lambda_pm_matrix, Mat2D Q, double gamma, char sign, int i, double epsilon)
{
    double norm_rho_i, norm_u_i, norm_tot_H_i, norm_rho_ip1, norm_u_ip1, norm_tot_H_ip1, ave_norm_rho, ave_norm_u, ave_norm_tot_H, ave_norm_a, R;

    norm_rho_i   = MAT2D_AT(Q, 0, i);
    norm_u_i     = MAT2D_AT(Q, 1, i) / norm_rho_i;
    norm_tot_H_i = calc_norm_total_enthalpy_at_i(Q, gamma, i);

    norm_rho_ip1   = MAT2D_AT(Q, 0, i+1);
    norm_u_ip1     = MAT2D_AT(Q, 1, i+1) / norm_rho_ip1;
    norm_tot_H_ip1 = calc_norm_total_enthalpy_at_i(Q, gamma, i+1);

    R = sqrt(norm_rho_ip1 / norm_rho_i);

    ave_norm_rho   = R * norm_rho_i;
    ave_norm_u     = (norm_u_i + R * norm_u_ip1) / (1 + R);
    ave_norm_tot_H = (norm_tot_H_i + R * norm_tot_H_ip1) / (1 + R);

    ave_norm_a = sqrt((gamma - 1) * (ave_norm_tot_H - 0.5 * ave_norm_u * ave_norm_u));
    
    mat2D_identity_mat(lambda_pm_matrix);
    mat2D_mult(lambda_pm_matrix, ave_norm_u);
    MAT2D_AT(lambda_pm_matrix, 1, 1) += ave_norm_a;
    MAT2D_AT(lambda_pm_matrix, 2, 2) -= ave_norm_a;

    for (int i = 0; i < 3; i++) {
        if ('p' == sign) {
            MAT2D_AT(lambda_pm_matrix, i, i) = (MAT2D_AT(lambda_pm_matrix, i, i) + sqrt(MAT2D_AT(lambda_pm_matrix, i, i) * MAT2D_AT(lambda_pm_matrix, i, i) + epsilon * epsilon)) / 2;
        } else if ('m' == sign) {
            MAT2D_AT(lambda_pm_matrix, i, i) = (MAT2D_AT(lambda_pm_matrix, i, i) - sqrt(MAT2D_AT(lambda_pm_matrix, i, i) * MAT2D_AT(lambda_pm_matrix, i, i) + epsilon * epsilon)) / 2;
        }
    }
}

void calc_T_matrix_roe_at_i(Mat2D T_matrix, Mat2D Q, double gamma, int i)
{
    double norm_rho_i, norm_u_i, norm_tot_H_i, norm_rho_ip1, norm_u_ip1, norm_tot_H_ip1, ave_norm_rho, ave_norm_u, ave_norm_tot_H, ave_norm_a, R;

    norm_rho_i   = MAT2D_AT(Q, 0, i);
    norm_u_i     = MAT2D_AT(Q, 1, i) / norm_rho_i;
    norm_tot_H_i = calc_norm_total_enthalpy_at_i(Q, gamma, i);

    norm_rho_ip1   = MAT2D_AT(Q, 0, i+1);
    norm_u_ip1     = MAT2D_AT(Q, 1, i+1) / norm_rho_ip1;
    norm_tot_H_ip1 = calc_norm_total_enthalpy_at_i(Q, gamma, i+1);

    R = sqrt(norm_rho_ip1 / norm_rho_i);

    ave_norm_rho   = R * norm_rho_i;
    ave_norm_u     = (norm_u_i + R * norm_u_ip1) / (1 + R);
    ave_norm_tot_H = (norm_tot_H_i + R * norm_tot_H_ip1) / (1 + R);

    ave_norm_a = sqrt((gamma - 1) * (ave_norm_tot_H - 0.5 * ave_norm_u * ave_norm_u));

    MAT2D_AT(T_matrix, 0, 0) = 1;
    MAT2D_AT(T_matrix, 0, 1) = ave_norm_rho / 2 / ave_norm_a;
    MAT2D_AT(T_matrix, 0, 2) = - ave_norm_rho / 2 / ave_norm_a;
    MAT2D_AT(T_matrix, 1, 0) = ave_norm_u;
    MAT2D_AT(T_matrix, 1, 1) = ave_norm_rho / 2 / ave_norm_a * (ave_norm_u + ave_norm_a);
    MAT2D_AT(T_matrix, 1, 2) = - ave_norm_rho / 2 / ave_norm_a * (ave_norm_u - ave_norm_a);
    MAT2D_AT(T_matrix, 2, 0) = ave_norm_u * ave_norm_u / 2;
    MAT2D_AT(T_matrix, 2, 1) = ave_norm_rho / 2 / ave_norm_a * (ave_norm_tot_H + ave_norm_u * ave_norm_a);
    MAT2D_AT(T_matrix, 2, 2) = - ave_norm_rho / 2 / ave_norm_a * (ave_norm_tot_H - ave_norm_u * ave_norm_a);
}

void calc_T_inverse_roe_matrix_at_i(Mat2D T_inverse_matrix, Mat2D Q, double gamma, int i)
{
    double norm_rho_i, norm_u_i, norm_tot_H_i, norm_rho_ip1, norm_u_ip1, norm_tot_H_ip1, ave_norm_rho, ave_norm_u, ave_norm_tot_H, ave_norm_a, R;

    norm_rho_i   = MAT2D_AT(Q, 0, i);
    norm_u_i     = MAT2D_AT(Q, 1, i) / norm_rho_i;
    norm_tot_H_i = calc_norm_total_enthalpy_at_i(Q, gamma, i);

    norm_rho_ip1   = MAT2D_AT(Q, 0, i+1);
    norm_u_ip1     = MAT2D_AT(Q, 1, i+1) / norm_rho_ip1;
    norm_tot_H_ip1 = calc_norm_total_enthalpy_at_i(Q, gamma, i+1);

    R = sqrt(norm_rho_ip1 / norm_rho_i);

    ave_norm_rho   = R * norm_rho_i;
    ave_norm_u     = (norm_u_i + R * norm_u_ip1) / (1 + R);
    ave_norm_tot_H = (norm_tot_H_i + R * norm_tot_H_ip1) / (1 + R);

    ave_norm_a = sqrt((gamma - 1) * (ave_norm_tot_H - 0.5 * ave_norm_u * ave_norm_u));

    MAT2D_AT(T_inverse_matrix, 0, 0) = 1 - (gamma - 1) / 2 * ave_norm_u * ave_norm_u / ave_norm_a / ave_norm_a;
    MAT2D_AT(T_inverse_matrix, 0, 1) = (gamma - 1) * ave_norm_u * ave_norm_u / ave_norm_a / ave_norm_a;
    MAT2D_AT(T_inverse_matrix, 0, 2) = - (gamma - 1) / ave_norm_a / ave_norm_a;
    MAT2D_AT(T_inverse_matrix, 1, 0) = 1 / ave_norm_rho / ave_norm_a * ((gamma - 1) * ave_norm_u * ave_norm_u - ave_norm_u * ave_norm_a);
    MAT2D_AT(T_inverse_matrix, 1, 1) = 1 / ave_norm_rho / ave_norm_a * (ave_norm_a - (gamma - 1) * ave_norm_u);
    MAT2D_AT(T_inverse_matrix, 1, 2) = (gamma - 1) / ave_norm_rho / ave_norm_a;
    MAT2D_AT(T_inverse_matrix, 2, 0) = - 1 / ave_norm_rho / ave_norm_a * ((gamma - 1) * ave_norm_u * ave_norm_u + ave_norm_u * ave_norm_a);
    MAT2D_AT(T_inverse_matrix, 2, 1) = 1 / ave_norm_rho / ave_norm_a * (ave_norm_a + (gamma - 1) * ave_norm_u);
    MAT2D_AT(T_inverse_matrix, 2, 2) = - (gamma - 1) / ave_norm_rho / ave_norm_a;
}

calc_A_roe_at_ip_half(Mat2D A_at_ip_half, Mat2D Q, double gamma, int i) 
{
    ;
}

void calc_norm_P_minus_Rx_matrix_at_i(Mat2D norm_P_minus_Rx, Mat2D Q, double norm_delta_x, double gamma, double T_inf, int i)
{
    double d_norm_mu_d_norm_x, norm_u_i, norm_T_i, norm_T_ipm1, norm_rho;

    norm_rho = MAT2D_AT(Q, 0, i);
    norm_u_i = MAT2D_AT(Q, 1, i)/norm_rho;
    norm_T_i = calc_norm_T(Q, gamma, i);

    if ((int)Q.cols-1 == i) {
        norm_T_ipm1 = calc_norm_T(Q, gamma, i-1);
        d_norm_mu_d_norm_x = (calc_norm_mu(norm_T_i, gamma, T_inf) - calc_norm_mu(norm_T_ipm1, gamma, T_inf)) / (2 * norm_delta_x);
    } else {
        norm_T_ipm1 = calc_norm_T(Q, gamma, i+1);
        d_norm_mu_d_norm_x = (calc_norm_mu(norm_T_ipm1, gamma, T_inf) - calc_norm_mu(norm_T_i, gamma, T_inf)) / (2 * norm_delta_x);
    }

    MAT2D_AT(norm_P_minus_Rx, 0, 0) = 0;
    MAT2D_AT(norm_P_minus_Rx, 0, 1) = 0;
    MAT2D_AT(norm_P_minus_Rx, 0, 2) = 0;
    MAT2D_AT(norm_P_minus_Rx, 1, 0) = norm_u_i * (double)4/3 * d_norm_mu_d_norm_x;
    MAT2D_AT(norm_P_minus_Rx, 1, 1) = - (double)4/3 * d_norm_mu_d_norm_x;
    MAT2D_AT(norm_P_minus_Rx, 1, 2) = 0;
    MAT2D_AT(norm_P_minus_Rx, 2, 0) = norm_u_i * norm_u_i * (double)4/3 * d_norm_mu_d_norm_x;
    MAT2D_AT(norm_P_minus_Rx, 2, 1) = - norm_u_i * (double)4/3 * d_norm_mu_d_norm_x;
    MAT2D_AT(norm_P_minus_Rx, 2, 2) = 0;

    mat2D_mult(norm_P_minus_Rx, (double)1/norm_rho);
}

void calc_norm_R_matrix_at_i(Mat2D norm_R, Mat2D Q, double gamma, double T_inf, double Pr_inf, double c_v, int i)
{
    double norm_rho_i, norm_u_i, norm_T_i, norm_mu_i, norm_kappa_i, norm_e_i, alpha;

    alpha = gamma / (Pr_inf * (gamma - 1));

    norm_rho_i   = MAT2D_AT(Q, 0, i);
    norm_u_i     = MAT2D_AT(Q, 1, i) / norm_rho_i;
    norm_e_i     = MAT2D_AT(Q, 2, i);
    norm_T_i     = calc_norm_T(Q, gamma, i);
    norm_mu_i    = calc_norm_mu(norm_T_i, gamma, T_inf);
    norm_kappa_i = calc_norm_kappa(norm_T_i, gamma, T_inf);

    MAT2D_AT(norm_R, 0, 0) = 0;
    MAT2D_AT(norm_R, 0, 1) = 0;
    MAT2D_AT(norm_R, 0, 2) = 0;
    MAT2D_AT(norm_R, 1, 0) = (double)4/3 * norm_u_i * norm_mu_i;
    MAT2D_AT(norm_R, 1, 1) = - (double)4/3 * norm_mu_i;
    MAT2D_AT(norm_R, 1, 2) = 0;
    MAT2D_AT(norm_R, 2, 0) = norm_u_i * norm_u_i * ((double)4/3 * norm_mu_i - alpha * norm_kappa_i / c_v) + alpha * norm_kappa_i / c_v * norm_e_i / norm_rho_i;
    MAT2D_AT(norm_R, 2, 1) = - norm_u_i * ((double)4/3 * norm_mu_i - alpha * norm_kappa_i / c_v);
    MAT2D_AT(norm_R, 2, 2) = - alpha * norm_kappa_i / c_v;

    mat2D_mult(norm_R, - (double)1/norm_rho_i);
}

double initialize_Q(char *init_conditions_file, Mat2D Q, double gamma, int N)
{
    char current_word[MAXWORD];
    float norm_rho, norm_u, norm_e, norm_p;
    double buffer, temp;
    int i;

    FILE *fp = fopen(init_conditions_file, "rt");
    if (fp == NULL) {
        fprintf(stderr, "%s:%d:[Error] problem opening input file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }
    
    buffer = (double)(N - 100) / 2;
    assert(!(buffer - (int)buffer));

    for (int i = 0; i < 5; i++) {
        fscanf(fp, "%s", current_word);
    }

    for (int i = 1; i <= N; i++) {
        fscanf(fp, "%s", current_word);
        fscanf(fp, "%s", current_word);
        norm_rho = atof(current_word);
        fscanf(fp, "%s", current_word);
        norm_u = atof(current_word);
        fscanf(fp, "%s", current_word);
        norm_p = atof(current_word);

        norm_e = calc_norm_energy(gamma, norm_rho, norm_u, norm_p);

        MAT2D_AT(Q, 0, i) = norm_rho;
        MAT2D_AT(Q, 1, i) = norm_rho * norm_u;
        MAT2D_AT(Q, 2, i) = norm_e;

        if (1 == i) {
            temp = buffer;
            while (1) {
                i++;
                MAT2D_AT(Q, 0, i) = norm_rho;
                MAT2D_AT(Q, 1, i) = norm_rho * norm_u;
                MAT2D_AT(Q, 2, i) = norm_e;
                if (!temp) {
                    break;
                }
                temp--;
            }
        }
    }
    temp = buffer;
    i = buffer + 100;
    while (1) {
        i++;
        MAT2D_AT(Q, 0, i) = MAT2D_AT(Q, 0, i-1);
        MAT2D_AT(Q, 1, i) = MAT2D_AT(Q, 1, i-1);
        MAT2D_AT(Q, 2, i) = MAT2D_AT(Q, 2, i-1);
        if (!temp) {
            break;
        }
        temp--;
    }

    fclose(fp);

    return mat2D_calc_norma(Q);
}

void apply_BC(Mat2D Q, int N)
{
    MAT2D_AT(Q, 0, 0) =   MAT2D_AT(Q, 0, 1);
    MAT2D_AT(Q, 1, 0) = - MAT2D_AT(Q, 1, 1);
    MAT2D_AT(Q, 2, 0) =   MAT2D_AT(Q, 2, 1);

    MAT2D_AT(Q, 0, N+1) =   MAT2D_AT(Q, 0, N);
    MAT2D_AT(Q, 1, N+1) = - MAT2D_AT(Q, 1, N);
    MAT2D_AT(Q, 2, N+1) =   MAT2D_AT(Q, 2, N);
}

void calc_vector_of_E(Mat2D E, Mat2D Q, double gamma, int N)
{
    double norm_rho, norm_u, norm_e, norm_p;

    for (int i = 0; i < N+2; i++) {
        // dprintINT(i);
        norm_rho = MAT2D_AT(Q, 0, i);
        norm_u   = MAT2D_AT(Q, 1, i) / norm_rho;
        norm_e   = MAT2D_AT(Q, 2, i);
        norm_p   = calc_norm_p(gamma, norm_rho, norm_u, norm_e);

        MAT2D_AT(E, 0, i) = norm_rho * norm_u;
        MAT2D_AT(E, 1, i) = norm_p + norm_rho * norm_u * norm_u;
        MAT2D_AT(E, 2, i) = (norm_e + norm_p) * norm_u;
    }
}

void calc_vector_of_tilde_norm_E_at_half_sw(Mat2D tilde_norm_E, Mat2D Q, Mat2D work_3_3_mat1, Mat2D work_3_3_mat2, Mat2D work_3_3_mat3, Mat2D work_3_3_mat4, Mat2D work_3_3_mat5, Mat2D work_3_1_mat1, Mat2D work_3_1_mat2, Mat2D work_3_1_mat3, double gamma, double epsilon, int N)
{
    Mat2D norm_A, Q_index, temp_3_by_1, E_at_index;

    norm_A         = work_3_3_mat5;
    Q_index        = work_3_1_mat1;
    temp_3_by_1    = work_3_1_mat2;
    E_at_index     = work_3_1_mat3;

    for (int i = 0; i <= N; i++) {
        mat2D_fill(E_at_index, 0);

        /* norm_A_plus */
        calc_A_plus_or_minus_at_i(Q, norm_A, work_3_3_mat1, work_3_3_mat2, work_3_3_mat3, work_3_3_mat4, gamma, epsilon, 'p', i);

        mat2D_fill(temp_3_by_1, 0);
        MAT2D_AT(Q_index, 0, 0) = MAT2D_AT(Q, 0, i);
        MAT2D_AT(Q_index, 1, 0) = MAT2D_AT(Q, 1, i);
        MAT2D_AT(Q_index, 2, 0) = MAT2D_AT(Q, 2, i);
        mat2D_dot(temp_3_by_1, norm_A, Q_index);
        mat2D_add(E_at_index, temp_3_by_1);

        /* norm_A_minus */
        calc_A_plus_or_minus_at_i(Q, norm_A, work_3_3_mat1, work_3_3_mat2, work_3_3_mat3, work_3_3_mat4, gamma, epsilon, 'm', i+1);

        mat2D_fill(temp_3_by_1, 0);
        MAT2D_AT(Q_index, 0, 0) = MAT2D_AT(Q, 0, i+1);
        MAT2D_AT(Q_index, 1, 0) = MAT2D_AT(Q, 1, i+1);
        MAT2D_AT(Q_index, 2, 0) = MAT2D_AT(Q, 2, i+1);
        mat2D_dot(temp_3_by_1, norm_A, Q_index);
        mat2D_add(E_at_index, temp_3_by_1);

        /* copying to vector of tilde norm E */
        MAT2D_AT(tilde_norm_E, 0, i) = MAT2D_AT(E_at_index, 0, 0);
        MAT2D_AT(tilde_norm_E, 1, i) = MAT2D_AT(E_at_index, 1, 0);
        MAT2D_AT(tilde_norm_E, 2, i) = MAT2D_AT(E_at_index, 2, 0);
    }
}

void calc_vector_of_norm_V1_at_half(Mat2D V1, Mat2D Q, double gamma, double Pr_inf, double T_inf, double norm_delta_x, int N)
{
    double norm_mu_at_half, norm_du_dx_at_half, norm_u_at_half, norm_kappa_at_half, norm_dT_dx_at_half, norm_T_i, norm_T_ip1, norm_u_i, norm_u_ip1;

    for (int i = 0; i <= N; i++) {
        norm_T_i = calc_norm_T(Q, gamma, i);
        norm_T_ip1 = calc_norm_T(Q, gamma, i+1);
        norm_u_i = MAT2D_AT(Q, 1, i) / MAT2D_AT(Q, 0, i);
        norm_u_ip1 = MAT2D_AT(Q, 1, i+1) / MAT2D_AT(Q, 0, i+1);
        
        norm_mu_at_half = (calc_norm_mu(norm_T_ip1, gamma, T_inf) + calc_norm_mu(norm_T_i, gamma, T_inf)) / 2;
        norm_kappa_at_half = (calc_norm_kappa(norm_T_ip1, gamma, T_inf) + calc_norm_kappa(norm_T_i, gamma, T_inf)) / 2;
        norm_u_at_half = (norm_u_ip1 - norm_u_i) / 2;
        norm_du_dx_at_half = (norm_u_ip1 - norm_u_i) / norm_delta_x;
        norm_dT_dx_at_half = (norm_T_ip1 - norm_T_i) / norm_delta_x;

        MAT2D_AT(V1, 0, i) = 0;
        MAT2D_AT(V1, 1, i) = (double)4/3 * norm_mu_at_half * norm_du_dx_at_half;
        MAT2D_AT(V1, 2, i) = (double)4/3 * norm_mu_at_half * norm_u_at_half * norm_du_dx_at_half + gamma / (Pr_inf * (gamma - 1)) * norm_kappa_at_half * norm_dT_dx_at_half;
    }
} 

void calc_vector_of_tilde_norm_E_at_half_roe()
{
    ;
}

double calc_delta_Q_explicit_steger_warming(Mat2D delta_Q, Mat2D Q, Mat2D work_3_Np1_mat1, Mat2D work_3_Np1_mat2, Mat2D work_3_3_mat1, Mat2D work_3_3_mat2, Mat2D work_3_3_mat3, Mat2D work_3_3_mat4, Mat2D work_3_3_mat5, Mat2D work_3_1_mat1, Mat2D work_3_1_mat2, Mat2D work_3_1_mat3, double gamma, double epsilon, double M_inf, double Re_inf, double Pr_inf, double T_inf, double norm_delta_t, double  norm_delta_x, int N)
{
    Mat2D tilde_norm_E, V1, norm_T, norm_inverse_T, norm_lambda, norm_A, m_temp, Q_index, temp_3_by_1, E_at_index, tilde_norm_E_at_half, V1_at_half, temp_3_1;

    tilde_norm_E   = work_3_Np1_mat1;
    V1             = work_3_Np1_mat2;
    norm_T         = work_3_3_mat1;
    norm_inverse_T = work_3_3_mat2;
    norm_lambda    = work_3_3_mat3;
    norm_A         = work_3_3_mat4;
    m_temp         = work_3_3_mat5;
    Q_index        = work_3_1_mat1;
    temp_3_by_1    = work_3_1_mat2;
    E_at_index     = work_3_1_mat3;

    calc_vector_of_tilde_norm_E_at_half_sw(tilde_norm_E, Q, norm_T, norm_inverse_T, norm_lambda, norm_A, m_temp, Q_index, temp_3_by_1, E_at_index, gamma, epsilon, N);
    calc_vector_of_norm_V1_at_half(V1, Q, gamma, Pr_inf, T_inf, norm_delta_x, N);

    tilde_norm_E_at_half = work_3_1_mat1;
    V1_at_half = work_3_1_mat2;
    temp_3_1 = work_3_1_mat3; 

    for (int i = 1; i <= N; i++) {
        mat2D_fill(temp_3_1, 0);

        mat2D_get_col(V1_at_half, V1, i);
        mat2D_add(temp_3_1, V1_at_half);

        mat2D_get_col(V1_at_half, V1, i-1);
        mat2D_sub(temp_3_1, V1_at_half);

        mat2D_mult(temp_3_1, - M_inf / Re_inf);

        mat2D_get_col(tilde_norm_E_at_half, tilde_norm_E, i);
        mat2D_add(temp_3_1, tilde_norm_E_at_half);

        mat2D_get_col(tilde_norm_E_at_half, tilde_norm_E, i-1);
        mat2D_sub(temp_3_1, tilde_norm_E_at_half);

        mat2D_mult(temp_3_1, - norm_delta_t / norm_delta_x);

        MAT2D_AT(delta_Q, 0, i) = MAT2D_AT(temp_3_1, 0, 0);
        MAT2D_AT(delta_Q, 1, i) = MAT2D_AT(temp_3_1, 1, 0);
        MAT2D_AT(delta_Q, 2, i) = MAT2D_AT(temp_3_1, 2, 0);
    }

    return mat2D_calc_norma(delta_Q);
}

double calc_delta_Q_implicit_steger_warming(Mat2D delta_Q, Mat2D Q, Mat2D work_3_Np2_mat1, Mat2D work_3_Np1_mat1, Mat2D work_3_Np1_mat2, Mat2D work_3_3_mat1, Mat2D work_3_3_mat2, Mat2D work_3_3_mat3, Mat2D work_3_3_mat4, Mat2D work_3_3_mat5, Mat2D work_3_1_mat1, Mat2D work_3_1_mat2, Mat2D work_3_1_mat3, double *work_Np2_3_3_array1, double *work_Np2_3_3_array2, double *work_Np2_3_3_array3, double *work_3_Np2_array1, double gamma, double epsilon, double M_inf, double Re_inf, double Pr_inf, double c_v, double T_inf, double norm_delta_t, double  norm_delta_x, int N)
{
    Mat2D temp_RHS, norm_A, norm_P_minus_Rx_matrix, norm_R_matrix, I;
    double *RHS, *theta, *phi, *psi;

    temp_RHS = work_3_Np2_mat1;
    RHS      = work_3_Np2_array1;
    theta    = work_Np2_3_3_array1;
    phi      = work_Np2_3_3_array2;
    psi      = work_Np2_3_3_array3;
    
    /* zeroing theta, phi, psi */
    for (int i = 0; i < N+2; i++) {
        for (int m = 0; m < 3; m++) {
            for (int  n = 0; n < 3; n++) {
                theta[offset3d(i, m, n, N+2, 3, 3)] = 0;
                phi[offset3d(i, m, n, N+2, 3, 3)] = 0;
                psi[offset3d(i, m, n, N+2, 3, 3)] = 0;
            }
        }
    }

    calc_delta_Q_explicit_steger_warming(temp_RHS, Q, work_3_Np1_mat1, work_3_Np1_mat2, work_3_3_mat1, work_3_3_mat2, work_3_3_mat3, work_3_3_mat4, work_3_3_mat5, work_3_1_mat1, work_3_1_mat2, work_3_1_mat3, gamma, epsilon, M_inf, Re_inf, Pr_inf, T_inf, norm_delta_t, norm_delta_x, N);

    /* copying temp_RHS to RHS according to right indexes */
    for (int i = 0; i < N+2; i++) {
        for (int j = 0; j < 3; j++) {
            RHS[offset2d(i, j, N+2, 3)] = MAT2D_AT(temp_RHS, j, i);    /* my i index is the row and the btri3s i index in the column */
        }
    }

    for (int i = 1; i <= N; i++) {
        /* Theta */
        norm_A = work_3_3_mat5;
        mat2D_fill(norm_A, 0);
        calc_A_plus_or_minus_at_i(Q, norm_A, work_3_3_mat1, work_3_3_mat2, work_3_3_mat3, work_3_3_mat4, gamma, epsilon, 'p', i-1);
        for (size_t m = 0; m < norm_A.rows; m++) {
            for (size_t n = 0; n < norm_A.cols; n++) {
                theta[offset3d(i, m, n, N+2, 3, 3)] += - norm_delta_t / norm_delta_x * MAT2D_AT(norm_A, m, n);
            }
        }
        norm_P_minus_Rx_matrix = work_3_3_mat5;
        mat2D_fill(norm_P_minus_Rx_matrix, 0);
        calc_norm_P_minus_Rx_matrix_at_i(norm_P_minus_Rx_matrix, Q, norm_delta_x, gamma, T_inf, i-1);
        for (size_t m = 0; m < norm_P_minus_Rx_matrix.rows; m++) {
            for (size_t n = 0; n < norm_P_minus_Rx_matrix.cols; n++) {
                theta[offset3d(i, m, n, N+2, 3, 3)] += norm_delta_t / (2 * norm_delta_x) * M_inf / Re_inf * MAT2D_AT(norm_P_minus_Rx_matrix, m, n);
            }
        }
        norm_R_matrix = work_3_3_mat5;
        mat2D_fill(norm_R_matrix, 0);
        calc_norm_R_matrix_at_i(norm_R_matrix, Q, gamma, T_inf, Pr_inf, c_v, i-1);
        for (size_t m = 0; m < norm_R_matrix.rows; m++) {
            for (size_t n = 0; n < norm_R_matrix.cols; n++) {
                theta[offset3d(i, m, n, N+2, 3, 3)] += - norm_delta_t / (norm_delta_x * norm_delta_x) * M_inf / Re_inf * MAT2D_AT(norm_R_matrix, m, n);
            }
        }

        /* Phi */
        norm_A = work_3_3_mat5;
        mat2D_fill(norm_A, 0);
        calc_A_plus_or_minus_at_i(Q, norm_A, work_3_3_mat1, work_3_3_mat2, work_3_3_mat3, work_3_3_mat4, gamma, epsilon, 'p', i);
        for (size_t m = 0; m < norm_A.rows; m++) {
            for (size_t n = 0; n < norm_A.cols; n++) {
                phi[offset3d(i, m, n, N+2, 3, 3)] += norm_delta_t / norm_delta_x * MAT2D_AT(norm_A, m, n);
            }
        }
        norm_A = work_3_3_mat5;
        mat2D_fill(norm_A, 0);
        calc_A_plus_or_minus_at_i(Q, norm_A, work_3_3_mat1, work_3_3_mat2, work_3_3_mat3, work_3_3_mat4, gamma, epsilon, 'm', i);
        for (size_t m = 0; m < norm_A.rows; m++) {
            for (size_t n = 0; n < norm_A.cols; n++) {
                phi[offset3d(i, m, n, N+2, 3, 3)] += - norm_delta_t / norm_delta_x * MAT2D_AT(norm_A, m, n);
            }
        }
        I = work_3_3_mat5;
        mat2D_identity_mat(I);
        for (size_t m = 0; m < I.rows; m++) {
            for (size_t n = 0; n < I.cols; n++) {
                phi[offset3d(i, m, n, N+2, 3, 3)] += MAT2D_AT(I, m, n);
            }
        }
        norm_R_matrix = work_3_3_mat5;
        mat2D_fill(norm_R_matrix, 0);
        calc_norm_R_matrix_at_i(norm_R_matrix, Q, gamma, T_inf, Pr_inf, c_v, i);
        for (size_t m = 0; m < norm_R_matrix.rows; m++) {
            for (size_t n = 0; n < norm_R_matrix.cols; n++) {
                phi[offset3d(i, m, n, N+2, 3, 3)] += 2 * norm_delta_t / (norm_delta_x * norm_delta_x) * M_inf / Re_inf * MAT2D_AT(norm_R_matrix, m, n);
            }
        }

        /* Psi */
        norm_A = work_3_3_mat5;
        mat2D_fill(norm_A, 0);
        calc_A_plus_or_minus_at_i(Q, norm_A, work_3_3_mat1, work_3_3_mat2, work_3_3_mat3, work_3_3_mat4, gamma, epsilon, 'm', i+1);
        for (size_t m = 0; m < norm_A.rows; m++) {
            for (size_t n = 0; n < norm_A.cols; n++) {
                psi[offset3d(i, m, n, N+2, 3, 3)] += norm_delta_t / norm_delta_x * MAT2D_AT(norm_A, m, n);
            }
        }
        norm_P_minus_Rx_matrix = work_3_3_mat5;
        mat2D_fill(norm_P_minus_Rx_matrix, 0);
        calc_norm_P_minus_Rx_matrix_at_i(norm_P_minus_Rx_matrix, Q, norm_delta_x, gamma, T_inf, i+1);
        for (size_t m = 0; m < norm_P_minus_Rx_matrix.rows; m++) {
            for (size_t n = 0; n < norm_P_minus_Rx_matrix.cols; n++) {
                psi[offset3d(i, m, n, N+2, 3, 3)] += - norm_delta_t / (2 * norm_delta_x) * M_inf / Re_inf * MAT2D_AT(norm_P_minus_Rx_matrix, m, n);
            }
        }
        norm_R_matrix = work_3_3_mat5;
        mat2D_fill(norm_R_matrix, 0);
        calc_norm_R_matrix_at_i(norm_R_matrix, Q, gamma, T_inf, Pr_inf, c_v, i+1);
        for (size_t m = 0; m < norm_R_matrix.rows; m++) {
            for (size_t n = 0; n < norm_R_matrix.cols; n++) {
                psi[offset3d(i, m, n, N+2, 3, 3)] += - norm_delta_t / (norm_delta_x * norm_delta_x) * M_inf / Re_inf * MAT2D_AT(norm_R_matrix, m, n);
            }
        }
    }

    btri3s(theta, phi, psi, RHS, N+2, 1, N);

    /* copying the solution (RHS) to delta Q according to right indexes */
    for (int i = 0; i < N+2; i++) {
        for (int j = 0; j < 3; j++) {
            MAT2D_AT(delta_Q, j, i) = RHS[offset2d(i, j, N+2, 3)];    /* my i index is the row and the btri3s i index in the column */
        }
    }

    return mat2D_calc_norma(delta_Q);
}

double calc_delta_Q_explicit_roe()
{
    return 0;
}

/* 3x3 block tri-diagonal matrix solver 
argument list:
a  - 1D array of 3D matrix 
b  - 1D array of 3D matrix 
c  - 1D array of 3D matrix 
f  - 1D array of 2D matrix 
kd - N+2
ks - 1
ke - N */
int btri3s(double *a, double *b, double *c, double *f, int kd, int ks, int ke)
{
  /* Local variables */
  int k, m, n, nd, md;

  double c1, d1, d2, d3, c2, c3, b11, b21, b22, b31, b32, b33, u12, u13, u23;
  
  
  /*   (A,B,C)F = F, F and B are overloaded, solution in F */

  md = 3;
  nd = 3;

  /*   Part 1. Forward block sweep */
  
  for (k = ks; k <= ke; k++)
    {
      
      /*      Step 1. Construct L in B */
      
      if (k != ks) 
	{
	  for (m = 0; m < md; m++) 
	    {
	      for (n = 0; n < nd; n++) 
		{
		  b[k + kd * (m + md * n)] = b[k + kd * (m + md * n)] 
		    - a[k + kd * (m + md * 0)] * b[k - 1 + kd * (0 + md * n)] 
		    - a[k + kd * (m + md * 1)] * b[k - 1 + kd * (1 + md * n)] 
		    - a[k + kd * (m + md * 2)] * b[k - 1 + kd * (2 + md * n)] ;
		}
	    }
	}
      
      /*      Step 2. Compute L inverse (block matrix) */
      
      /*          A. Decompose L into L and U */
      
      b11 = 1. / b[k + kd * (0 + md * 0)];
      u12 = b[k + kd * (0 + md * 1)] * b11;
      u13 = b[k + kd * (0 + md * 2)] * b11;
      b21 = b[k + kd * (1 + md * 0)];
      b22 = 1. / (b[k + kd * (1 + md * 1)] - b21 * u12);
      u23 = (b[k + kd * (1 + md * 2)] - b21 * u13) * b22;
      b31 = b[k + kd * (2 + md * 0)];
      b32 = b[k + kd * (2 + md * 1)] - b31 * u12;
      b33 = 1. / (b[k + kd * (2 + md * 2)] - b31 * u13 - b32 * u23);
      
      /*      Step 3. Solve for intermediate vector */
      
      /*          A. Construct RHS */
      if (k != ks) 
	{
	  for (m = 0; m < md; m++) 
	    {
	      f[k + kd * m] = f[k + kd * m] 
		- a[k + kd * (m + md * 0)] * f[k - 1 + kd * 0] 
		- a[k + kd * (m + md * 1)] * f[k - 1 + kd * 1] 
		- a[k + kd * (m + md * 2)] * f[k - 1 + kd * 2] ;
	    }
	}
      
      /*          B. Intermediate vector */
      
      /*          Forward substitution */
      
      d1 = f[k + kd * 0] * b11;
      d2 = (f[k + kd * 1] - b21 * d1) * b22;
      d3 = (f[k + kd * 2] - b31 * d1 - b32 * d2) * b33;
      
      /*          Backward substitution */
      
      f[k + kd * 2] = d3;
      f[k + kd * 1] = d2 - u23 * f[k + kd * 2];
      f[k + kd * 0] = d1 - u12 * f[k + kd * 1] - u13 * f[k + kd * 2];
      
      /*      Step 4. Construct U = L ** (-1) * C */
      /*              by columns and store in B */
      
      if (k != ke) 
	{
	  for (n = 0; n < nd; n++) 
	    {

	      /*          Forward substitution */
	      
	      c1 = c[k + kd * (0 + md * n)] * b11;
	      c2 = (c[k + kd * (1 + md * n)] - b21 * c1) * b22;
	      c3 = (c[k + kd * (2 + md * n)] - b31 * c1 - b32 * c2) * b33;
	      
	      /*          Backward substitution */
	      
	      b[k + kd * (2 + md * n)] = c3;
	      b[k + kd * (1 + md * n)] = c2 - u23 * b[k + kd * (2 + md * n)];
	      b[k + kd * (0 + md * n)] = c1 - u12 * b[k + kd * (1 + md * n)] 
		- u13 * b[k + kd * (2 + md * n)];
	    }
	}
    }
  
  /*   Part 2. Backward block sweep */
  
  if (ke == ks) 
    {
      return 0;
    }

  for (k = ke - 1; k >= ks; --k) 
    {
      for (m = 0; m < md; m++) 
	{
	  f[k + kd * m] = f[k + kd * m] 
	    - b[k + kd * (m + md * 0)] * f[k + 1 + kd * 0] 
	    - b[k + kd * (m + md * 1)] * f[k + 1 + kd * 1] 
	    - b[k + kd * (m + md * 2)] * f[k + 1 + kd * 2] ;
	}
    }
  
  return 0;
  
}

double calc_delta_Q(Mat2D delta_Q, Mat2D Q, Mat2D work_3_Np2_mat1, Mat2D work_3_Np1_mat1, Mat2D work_3_Np1_mat2, Mat2D work_3_3_mat1, Mat2D work_3_3_mat2, Mat2D work_3_3_mat3, Mat2D work_3_3_mat4, Mat2D work_3_3_mat5, Mat2D work_3_1_mat1, Mat2D work_3_1_mat2, Mat2D work_3_1_mat3, double *work_Np2_3_3_array1, double *work_Np2_3_3_array2, double *work_Np2_3_3_array3, double *work_3_Np2_array1, double gamma, double epsilon, double M_inf, double Re_inf, double Pr_inf, double c_v, double T_inf, double norm_delta_t, double  norm_delta_x, int N, t_flag flags)
{
    double current_norma;

    if (flags & EXPLICIT_SW) {
        current_norma = calc_delta_Q_explicit_steger_warming(delta_Q, Q, work_3_Np1_mat1, work_3_Np1_mat2, work_3_3_mat1, work_3_3_mat2, work_3_3_mat3, work_3_3_mat4, work_3_3_mat5, work_3_1_mat1, work_3_1_mat2, work_3_1_mat3, gamma, epsilon, M_inf, Re_inf, Pr_inf, T_inf, norm_delta_t, norm_delta_x, N);
    } else if (flags & IMPLICIT_SW) {
        current_norma = calc_delta_Q_implicit_steger_warming(delta_Q, Q, work_3_Np2_mat1, work_3_Np1_mat1, work_3_Np1_mat2, work_3_3_mat1, work_3_3_mat2, work_3_3_mat3, work_3_3_mat4, work_3_3_mat5, work_3_1_mat1, work_3_1_mat2, work_3_1_mat3, work_Np2_3_3_array1, work_Np2_3_3_array2, work_Np2_3_3_array3, work_3_Np2_array1, gamma, epsilon, M_inf, Re_inf, Pr_inf, c_v, T_inf, norm_delta_t, norm_delta_x, N);
    } else if (flags & EXPLICIT_ROE) {
        ;
    }
    
    return current_norma;
}

