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
double calc_norm_energy(double gamma, double norm_rho, double norm_u, double norm_p);
double calc_norm_p(double gamma, double norm_rho, double norm_u, double norm_e);
void calc_T_matrix(Mat2D T_matrix, double gamma, double norm_rho, double norm_a, double norm_u);
void calc_T_inverse_matrix(Mat2D T_inverse_matrix, double gamma, double norm_rho, double norm_p, double norm_u);
void calc_lambda_plus_matrix(Mat2D lambda_plus_matrix, double gamma, double norm_rho, double norm_p, double norm_u, double epsilon);
void calc_lambda_minus_matrix(Mat2D lambda_minus_matrix, double gamma, double norm_rho, double norm_p, double norm_u, double epsilon);
void initialize_Q(char *init_conditions_file, Mat2D Q, double gamma, int N);
void calc_vector_of_E(Mat2D E, Mat2D Q, double gamma, int N);

int main(int argc, char const *argv[])
{
/* declarations */
    char input_file[MAXDIR], init_conditions_file[MAXDIR], output_dir[MAXDIR], temp_word[MAXWORD];
    int N, max_iterations;
    double x_min, x_max, L, delta_x, Re_inf, M_inf, CFL, gamma, T_inf, Pr_inf, final_time, epsilon = 1e-10;
    t_flag flags = 0;
    Mat2D init_Q, current_Q, next_Q, work_mat_3_by_N_2;

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
    mat2D_alloc(&init_Q         , 3, N+2);
    mat2D_alloc(&current_Q      , 3, N+2);
    mat2D_alloc(&next_Q         , 3, N+2);
    mat2D_alloc(&work_mat_3_by_N_2, 3, N+2);

/*------------------------------------------------------------*/
/* initialization */
    initialize_Q(init_conditions_file, init_Q, gamma, N);
    // MAT2D_PRINT(init_Q);

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

double calc_norm_energy(double gamma, double norm_rho, double norm_u, double norm_p)
{
    return (norm_p) / (gamma - 1) + 0.5 * norm_rho * norm_u * norm_u;
}

double calc_norm_p(double gamma, double norm_rho, double norm_u, double norm_e)
{
    return (gamma - 1) * (norm_e - 0.5 * norm_rho * norm_u * norm_u);
}

void calc_T_matrix(Mat2D T_matrix, double gamma, double norm_rho, double norm_p, double norm_u)
{
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

void calc_T_inverse_matrix(Mat2D T_inverse_matrix, double gamma, double norm_rho, double norm_p, double norm_u)
{
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

void calc_lambda_plus_matrix(Mat2D lambda_plus_matrix, double gamma, double norm_rho, double norm_p, double norm_u, double epsilon)
{
    double norm_a = sqrt(gamma * norm_p / norm_rho);

    mat2D_identity_mat(lambda_plus_matrix);
    mat2D_mult(lambda_plus_matrix, norm_u);
    MAT2D_AT(lambda_plus_matrix, 1, 1) += norm_a;
    MAT2D_AT(lambda_plus_matrix, 2, 2) -= norm_a;

    for (int i = 0; i < 3; i++) {
        MAT2D_AT(lambda_plus_matrix, i, i) = (MAT2D_AT(lambda_plus_matrix, i, i) + sqrt(MAT2D_AT(lambda_plus_matrix, i, i) * MAT2D_AT(lambda_plus_matrix, i, i) + epsilon * epsilon)) / 2;
    }
}

void calc_lambda_minus_matrix(Mat2D lambda_minus_matrix, double gamma, double norm_rho, double norm_p, double norm_u, double epsilon)
{
    double norm_a = sqrt(gamma * norm_p / norm_rho);

    mat2D_identity_mat(lambda_minus_matrix);
    mat2D_mult(lambda_minus_matrix, norm_u);
    MAT2D_AT(lambda_minus_matrix, 1, 1) += norm_a;
    MAT2D_AT(lambda_minus_matrix, 2, 2) -= norm_a;

    for (int i = 0; i < 3; i++) {
        MAT2D_AT(lambda_minus_matrix, i, i) = (MAT2D_AT(lambda_minus_matrix, i, i) - sqrt(MAT2D_AT(lambda_minus_matrix, i, i) * MAT2D_AT(lambda_minus_matrix, i, i) + epsilon * epsilon)) / 2;
    }
}

void initialize_Q(char *init_conditions_file, Mat2D Q, double gamma, int N)
{
    char current_word[MAXWORD], current_line[MAXWORD];
    float norm_rho, norm_u, norm_e, norm_p;
    float temp_d;

    FILE *fp = fopen(init_conditions_file, "rt");
    if (fp == NULL) {
        fprintf(stderr, "%s:%d:[Error] problem opening input file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }

    for (int i = 0; i < 5; i++) {
        fscanf(fp, "%s", current_word);
    }

    for (int i = 1; i <= N; i++) {
        // fscanf(fp, "%g %g %g %g", &temp_d, &norm_rho, &norm_u, &norm_p);
        // printf("%g, %g, %g, %g\n", temp_d, norm_rho, norm_u, norm_p);


        fscanf(fp, "%s", current_word);
        double x = atof(current_word);
        fscanf(fp, "%s", current_word);
        norm_rho = atof(current_word);
        fscanf(fp, "%s", current_word);
        norm_u = atof(current_word);
        fscanf(fp, "%s", current_word);
        norm_p = atof(current_word);
        printf("%g, %g, %g, %g\n", x, norm_rho, norm_u, norm_p);

        norm_e = calc_norm_energy(gamma, norm_rho, norm_u, norm_p);

        Q.elements[i] = norm_rho;
        (Q).elements[mat2D_offset2d((Q), (0), (i))];
        MAT2D_AT(Q, 0, i) = norm_rho;
        MAT2D_AT(Q, 1, i) = norm_rho * norm_u;
        MAT2D_AT(Q, 2, i) = norm_e;
    }

    fclose(fp);
}

void calc_vector_of_E(Mat2D E, Mat2D Q, double gamma, int N)
{
    double norm_rho, norm_u, norm_e, norm_p;

    for (int i = 0; i < N+2; i++) {
        norm_rho = MAT2D_AT(Q, 0, i);
        norm_u   = MAT2D_AT(Q, 1, i) / norm_rho;
        norm_e   = MAT2D_AT(Q, 2, i);
        norm_p   = calc_norm_p(gamma, norm_rho, norm_u, norm_e);

        MAT2D_AT(E, 0, i) = norm_rho * norm_u;
        MAT2D_AT(E, 1, i) = norm_p + norm_rho * norm_u * norm_u;
        MAT2D_AT(E, 2, i) = (norm_e + norm_p) * norm_u;
    }
}
