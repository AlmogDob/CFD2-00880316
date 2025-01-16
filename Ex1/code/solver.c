#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>

#define MAXDIR 1000
#define MAXWORD MAXDIR
#define PI M_PI
#define dprintSTRING(expr) do{printf("%s:%d: ", __FILE__, __LINE__); printf(#expr " = %s\n", expr);} while(0)   /* macro for easy debuging*/
#define dprintINT(expr) do{printf("%s:%d: ", __FILE__, __LINE__); printf(#expr " = %d\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintF(expr) do{printf("%s:%d: ", __FILE__, __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintD(expr) do{printf("%s:%d: ", __FILE__, __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/

typedef enum {
    ROE_FIRST        = (1 << 0),
    ROE_SECOND       = (1 << 1),
    NO_LIMITER       = (1 << 2),
    VAN_ALBADA       = (1 << 3),
    SUPERBEE         = (1 << 4),
    VAN_LEER         = (1 << 5),
    MINMOD           = (1 << 6),
    INVISCID         = (1 << 7),
    GENERAL          = (1 << 8), 
    MACCORMACK       = (1 << 9),
    BEAM_AND_WARMING = (1 << 10),
} t_flag;

#ifdef __linux__
    #include <sys/stat.h>
    #include <dirent.h>
    #define ON_LINUX 1
    int create_empty_dir(char *parent_directory);
#endif

void read_input(char *input_file, int *N, double *u0, double *u1, double *x_max, double *x_min, double *k, double *b, double *c, double *mu, t_flag *flags, double *CFL, double *w, double *theta, double *delta_x, double *delta_time, int *iterations, double *final_time);
void init(double *u, double u0, double u1, double x_max, double x_min, double delta_x, int N, t_flag flags);
void set_ghost_cell(double *u, double u0, double u1, int N);
void print_array(double *array, int start, int end);
void print_array_to_file(FILE *fp, double *array, int start, int end);
void make_u_physical(double *u_bar_i_plus_half, double u_i, double u_bar_i_plus_1);
double calculate_F(double u, double b, double c);
double calculate_Roe_first_f_i_plus_half(double u_i, double u_i_plus_1, double b, double c);
double limiter(double r, t_flag flags);
double calculate_Roe_second_f_i_plus_half(double u_i_minus_1, double u_i, double u_i_plus_1, double u_i_plus_2, double k, double b, double c, t_flag flags);
void calculate_vec_f(double *f, double *u, int N, double k, double b, double c, t_flag flags);
void calculate_delta_u_MacCormack(double *delta_u, double *u, double *u_bar, double delta_time, double delta_x, double mu, double b, double c, int N);
double calculate_delta_time(double *u, double delta_x, double CFL, int N);
void RHS(double *D, double *u, double delta_x, double delta_time, double mu, double b, double c, double w, int N);
void LHS(double *A, double *B, double *C, double *u, double delta_time, double delta_x, double mu, double b, double c, double theta, int N);
void BC(double *A, double *C, double *D, double *u, double u0, double u1, int N);
int tridiag(double *a, double *b, double *c, double *d, double *u, int is, int ie);
void calculate_delta_u(double *f, double *u, double *delta_u, double *u_bar, double *A, double *B, double *C, double *D, double delta_time, double delta_x, double b, double c, double mu, double w, double theta, int N, t_flag flags);
void update_u(double *u, double *delta_u, int N);
double calculate_norm(double *delta_u, int start, int end);
void output(char *output_dir, int N, double u0, double u1, double x_max, double x_min, double delta_x, double delta_time, double k, double b, double c, double mu,  double CFL, double w, double theta, int iterations, double final_time, t_flag flags);

int main(int argc, char const *argv[])
{
/* declerations */
    char input_file[MAXDIR], output_dir[MAXDIR], temp_word[MAXWORD];
    double *u, *f, *delta_u, *u_bar, *A, *B, *C, *D, u0, u1, x_max, x_min, CFL, delta_x, delta_time, first_delta_u_norm, current_delta_u_norm, k, c, b, mu, w, theta, ellapsed_time=0, final_time;
    int N, iterations;
    FILE *output_u_file, *output_iter_file;  
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
    read_input(input_file, &N, &u0, &u1, &x_max, &x_min, &k, &b, &c, &mu, &flags, &CFL, &w, &theta, &delta_x, &delta_time, &iterations, &final_time);

/* Checking the input */
    dprintINT(N);
    dprintD(u0);
    dprintD(u1);
    dprintD(x_max);
    dprintD(x_min);
    dprintD(k);
    dprintD(c);
    dprintD(b);
    dprintD(mu);
    dprintD(delta_x);
    dprintD(delta_time);
    dprintD(CFL);
    dprintD(w);
    dprintD(theta);
    dprintINT(iterations);
    dprintD(final_time);
    dprintINT(flags);
    printf("--------------------\n");

/*------------------------------------------------------------*/
/* creating output directory */
    if (ON_LINUX) {
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
            return -1;
        }
        if (flags & INVISCID) {
            strcat(output_dir, "/inviscid");
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
            sprintf(temp_word, "/u0_%g_u1_%g", u0, u1);
            strcat(output_dir, temp_word);
            if (create_empty_dir(output_dir) != 0) {
                fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
                return -1;
            }
            if (flags & ROE_FIRST) {
                strcat(output_dir, "/Roe_first");
            } else if (flags & ROE_SECOND) {
                if (flags & NO_LIMITER) {
                    strcat(output_dir, "/Roe_second_no_limiter");
                } else if (flags & VAN_ALBADA) {
                    strcat(output_dir, "/Roe_second_van_Albada");
                } else if (flags & SUPERBEE) {
                    strcat(output_dir, "/Roe_second_superbee");
                } else if (flags & VAN_LEER) {
                    strcat(output_dir, "/Roe_second_van_Leer");
                } else if (flags & MINMOD) {
                    strcat(output_dir, "/Roe_second_minmod");
                }
            }
            sprintf(temp_word, "_CFL%g", CFL);
            strcat(output_dir, temp_word);
        } else if (flags & GENERAL) {
            strcat(output_dir, "/general");
            if (create_empty_dir(output_dir) != 0) {
                fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
                return -1;
            }
            if (flags & MACCORMACK) {
                sprintf(temp_word, "/MacCormack");
                strcat(output_dir, temp_word);
                if (create_empty_dir(output_dir) != 0) {
                    fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
                    return -1;
                }
            } else if (flags & ROE_FIRST) {
                sprintf(temp_word, "/Roe_first");
                strcat(output_dir, temp_word);
                if (create_empty_dir(output_dir) != 0) {
                    fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
                    return -1;
                }
            } else if (flags & BEAM_AND_WARMING) {
                sprintf(temp_word, "/Beam_and_Warming");
                strcat(output_dir, temp_word);
                if (create_empty_dir(output_dir) != 0) {
                    fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
                    return -1;
                }
                sprintf(temp_word, "/theta%g", theta);
                strcat(output_dir, temp_word);
                if (create_empty_dir(output_dir) != 0) {
                    fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
                    return -1;
                }
                sprintf(temp_word, "/w%g", w);
                strcat(output_dir, temp_word);
                if (create_empty_dir(output_dir) != 0) {
                    fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
                    return -1;
                }
            }
            sprintf(temp_word, "/mu%g", mu);
            strcat(output_dir, temp_word);
            if (create_empty_dir(output_dir) != 0) {
                fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
                return -1;
            }
            sprintf(temp_word, "/delta_time%g", delta_time);
            strcat(output_dir, temp_word);
            if (create_empty_dir(output_dir) != 0) {
                fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
                return -1;
            }
            sprintf(temp_word, "/u0_%g_u1_%g", u0, u1);
            strcat(output_dir, temp_word);
            if (create_empty_dir(output_dir) != 0) {
                fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
                return -1;
            }
        }
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
            return -1;
        }

        strcpy(temp_word, output_dir);
        strcat(temp_word, "/output_u.txt");
        output_u_file = fopen(temp_word, "wt");

        strcpy(temp_word, output_dir);
        strcat(temp_word, "/output_iter.txt");
        output_iter_file = fopen(temp_word, "wt");
    }

/*------------------------------------------------------------*/
/* allocating the matrices */

    u = (double *)malloc(sizeof(double) * (N+2));
    for (int i = 0; i < N+2; i++) {
        u[i] = 0;
    }
    delta_u = (double *)malloc(sizeof(double) * (N+2));
    for (int i = 0; i < N+2; i++) {
        delta_u[i] = 0;
    }
    u_bar = (double *)malloc(sizeof(double) * (N+2));
    for (int i = 0; i < N+2; i++) {
        u_bar[i] = 0;
    }
    A = (double *)malloc(sizeof(double) * (N+2));
    for (int i = 0; i < N+2; i++) {
        A[i] = 0;
    }
    B = (double *)malloc(sizeof(double) * (N+2));
    for (int i = 0; i < N+2; i++) {
        B[i] = 0;
    }
    C = (double *)malloc(sizeof(double) * (N+2));
    for (int i = 0; i < N+2; i++) {
        C[i] = 0;
    }
    D = (double *)malloc(sizeof(double) * (N+2));
    for (int i = 0; i < N+2; i++) {
        D[i] = 0;
    }
    f = (double *)malloc(sizeof(double) * (N+1));
    for (int i = 0; i < N+1; i++) {
        f[i] = 0;
    }

/*------------------------------------------------------------*/
/* initializtion */
    init(u, u0, u1, x_max, x_min, delta_x, N, flags);
    print_array(u, 0, N+1);
    print_array_to_file(output_u_file, u, 0, N+1);
    fprintf(output_iter_file, "%s, %s, %s\n", "No", "norm", "delta_time");
/*------------------------------------------------------------*/
/* the loop */
    dprintINT(flags & BEAM_AND_WARMING);

    for (int iter = 0; iter < iterations; iter++) {
        set_ghost_cell(u, u0, u1, N);
        if (flags & ROE_FIRST || flags & ROE_SECOND) {
            calculate_vec_f(f, u, N, k, b, c, flags);
        }
        if (flags & INVISCID) {
            delta_time = calculate_delta_time(u, delta_x, CFL, N);
        }
        ellapsed_time += delta_time;

        calculate_delta_u(f, u, delta_u, u_bar, A, B, C, D, delta_time, delta_x, b, c, mu, w, theta, N, flags);

        if (iter == 0) {
            first_delta_u_norm = calculate_norm(delta_u, 1, N);
            current_delta_u_norm = calculate_norm(delta_u, 1, N);
        }
        current_delta_u_norm = calculate_norm(delta_u, 1, N);

        update_u(u, delta_u, N);

        print_array_to_file(output_u_file, u, 0, N+1);
        fprintf(output_iter_file, "%d, %g, %g\n", iter+1, current_delta_u_norm, delta_time);
        printf("%4d: %14g\t, time: %g\n", iter, (float) current_delta_u_norm, ellapsed_time);

        if (current_delta_u_norm/first_delta_u_norm < 1e-6 || isinf(current_delta_u_norm) || isnan(current_delta_u_norm) || (flags & GENERAL && ellapsed_time > final_time)) {
            break;
        }
    }
    // dprintINT(flags & ROE_FIRST);
    
/*------------------------------------------------------------*/
/* output */
    output(output_dir, N, u0, u1, x_max, x_min, delta_x, delta_time, k, b, c, mu, CFL, w, theta, iterations, final_time, flags);

/*------------------------------------------------------------*/
/* freeing the memory */
    free(u); 
    free(u_bar); 
    free(delta_u); 
    free(f); 
    free(A);
    free(B);
    free(C);
    free(D);

    fclose(output_u_file);
    fclose(output_iter_file);

    return 0;
}

#if ON_LINUX
/* create empty dir at 'parent directory'.
if allready exisest, delet all the files inside
returns 0 on success
this functin handls the errors so on fail just quit
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
input-file - file pointer to input file
N          - int pointer
u1         - double pointer
x_max      - double pointer
x_min      - double pointer
flags      - bit flags
CFL        - double pointer
delta_x    - double pointer 
iterations - max number of desierd iterations */
void read_input(char *input_file, int *N, double *u0, double *u1, double *x_max, double *x_min, double *k, double *b, double *c, double *mu, t_flag *flags, double *CFL, double *w, double *theta, double *delta_x, double *delta_time, int *iterations, double *final_time)
{
    char current_word[MAXWORD];
    float temp_f;
    FILE *fp = fopen(input_file, "rt");
    if (fp == NULL) {
        fprintf(stderr, "%s:%d:[Error] problem opening input file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }

    *flags |= NO_LIMITER;

    while(fscanf(fp, "%s", current_word) != EOF) {
        if (!strcmp(current_word, "N")) {
            fscanf(fp, "%d", N);
            *N = *N;
        } else if (!strcmp(current_word, "u0")) {
            fscanf(fp, "%g", &temp_f);
            *u0 = (double)temp_f;
        } else if (!strcmp(current_word, "u1")) {
            fscanf(fp, "%g", &temp_f);
            *u1 = (double)temp_f;
        } else if (!strcmp(current_word, "x_max")) {
            fscanf(fp, "%g", &temp_f);
            *x_max = (double)temp_f;
        } else if (!strcmp(current_word, "x_min")) {
            fscanf(fp, "%g", &temp_f);
            *x_min = (double)temp_f;
        } else if (!strcmp(current_word, "delta_x")) {
            fscanf(fp, "%g", &temp_f);
            *delta_x = (double)temp_f;
        } else if (!strcmp(current_word, "delta_time")) {
            fscanf(fp, "%g", &temp_f);
            *delta_time = (double)temp_f;
        } else if (!strcmp(current_word, "k")) {
            fscanf(fp, "%g", &temp_f);
            *k = (double)temp_f;
        } else if (!strcmp(current_word, "c")) {
            fscanf(fp, "%g", &temp_f);
            *c = (double)temp_f;
        } else if (!strcmp(current_word, "b")) {
            fscanf(fp, "%g", &temp_f);
            *b = (double)temp_f;
        } else if (!strcmp(current_word, "mu")) {
            fscanf(fp, "%g", &temp_f);
            *mu = (double)temp_f;
            if (*mu != 0) {
                *flags |= GENERAL;
            } else if (*mu == 0) {
                *flags |= INVISCID;
            }
        } else if (!strcmp(current_word, "method")) {
            fscanf(fp, "%s", current_word);
                if (!strcmp(current_word, "Roe_first")) {
                    *flags |= ROE_FIRST;
                } else if (!strcmp(current_word, "Roe_second")) {
                    *flags |= ROE_SECOND;
                } else if (!strcmp(current_word, "MacCormack")) {
                    *flags |= MACCORMACK;
                } else if (!strcmp(current_word, "Beam_and_Warming")) {
                    *flags |= BEAM_AND_WARMING;
                }
        } else if (!strcmp(current_word, "limiter")) {
            fscanf(fp, "%s", current_word);
                if (!strcmp(current_word, "van_Albada")) {
                    *flags &= ~ NO_LIMITER;
                    *flags |= VAN_ALBADA;
                } else if (!strcmp(current_word, "superbee")) {
                    *flags &= ~ NO_LIMITER;
                    *flags |= SUPERBEE;
                } else if (!strcmp(current_word, "van_Leer")) {
                    *flags &= ~ NO_LIMITER;
                    *flags |= VAN_LEER;
                } else if (!strcmp(current_word, "minmod")) {
                    *flags &= ~ NO_LIMITER;
                    *flags |= MINMOD;
                }
        } else if (!strcmp(current_word, "CFL")) {
            fscanf(fp, "%g", &temp_f);
            *CFL = (double)temp_f;
        } else if (!strcmp(current_word, "w")) {
            fscanf(fp, "%g", &temp_f);
            *w = (double)temp_f;
        } else if (!strcmp(current_word, "theta")) {
            fscanf(fp, "%g", &temp_f);
            *theta = (double)temp_f;
        } else if (!strcmp(current_word, "iterations")) {
            fscanf(fp, "%g", &temp_f);
            *iterations = (int)temp_f;
        } else if (!strcmp(current_word, "final_time")) {
            fscanf(fp, "%g", &temp_f);
            *final_time = (double)temp_f;
            if (*final_time <= 0.1) {
                *final_time = __DBL_MAX__;
            }
        }
    }

    if (*flags & INVISCID) { 
        *delta_x = (double)(*x_max-*x_min)/(double)(*N-1);
    }
    if (*flags & ROE_FIRST && *flags & ROE_SECOND && *flags & MACCORMACK) {
        fprintf(stderr, "%s:%d: [Error] unable to do more than one method\n", __FILE__, __LINE__);
        exit(1);
    }
    if (*flags & INVISCID && *flags & GENERAL) {
        fprintf(stderr, "%s:%d: [Error] defined both invicsid and general Burger's equation\n", __FILE__, __LINE__);
        exit(2);
    }

    fclose(fp);
}

/* initializing the solution vector with ones and the boundry conditions
argument list:
u       - pointer to the solution array
u1      - boundry condition at x_max
delta_x - double value of delta_x
N       - number of grid points */
void init(double *u, double u0, double u1, double x_max, double x_min, double delta_x, int N, t_flag flags)
{
    if (flags & INVISCID) {
        for (int i = 1; i <= N; i++) {
            u[i] = (x_max - x_min) - (u0 - u1) * delta_x * (i-1);
        }
        u[0] = u[1];
        u[N+1] = u[N];
    } else if (flags & GENERAL) {
        for (int i = 1; i <= N; i++) {
            u[i] = u0 + 0.5 * ((u1-u0) + (u1 - u0) * tanh(250 * (delta_x * (i-1) - 20)));
        }
        u[0] = u[1];
        u[N+1] = u[N];
    }
}

void set_ghost_cell(double *u, double u0, double u1, int N)
{
    u[0] = -u[1] + u0 * 2;
    u[N+1] = -u[N] + u1 * 2;
}

/* printing a double array with length 'len' to 'stdin'
argument list:
array - pointer to the array
start - start index in the array
end   - end index in the array */
void print_array(double *array, int start, int end)
{
    printf("--------------------\n");
    for (int i = start; i <= end; i++) {
        printf("%g ", array[i]);
    }
    printf("\n");
    printf("--------------------\n");
}

/* checking if the inputed delta t will result in converging solution
arugemnt list: 
delta_time - time step between iterations
delta_x    - the distance between grid points */
void check_delta_time(double delta_time, double delta_x)
{
    if (delta_time/delta_x > 1) {
        printf("%s:%d:[Warning] delta_time too big\n", __FILE__, __LINE__);
    }
}

/* printing a double array with length 'len' to file 'fp' 
argument list:
fp    - file pointer
array - pointer to the array
start - start index in the array
end   - end index in the array */
void print_array_to_file(FILE *fp, double *array, int start, int end)
{
    for (int i = start; i <= end; i++) {
        fprintf(fp, "%g ", array[i]);
    }
    fprintf(fp, "\n");
}

void make_u_physical(double *u_bar_i_plus_half, double u_i, double u_i_plus_1)
{
    double epsilon = fmax(0, (u_i_plus_1 - u_i) / 2); 

    if (*u_bar_i_plus_half < epsilon) {
        *u_bar_i_plus_half = epsilon;
    }
}

/* calculating F
argument list: */
double calculate_F(double u, double b, double c)
{
    return c * u + b * u * u * 0.5;
}

/* calculating the flax for a face
argument list: */
double calculate_Roe_first_f_i_plus_half(double u_i, double u_i_plus_1, double b, double c)
{
    double A_bar_i_plus_half, A_bar_i_plus_half_plus, A_bar_i_plus_half_minus, F_i, F_i_plus_1,  f_i_plus_half;

    if (u_i == u_i_plus_1) {
        A_bar_i_plus_half = c + b * u_i;
    } else {
        A_bar_i_plus_half = (c * (u_i_plus_1 - u_i) + 0.5 * b * (u_i_plus_1 * u_i_plus_1 - u_i * u_i)) / (u_i_plus_1 - u_i);
    }
    A_bar_i_plus_half_plus = (A_bar_i_plus_half + fabs(A_bar_i_plus_half)) * 0.5;
    A_bar_i_plus_half_minus = (A_bar_i_plus_half - fabs(A_bar_i_plus_half)) * 0.5;

    F_i = calculate_F(u_i, b, c);
    F_i_plus_1 = calculate_F(u_i_plus_1, b, c);

    f_i_plus_half = (F_i + F_i_plus_1) * 0.5 - 0.5 * (A_bar_i_plus_half_plus - A_bar_i_plus_half_minus) * ( u_i_plus_1 - u_i);

    return f_i_plus_half;
}

/* calculating the van Albada limiter
argument list:*/
double limiter(double r, t_flag flags)
{
    if (!isfinite(r)) {
        return 1;
    } else if (flags & NO_LIMITER) {
        return 1;
    } else if (flags & VAN_ALBADA) {
        return (r + r * r) / (1 + r * r);
    } else if (flags & SUPERBEE) {
        return fmax(fmax(0,fmin(2*r,1)),fmin(r,2));
    } else if (flags & VAN_LEER) {
        return (r + fabs(r)) / (1 + r * r);
    } else if (flags & MINMOD) {
        if (1 * r > 0) {
            return 1 > fabs(r)?r:1;
        } else {
            return 0;
        }
    }
    return 1;
}

/* calculating the flax for a face
argument list:*/
double calculate_Roe_second_f_i_plus_half(double u_i_minus_1, double u_i, double u_i_plus_1, double u_i_plus_2, double k, double b, double c, t_flag flags)
{
    double r_plus_i_minus_half = (u_i_plus_1 - u_i) / (u_i - u_i_minus_1);

    double r_minus_i_plus_half = (u_i - u_i_minus_1) / (u_i_plus_1 - u_i);
    double r_plus_i_plus_half = (u_i_plus_2 - u_i_plus_1) / (u_i_plus_1 - u_i);

    double r_minus_i_plus_three_half = (u_i_plus_1 - u_i) / (u_i_plus_2 - u_i_plus_1);


    double psi_plus_i_minus_half = limiter(r_plus_i_minus_half, flags);

    double psi_minus_i_plus_half = limiter(r_minus_i_plus_half, flags);
    double psi_plus_i_plus_half = limiter(r_plus_i_plus_half, flags);

    double psi_minus_i_plus_three_half = limiter(r_minus_i_plus_three_half, flags);

    double delta_u_i_minus_half = u_i - u_i_minus_1;
    double delta_u_i_plus_half = u_i_plus_1 - u_i;
    double delta_u_i_plus_three_half = u_i_plus_2 - u_i_plus_1;

    double u_left_i_plus_half = u_i + (1-k)/4 * psi_plus_i_minus_half * delta_u_i_minus_half + (1+k)/4 * psi_minus_i_plus_half * delta_u_i_plus_half;
    double u_right_i_plus_half = u_i_plus_1 - (1+k)/4 * psi_plus_i_plus_half * delta_u_i_plus_half - (1-k)/4 * psi_minus_i_plus_three_half * delta_u_i_plus_three_half;

    double F_left = calculate_F(u_left_i_plus_half, b, c);
    double F_right = calculate_F(u_right_i_plus_half, b, c);
    
    double u_bar_i_plus_half = (u_left_i_plus_half + u_right_i_plus_half) * 0.5;

    // make_u_physical(&u_bar_i_plus_half, u_i, u_i_plus_1);

    double u_bar_plus_i_plus_half = 0.5 * (u_bar_i_plus_half + fabs(u_bar_i_plus_half));
    double u_bar_minus_i_plus_half = 0.5 * (u_bar_i_plus_half - fabs(u_bar_i_plus_half));

    double f_i_plus_half = 0.5 * (F_left + F_right - (u_bar_plus_i_plus_half - u_bar_minus_i_plus_half) * ( u_right_i_plus_half - u_left_i_plus_half));
    return f_i_plus_half;
}

/* calculating the flax vector
argument list:
f     - pointer to the flax array
u     - pointer to the u array
N     - number of grid points
flags - bit flags */
void calculate_vec_f(double *f, double *u, int N, double k, double b, double c, t_flag flags)
{
    if (flags & ROE_FIRST) {
        for (int i = 0; i < N+1; i++) {
            f[i] = calculate_Roe_first_f_i_plus_half(u[i], u[i+1], b, c);
        }
    } else if (flags & ROE_SECOND) {
        for (int i = 1; i <= N-1; i++) {
            f[i] = calculate_Roe_second_f_i_plus_half(u[i-1], u[i], u[i+1], u[i+2], k, b, c, flags);
        }
        f[0] = calculate_Roe_second_f_i_plus_half(u[0], u[0], u[1], u[2], k, b, c, flags);
        f[N] = calculate_Roe_second_f_i_plus_half(u[N-1], u[N], u[N+1], u[N+1], k, b, c, flags);
    }
}

void calculate_delta_u_MacCormack(double *delta_u, double *u, double *u_bar, double delta_time, double delta_x, double mu, double b, double c, int N)
{
    for (int i = 1; i <= N; i++) {
        u_bar[i] = u[i] - delta_time * (calculate_F(u[i+1], b, c) - calculate_F(u[i], b, c)) / (delta_x) + mu * delta_time / delta_x / delta_x * (u[i+1] - 2 * u[i] + u[i-1]);
    }

    for (int i = 2; i <= N-1; i++) {
        delta_u[i] = -u[i] + 0.5 * (u[i] + u_bar[i] - delta_time * (calculate_F(u_bar[i], b, c) - calculate_F(u_bar[i-1], b, c)) / (delta_x) + mu * delta_time / delta_x / delta_x * (u_bar[i+1] - 2 * u_bar[i] + u_bar[i-1]));
    }
}

/* calculating the delta time according to the CFL number argument list:
*/
double calculate_delta_time(double *u, double delta_x, double CFL, int N)
{
    double delta_time = __DBL_MAX__;
    for (int i = 1; i <= N+2; i++) {
        if (delta_time > CFL / fabs(u[i]) * delta_x) {
            delta_time = CFL / fabs(u[i]) * delta_x;
        }
    }

    return delta_time;
}

/* calculate the RHS (vector D) 
argumet list:
D - pointer to the array
u - pointer to the flow solution
N - number of grid points */
void RHS(double *D, double *u, double delta_x, double delta_time, double mu, double b, double c, double w, int N)
{
    for (int i = 1; i <= N; i++) {
        D[i] = - delta_time * (calculate_F(u[i+1], b, c) - calculate_F(u[i-1], b, c)) / (2 * delta_x) + delta_time * mu * (u[i+1] - 2*u[i] + u[i-1]) / delta_x / delta_x;
    }
    for (int i = 2; i <= N-1; i++) {
        D[i] += -w / (double)8 * (u[i+2] - 4 * u[i+1] + 6 * u[i] - 4 * u[i-1] + u[i-2]);
    }
    // D[0] += -w / (double)8 * (u[0+2] - 4 * u[0+1] + 6 * u[0] - 4 * u[0] + u[0]);
    // D[1] += -w / (double)8 * (u[1+2] - 4 * u[1+1] + 6 * u[1] - 4 * u[1-1] + u[0]);
    // D[N] += -w / (double)8 * (u[N] - 4 * u[N] + 6 * u[N] - 4 * u[N-1] + u[N-2]);
    // D[N+1] += -w / (double)8 * (u[N+1] - 4 * u[N+1] + 6 * u[N+1] - 4 * u[N+1-1] + u[N+1-2]);
}

/* calculate the LHS (vectors A, B, C)
A          - pointer to the A array
B          - pointer to the B array
C          - pointer to the C array 
delta_time - time step between iterations
delta_y    - distance between to grid points
mu         - viscosity
alpha      - parameter that control the scheme
    - number of grid points */
void LHS(double *A, double *B, double *C, double *u, double delta_time, double delta_x, double mu, double b, double c, double theta, int N)
{
    for (int i = 1; i <= N; i++) {
        A[i] = -theta * (mu * delta_time) / (delta_x * delta_x) - 0.5 * theta * delta_time / delta_x * (c + b * u[i-1]);
        B[i] = 1 + 2 * theta * mu * delta_time / delta_x / delta_x;
        C[i] = -theta * (mu * delta_time) / delta_x / delta_x + 0.5 * theta * delta_time / delta_x * (c + b * u[i+1]);
    }
}

/*   a, b, c, are the vectors of the diagonal and the two off-diagonals.
The vector d is the RHS vector, the vector u is the solution
vector, "is" is the starting point, and ie is
the last point.
*/
int tridiag(double *a, double *b, double *c, double *d, double *u, int is, int ie)
{

  int i;
  double beta;

  for (i = is + 1; i <= ie; i++)
    {
      if(b[i-1] == 0.) return(1);
      beta = a[i] / b[i-1];
      b[i] = b[i] - c[i-1] * beta;
      d[i] = d[i] - d[i-1] * beta;
    }

  u[ie] = d[ie] / b[ie];
  for (i = ie - 1; i >= is; i--)
    {
      u[i] = (d[i] - c[i] * u[i+1]) / b[i];
    }
  return(0);
}

void calculate_delta_u(double *f, double *u, double *delta_u, double *u_bar, double *A, double *B, double *C, double *D, double delta_time, double delta_x, double b, double c, double mu, double w, double theta, int N, t_flag flags)
{
    if (flags & ROE_FIRST || flags & ROE_SECOND) {
        for (int i = 1; i <= N; i++) {
            delta_u[i] = -delta_time / delta_x * (f[i] - f[i-1]) + mu * (delta_time) / (delta_x * delta_x) * (u[i+1] - 2 * u[i] + u[i-1]);
        }
    } else if (flags & MACCORMACK) {
        calculate_delta_u_MacCormack(delta_u, u, u_bar, delta_time, delta_x, mu, b, c, N);
    } else if (flags & BEAM_AND_WARMING) {
        RHS(D, u, delta_x, delta_time, mu, b, c, w, N);
        LHS(A, B, C, u, delta_time, delta_x, mu, b, c, theta, N);
        tridiag(A, B, C, D, delta_u, 1, N);
    }
}

void update_u(double *u, double *delta_u, int N)
{
    for (int i = 1; i <= N; i++) {
        u[i] = u[i] + delta_u[i];
    }
}

/* calculating the second norma of the vector 'delta_u'
argument list:
delta_u - pointer to the vector elements array
start   - start index in the array
end     - end index in the array */
double calculate_norm(double *delta_u, int start, int end)
{
    double sum = 0;
    for (int i = start; i <= end; i++) {
        sum += delta_u[i] * delta_u[i];
    }
    return sqrt(sum);
}

/* ouputing metadata of the solution 
argumetn list:
output_dir - name of the output directory
N          - number of grid points
*/
void output(char *output_dir, int N, double u0, double u1, double x_max, double x_min, double delta_x, double delta_time, double k, double b, double c, double mu,  double CFL, double w, double theta, int iterations, double final_time, t_flag flags)
{
    char temp_word[MAXWORD], method[MAXWORD], limiter_name[MAXWORD];
    FILE *mata_data_file;

    strcpy(temp_word, output_dir);
    strcat(temp_word, "/mata_data.txt");
    mata_data_file = fopen(temp_word, "wt");

    if (flags & ROE_FIRST) {
        strcpy(method, "Roe_first");
    } else if(flags & ROE_SECOND) {
        strcpy(method, "Roe_second");
    } else if(flags & MACCORMACK) {
        strcpy(method, "MacCormack");
    } else if(flags & BEAM_AND_WARMING) {
        strcpy(method, "Beam_and_Warming");
    }
    if (flags & NO_LIMITER) {
        strcpy(limiter_name, "no_limiter");
    } else if (flags & VAN_ALBADA) {
        strcpy(limiter_name, "van_Albada");
    } else if (flags & SUPERBEE) {
        strcpy(limiter_name, "superbee");
    } else if (flags & VAN_LEER) {
        strcpy(limiter_name, "van_Leer");
    } else if (flags & MINMOD) {
        strcpy(limiter_name, "minmod");
    }
    // ROE_FIRST        = (1 << 0),
    // ROE_SECOND       = (1 << 1),
    // NO_LIMITER       = (1 << 2),
    // VAN_ALBADA       = (1 << 3),
    // SUPERBEE         = (1 << 4),
    // VAN_LEER         = (1 << 5),
    // MINMOD           = (1 << 6),
    // INVISCID         = (1 << 7),
    // GENERAL          = (1 << 8), 
    // MACCORMACK       = (1 << 9),
    // BEAM_AND_WARMING = (1 << 10),

    fprintf(mata_data_file, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n" , "N", "u0", "u1", "x_max", "x_min", "delta_x", "delta_time", "k", "b", "c", "mu", "method", "limiter", "CFL", "w", "theta", "iterations", "final_time");
    fprintf(mata_data_file, "%d, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %s, %s, %g, %g, %g, %d, %g" , N, u0, u1, x_max, x_min, delta_x, delta_time, k, b, c, mu, method, limiter_name, CFL, w, theta, iterations, final_time);


    fclose(mata_data_file);
}
