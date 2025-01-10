#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>

#define MAXDIR 1000
#define MAXWORD MAXDIR
#define PI M_PI
#define dprintSTRING(expr) do{printf("%d: ", __LINE__); printf(#expr " = %s\n", expr);} while(0)   /* macro for easy debuging*/
#define dprintINT(expr) do{printf("%d: ", __LINE__); printf(#expr " = %d\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintF(expr) do{printf("%d: ", __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintD(expr) do{printf("%d: ", __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/

typedef enum {
    ROE_FIRST  = (1 << 0),
    ROE_SECOND = (1 << 1),
    NO_LIMITER = (1 << 2),
    VAN_ALBADA = (1 << 3),
    SUPERBEE   = (1 << 4),
    VAN_LEER   = (1 << 5),
    MINMOD     = (1 << 6),
} t_flag;

#ifdef __linux__
    #include <sys/stat.h>
    #include <dirent.h>
    #define ON_LINUX 1
    int create_empty_dir(char *parent_directory);
#endif

void read_input(char *input_file, int *N, double *u1, double *x_max, double *x_min, double *k, t_flag *flags, double *CFL, double *delta_x, int *iterations);
void init(double *u, double u1, double delta_x, int N);
void print_array(double *array, int start, int end);
void print_array_to_file(FILE *fp, double *array, int start, int end);
void make_u_physical(double *u_bar_i_plus_half, double u_i, double u_bar_i_plus_1);
double calculate_F(double u);
double calculate_Roe_first_f_i_plus_half(double u_i, double u_i_plus_1);
double limiter(double r, t_flag flags);
double calculate_Roe_second_f_i_plus_half(double u_i_minus_1, double u_i, double u_i_plus_1, double u_i_plus_2, double k, t_flag flags);
void calculate_vec_f(double *f, double *u, int N, double k, t_flag flags);
double calculate_delta_time(double *u, double delta_x, double CFL, int N);
void calculate_delta_u(double *f, double *delta_u, double delta_time, double delta_x, int N);
void update_u(double *u, double *delta_u, int N);
double calculate_norm(double *delta_u, int start, int end);
void output(char *output_dir, int N, double u1, double x_max, double x_min, t_flag flags, double CFL, double delta_x);

int main(int argc, char const *argv[])
{
/* declerations */
    char input_file[MAXDIR], output_dir[MAXDIR], temp_word[MAXWORD];
    double *u, *f, *delta_u, u1, x_max, x_min, CFL, delta_x, delta_time, first_delta_u_norm, current_delta_u_norm, k;
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
    read_input(input_file, &N, &u1, &x_max, &x_min, &k, &flags, &CFL, &delta_x, &iterations);

/* Checking the input */
    dprintINT(N);
    dprintD(u1);
    dprintD(x_max);
    dprintD(x_min);
    dprintD(k);
    dprintD(delta_x);
    dprintD(CFL);
    dprintINT(iterations);
    dprintINT(flags);
    printf("--------------------\n");

    if (ON_LINUX) {
    /* creating output directory */
        if (create_empty_dir(output_dir) != 0) {
            fprintf(stderr, "%s:%d: [Error] creating ouput directory\n", __FILE__, __LINE__);
            return -1;
        }
        sprintf(temp_word, "_CFL%g", CFL);
        if (flags & ROE_FIRST) {
            strcat(output_dir, "/Roe_first");
        }
        if (flags & ROE_SECOND) {
            if (flags & NO_LIMITER) {
                strcat(output_dir, "/Roe_second_no_limiter");
            }
            if (flags & VAN_ALBADA) {
                strcat(output_dir, "/Roe_second_van_Albada");
            }
            if (flags & SUPERBEE) {
                strcat(output_dir, "/Roe_second_superbee");
            }
            if (flags & VAN_LEER) {
                strcat(output_dir, "/Roe_second_van_Leer");
            }
            if (flags & MINMOD) {
                strcat(output_dir, "/Roe_second_minmod");
            }
        }
        strcat(output_dir, temp_word);
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
    f = (double *)malloc(sizeof(double) * (N+1));
    for (int i = 0; i < N+1; i++) {
        f[i] = 0;
    }

/*------------------------------------------------------------*/
/* initializtion */
    init(u, u1, delta_x, N);
    print_array(u, 0, N+1);
    print_array_to_file(output_u_file, u, 0, N+1);
    fprintf(output_iter_file, "%s, %s, %s\n", "No", "norm", "delta_time");
/*------------------------------------------------------------*/
/* the loop */


    for (int iter = 0; iter < iterations; iter++) {
        calculate_vec_f(f, u, N, k, flags);
        delta_time = calculate_delta_time(u, delta_x, CFL, N);
        calculate_delta_u(f, delta_u, delta_time, delta_x, N);
        if (iter == 0) {
            first_delta_u_norm = calculate_norm(delta_u, 1, N);
        }
        current_delta_u_norm = calculate_norm(delta_u, 1, N);

        // print_array(delta_u, 1, N);

        update_u(u, delta_u, N);

        // print_array(u, 0, N+1);

        print_array_to_file(output_u_file, u, 0, N+1);
        fprintf(output_iter_file, "%d, %g, %g\n", iter, current_delta_u_norm+1, delta_time);

        printf("%d: %g\n", iter, current_delta_u_norm);

        if (current_delta_u_norm/first_delta_u_norm < 1e-6 || isinf(current_delta_u_norm) || isnan(current_delta_u_norm)) {
            break;
        }
    }

    
/*------------------------------------------------------------*/
/* output */
    output(output_dir, N, u1, x_max, x_min, flags, CFL, delta_x);

/*------------------------------------------------------------*/
/* freeing the memory */
    free(u); 

    fclose(output_u_file);

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
void read_input(char *input_file, int *N, double *u1, double *x_max, double *x_min, double *k, t_flag *flags, double *CFL, double *delta_x, int *iterations)
{
    char current_word[MAXWORD];
    float temp_f;
    int temp_i;
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
        } else if (!strcmp(current_word, "u1")) {
            fscanf(fp, "%g", &temp_f);
            *u1 = (double)temp_f;
        } else if (!strcmp(current_word, "x_max")) {
            fscanf(fp, "%g", &temp_f);
            *x_max = (double)temp_f;
        } else if (!strcmp(current_word, "x_min")) {
            fscanf(fp, "%g", &temp_f);
            *x_min = (double)temp_f;
        } else if (!strcmp(current_word, "k")) {
            fscanf(fp, "%g", &temp_f);
            *k = (double)temp_f;
        } else if (!strcmp(current_word, "Roe_first")) {
            fscanf(fp, "%d", &temp_i);
            if (temp_i)
                *flags |= ROE_FIRST;
        } else if (!strcmp(current_word, "Roe_second")) {
            fscanf(fp, "%d", &temp_i);
            if (temp_i)
                *flags |= ROE_SECOND;
        } else if (!strcmp(current_word, "Limiter")) {
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
        } else if (!strcmp(current_word, "iterations")) {
            fscanf(fp, "%g", &temp_f);
            *iterations = (int)temp_f;
        }
    }

    if (*flags & ROE_SECOND || *flags & ROE_FIRST) { 
        *delta_x = (double)(*x_max-*x_min)/(double)(*N-1);
    }
    if (*flags & ROE_FIRST && *flags & ROE_SECOND) {
        exit(1);
    }

    fclose(fp);
}

/* initializing the solution vector with ones and the boundry conditions
argument list:
u       - pointer to the solution array
u1      - boundry condition at x_max
delta_x - double value of delta_x
N       - number of grid points */
void init(double *u, double u1, double delta_x, int N)
{
    for (int i = 1; i <= N; i++) {
        u[i] = 1 - (1 - u1) * delta_x * (i-1);
    }
    u[0] = 1;
    u[N+1] = u1;
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
double calculate_F(double u)
{
    return u * u / 2;
}

/* calculating the flax for a face
argument list: */
double calculate_Roe_first_f_i_plus_half(double u_i, double u_i_plus_1)
{
    double u_bar_i_plus_half = (u_i + u_i_plus_1) * 0.5;
    // make_u_physical(&u_bar_i_plus_half, u_i, u_i_plus_1);

    double u_bar_i_plus_half_plus = (u_bar_i_plus_half + fabs(u_bar_i_plus_half)) * 0.5;
    double u_bar_i_plus_half_minus = (u_bar_i_plus_half - fabs(u_bar_i_plus_half)) * 0.5;
    
    double F_i = calculate_F(u_i);
    double F_i_plus_1 = calculate_F(u_i_plus_1);

    double f_i_plus_half = (F_i + F_i_plus_1) * 0.5 - 0.5 * (u_bar_i_plus_half_plus - u_bar_i_plus_half_minus) * ( u_i_plus_1 - u_i);

    return f_i_plus_half;
}

/* calculating the van Albada limiter
argument list:*/
double limiter(double r, t_flag flags)
{
    if (isnan(r) || isinf(r)) {
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
double calculate_Roe_second_f_i_plus_half(double u_i_minus_1, double u_i, double u_i_plus_1, double u_i_plus_2, double k, t_flag flags)
{
    double r_plus = (u_i_plus_2 - u_i_plus_1) / (u_i_plus_1 - u_i);
    double r_minus = (u_i - u_i_minus_1) / (u_i_plus_1 - u_i);
    double psi_plus = limiter(r_plus, flags);
    double psi_minus = limiter(r_minus, flags);

    // printf("%g, %g, %g, %g\n", r_plus, r_minus, psi_plus, psi_minus);

    double delta_u_i_minus_half = u_i - u_i_minus_1;
    double delta_u_i_plus_half = u_i_plus_1 - u_i;
    double delta_u_i_plus_three_half = u_i_plus_2 - u_i_plus_1;

    double u_left_i_plus_half = u_i + (1-k)/4 * psi_minus * delta_u_i_minus_half + (1+k)/4 * psi_plus * delta_u_i_plus_half;
    double u_right_i_plus_half = u_i_plus_1 - (1+k)/4 * psi_plus * delta_u_i_plus_half - (1-k)/4 * psi_minus * delta_u_i_plus_three_half;

    // double u_left_i_plus_half = u_i + (1-k)/4 * psi_plus * delta_u_i_minus_half;
    // double u_right_i_plus_half = u_i_plus_1 - (1-k)/4 * psi_minus * delta_u_i_plus_three_half;

    double F_left = calculate_F(u_left_i_plus_half);
    double F_right = calculate_F(u_right_i_plus_half);
    
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
void calculate_vec_f(double *f, double *u, int N, double k, t_flag flags)
{
    if (flags & ROE_FIRST) {
        for (int i = 0; i < N+1; i++) {
            f[i] = calculate_Roe_first_f_i_plus_half(u[i], u[i+1]);
        }
    }
    if (flags & ROE_SECOND) {
        for (int i = 1; i <= N-1; i++) {
            f[i] = calculate_Roe_second_f_i_plus_half(u[i-1], u[i], u[i+1], u[i+2], k, flags);
        }
        f[0] = calculate_Roe_second_f_i_plus_half(u[0], u[0], u[1], u[2], k, flags);
        f[N] = calculate_Roe_second_f_i_plus_half(u[N-1], u[N], u[N+1], u[N+1], k, flags);

        // f[0] = u[0];
        // f[N] = u[N+1];
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

void calculate_delta_u(double *f, double *delta_u, double delta_time, double delta_x, int N)
{
    for (int i = 1; i <= N; i++) {
        delta_u[i] = delta_time / delta_x * (f[i] - f[i-1]);
    }
}

void update_u(double *u, double *delta_u, int N)
{
    for (int i = 1; i < N+1; i++) {
        u[i] = u[i] - delta_u[i];
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
void output(char *output_dir, int N, double u1, double x_max, double x_min, t_flag flags, double CFL, double delta_x)
{
    char temp_word[MAXWORD];
    FILE *meta_data_file;

    strcpy(temp_word, output_dir);
    strcat(temp_word, "/mata_data.txt");
    meta_data_file = fopen(temp_word, "wt");

    fprintf(meta_data_file, "%s, %s, %s, %s, %s, %s, %s, %s\n", "N", "u1", "x_max", "x_min", "Roe_first", "Roe_second", "CFL", "delta_x");
    fprintf(meta_data_file, "%d, %g, %g, %g, %d, %d, %g, %g", N, u1, x_max, x_min, flags & ROE_FIRST, flags & ROE_SECOND, CFL, delta_x);


    fclose(meta_data_file);
}
