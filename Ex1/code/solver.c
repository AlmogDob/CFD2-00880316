#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <limits.h>

#define MAXDIR 1000
#define MAXWORD MAXDIR
#define PI M_PI
#define dprintSTRING(expr) do{printf("%d: ", __LINE__); printf(#expr " = %s\n", expr);} while(0)   /* macro for easy debuging*/
#define dprintINT(expr) do{printf("%d: ", __LINE__); printf(#expr " = %d\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintF(expr) do{printf("%d: ", __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintD(expr) do{printf("%d: ", __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/

int create_empty_dir(char *parent_directory);
void read_input(char *input_file, int *N, double *u1, double *x_max, double *x_min, int *Roe_first, int *Roe_second, double *CFL, double *delta_x, int *iterations);
void init(double *u, int u1, double delta_x, int N);
void print_array(double *array, int start, int end);
void print_array_to_file(FILE *fp, double *array, int start, int end);
double claculate_f_i_plus_half(double u_i, double u_i_plus_1);
void calculate_vec_f(double *f, double *u, int N);
double calculate_delta_time(double *u, double delta_x, double CFL, int N);
void make_u_physical(double *u_i_plus_half, double u_i, double u_i_plus_1);
void calculate_delta_u(double *f, double *delta_u, double delta_time, double delta_x, int N);
void update_u(double *u, double *delta_u, int N);
double calculate_norm(double *delta_u, int start, int end);
void output(char *output_dir, int N, double u1, double x_max, double x_min, int Roe_first, int Roe_second, double CFL, double delta_x);

int main(int argc, char const *argv[])
{
/* declerations */
    char input_file[MAXDIR], output_dir[MAXDIR], temp_word[MAXWORD];
    double *u, *f, *delta_u, u1, x_max, x_min, CFL, delta_x, delta_time, first_delta_u_norm, current_delta_u_norm;
    int N, Roe_first, Roe_second, iterations;
    FILE *output_u_file, *output_iter_file;  

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

/* creating output directory */
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

/*------------------------------------------------------------*/
/* reading the input */
    read_input(input_file, &N, &u1, &x_max, &x_min, &Roe_first, &Roe_second, &CFL, &delta_x, &iterations);

/* Checking the input */
    dprintINT(N);
    dprintD(u1);
    dprintD(x_max);
    dprintD(x_min);
    dprintD(delta_x);
    dprintINT(Roe_first);
    dprintINT(Roe_second);
    dprintD(CFL);
    dprintINT(iterations);
    printf("--------------------\n");

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
    print_array_to_file(output_u_file, u, 1, N);
/*------------------------------------------------------------*/
/* the loop */

    for (int iter = 0; iter < iterations; iter++) {
        calculate_vec_f(f, u, N);
        delta_time = calculate_delta_time(u, delta_x, CFL, N);
        calculate_delta_u(f, delta_u, delta_time, delta_x, N);
        if (iter == 0) {
            first_delta_u_norm = calculate_norm(delta_u, 1, N);
        }
        current_delta_u_norm = calculate_norm(delta_u, 1, N);

        print_array(delta_u, 1, N);

        update_u(u, delta_u, N);
        print_array_to_file(output_u_file, u, 1, N);
        fprintf(output_iter_file, "%d, %g\n", iter, current_delta_u_norm);

        printf("%d: %g\n", iter, current_delta_u_norm);

        if (current_delta_u_norm/first_delta_u_norm < 1e-6 || isinf(current_delta_u_norm) || isnan(current_delta_u_norm)) {
            break;
        }
    }

    // print_array_to_file(output_u_file, u, 1, N);
    // fprintf(output_iter_file, "%d, %g\n", iteration, current_L2Norm);


/*------------------------------------------------------------*/
/* output */
    output(output_dir, N, u1, x_max, x_min, Roe_first, Roe_second, CFL, delta_x);

/*------------------------------------------------------------*/
/* freeing the memory */
    free(u); 

    fclose(output_u_file);

    return 0;
}

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

/* read input parameters from input file
argument list:
input-file - file pointer to input file
N          - int pointer
u1         - double pointer
x_max      - double pointer
x_min      - double pointer
Roe_first  - int pointer
Roe_second - int pointer
CFL        - double pointer
delta_x    - double pointer 
iterations - max number of desierd iterations */
void read_input(char *input_file, int *N, double *u1, double *x_max, double *x_min, int *Roe_first, int *Roe_second, double *CFL, double *delta_x, int *iterations)
{
    char current_word[MAXWORD];
    float temp;
    FILE *fp = fopen(input_file, "rt");
    if (fp == NULL) {
        fprintf(stderr, "%s:%d:[Error] problem opening input file: %s\n", __FILE__, __LINE__, strerror(errno));
        exit(1);
    }

    while(fscanf(fp, "%s", current_word) != EOF) {
        if (!strcmp(current_word, "N")) {
            fscanf(fp, "%d", N);
            *N = *N;
        } else if (!strcmp(current_word, "u1")) {
            fscanf(fp, "%g", &temp);
            *u1 = (double)temp;
        } else if (!strcmp(current_word, "x_max")) {
            fscanf(fp, "%g", &temp);
            *x_max = (double)temp;
        } else if (!strcmp(current_word, "x_min")) {
            fscanf(fp, "%g", &temp);
            *x_min = (double)temp;
        } else if (!strcmp(current_word, "Roe_first")) {
            fscanf(fp, "%d", Roe_first);
        } else if (!strcmp(current_word, "Roe_second")) {
            fscanf(fp, "%d", Roe_second);
        } else if (!strcmp(current_word, "CFL")) {
            fscanf(fp, "%g", &temp);
            *CFL = (double)temp;
        } else if (!strcmp(current_word, "iterations")) {
            fscanf(fp, "%g", &temp);
            *iterations = (int)temp;
        }
    }

    *delta_x = (double)(*x_max-*x_min)/(double)(*N-1);

    fclose(fp);
}

/* initializing the solution vector with ones and the boundry conditions
argument list:
u       - pointer to the solution array
u1      - boundry condition at x_max
delta_x - double value of delta_x
N       - number of grid points */
void init(double *u, int u1, double delta_x, int N)
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

void make_u_physical(double *u_i_plus_half, double u_i, double u_i_plus_1)
{
    double epsilon = fmax(0, (u_i_plus_1 - u_i) / 2); 

    if (*u_i_plus_half < epsilon) {
        *u_i_plus_half = epsilon;
    }
}

/* calculating the flax for a face
argument list:*/
double claculate_f_i_plus_half(double u_i, double u_i_plus_1)
{
    double u_bar_i_plus_half = (u_i + u_i_plus_1) / 2;
    make_u_physical(&u_bar_i_plus_half, u_i, u_i_plus_1);

    double u_bar_i_plus_half_plus = (u_bar_i_plus_half + fabs(u_bar_i_plus_half)) / 2;
    double u_bar_i_plus_half_minus = (u_bar_i_plus_half - fabs(u_bar_i_plus_half)) / 2;
    
    double F_i = u_i * u_i / 2;
    double F_i_plus_1 = u_i_plus_1 * u_i_plus_1 / 2;

    double f_i_plus_half = (F_i + F_i_plus_1) / 2 - 0.5 * (u_bar_i_plus_half_plus - u_bar_i_plus_half_minus) * ( u_i_plus_1 - u_i);

    return f_i_plus_half;
}

/* calculating the flax vector
argument list:
f - pointer to the flax array
u - pointer to the u array
N - number of grid points */
void calculate_vec_f(double *f, double *u, int N)
{
    for (int i = 0; i < N+1; i++) {
        f[i] = claculate_f_i_plus_half(u[i], u[i+1]);
    }
}

/* calculating the delta time according to the CFL number argument list:
*/
double calculate_delta_time(double *u, double delta_x, double CFL, int N)
{
    double delta_time = __DBL_MAX__;
    for (int i = 1; i < N+1; i++) {
        if (delta_time > CFL / u[i] * delta_x) {
            delta_time = CFL / u[i] * delta_x;
        }
    }

    return delta_time;
}

void calculate_delta_u(double *f, double *delta_u, double delta_time, double delta_x, int N)
{
    for (int i = 1; i < N+1; i++) {
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
void output(char *output_dir, int N, double u1, double x_max, double x_min, int Roe_first, int Roe_second, double CFL, double delta_x)
{
    char temp_word[MAXWORD];
    FILE *meta_data_file;

    strcpy(temp_word, output_dir);
    strcat(temp_word, "/mata_data.txt");
    meta_data_file = fopen(temp_word, "wt");

    fprintf(meta_data_file, "%s, %s, %s, %s, %s, %s, %s, %s\n", "N", "u1", "x_max", "x_min", "Roe_first", "Roe_second", "CFL", "delta_x");
    fprintf(meta_data_file, "%d, %g, %g, %g, %d, %d, %g, %g", N, u1, x_max, x_min, Roe_first, Roe_second, CFL, delta_x);


    fclose(meta_data_file);
}
