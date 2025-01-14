#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>

#define MAXDIR 1000
#define MAXWORD MAXDIR
#define PI M_PI
#define dprintSTRING(expr) do{printf("%d: ", __LINE__); printf(#expr " = %s\n", expr);} while(0)   /* macro for easy debuging*/
#define dprintINT(expr) do{printf("%d: ", __LINE__); printf(#expr " = %d\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintF(expr) do{printf("%d: ", __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/
#define dprintD(expr) do{printf("%d: ", __LINE__); printf(#expr " = %g\n", expr);} while(0)     /* macro for easy debuging*/

int create_empty_dir(char *parent_directory);
void read_input(char *input_file, int *N, double *alpha, double *y_0, double *u_0, double *y_N, double *u_N, double *delta_time, double *delta_y, double *mu);
void init(double *u, int u_0, int u_N, int N);
void print_array(double *array, int len);
void check_delta_time(double delta_time, double delta_y, double mu);
void RHS(double *D, double *u, int N);
void LHS(double *A, double *B, double *C, double delta_time, double delta_y, double mu, double alpha, int N);
void BC(double *A, double *C, double *D, double u_0, double u_N, int N);
int tridiag(double *a, double *b, double *c, double *d, double *u, int is, int ie);
double calculate_norm(double *delta_u, int N);
double step(double *A, double *B, double *C, double *D, double *u, double *delta_u, double delta_time, double delta_y, double mu, double alpha, double u_0, double u_N, int N);
void print_array_to_file(FILE *fp, double *array, int len);
void output(char *output_dir, int N, double y_0, double u_0, double y_N, double u_N, double delta_time, double alpha, double delta_y, double mu);

int main(int argc, char const *argv[])
{
/* declerations */
    char input_fill[MAXDIR], output_dir[MAXDIR], temp_word[MAXWORD];
    double *A, *B, *C, *D, *u, *delta_u,
           y_0, u_0, y_N, u_N, delta_time, alpha, delta_y, mu, first_L2Norm, current_L2Norm;
    int N;
    FILE *output_u_file, *output_iter_file;  

/* getting the input file and output file */
    if (--argc != 2 && argc != 4) {
        fprintf(stderr, "%s:%d: [Error] not right usage... Usage: main 'input file' 'output directory'\n", __FILE__, __LINE__);
        return -1;
    }

    strncpy(input_fill, (*(++argv)), MAXDIR);

    if (input_fill[MAXDIR-1] != '\0') {
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
    read_input(input_fill, &N, &alpha, &y_0, &u_0, &y_N, &u_N, &delta_time, &delta_y, &mu);

/* Checking the input */
    dprintINT(N);
    dprintD(y_0);
    dprintD(u_0);
    dprintD(y_N);
    dprintD(u_N);
    dprintD(mu);
    dprintD(alpha);
    dprintD(delta_time);
    dprintD(delta_y);
    printf("--------------------\n");

/*------------------------------------------------------------*/
/* allocating the matrices */

    A = (double *)malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {
        A[i] = 0;
    }
    B = (double *)malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {
        B[i] = 0;
    }
    C = (double *)malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {
        C[i] = 0;
    }
    D = (double *)malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {
        D[i] = 0;
    }
    u = (double *)malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {
        u[i] = 0;
    }
    delta_u = (double *)malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {
        delta_u[i] = 0;
    }

/*------------------------------------------------------------*/
/* initializtion */
    init(u, u_0, u_N, N);

    print_array_to_file(output_u_file, u, N);
    
/*------------------------------------------------------------*/
/* the loop */
    for (int iteration = 0; iteration < 2e5; iteration++) {
        if (iteration == 0) {
            first_L2Norm = step(A, B, C, D, u, delta_u, delta_time, delta_y, mu, alpha, u_0, u_N, N);
            current_L2Norm = first_L2Norm;
        } else {
            current_L2Norm = step(A, B, C, D, u, delta_u, delta_time, delta_y, mu, alpha, u_0, u_N, N);
        }

        printf("%d: %g\n", iteration, current_L2Norm);
        print_array_to_file(output_u_file, u, N);
        fprintf(output_iter_file, "%d, %g\n", iteration, current_L2Norm);

        if (current_L2Norm/first_L2Norm < 1e-6 || isinf(current_L2Norm) || isnan(current_L2Norm)) {
            break;
        }
    }

    print_array(u, N);

/*------------------------------------------------------------*/
/* output */
    output(output_dir, N, y_0, u_0, y_N, u_N, delta_time, alpha, delta_y, mu);

/*------------------------------------------------------------*/
/* freeing the memory */
    free(A);
    free(B);
    free(C);
    free(D);
    free(u); 
    free(delta_u);

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
alpha      - double pointer
y_0        - double pointer
u_0        - double pointer
y_N        - double pointer
u_N        - double pointer
delta_time - double pointer
delta_y    - double pointer
mu         - double pointer */
void read_input(char *input_file, int *N, double *alpha, double *y_0, double *u_0, double *y_N, double *u_N, double *delta_time, double *delta_y, double *mu)
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
        } else if (!strcmp(current_word, "y_0")) {
            fscanf(fp, "%g", &temp);
            *y_0 = (double)temp;
        } else if (!strcmp(current_word, "u_0")) {
            fscanf(fp, "%g", &temp);
            *u_0 = (double)temp;
        } else if (!strcmp(current_word, "y_N")) {
            fscanf(fp, "%g", &temp);
            *y_N = (double)temp;
        } else if (!strcmp(current_word, "u_N")) {
            fscanf(fp, "%g", &temp);
            *u_N = (double)temp;
        } else if (!strcmp(current_word, "mu")) {
            fscanf(fp, "%g", &temp);
            *mu = (double)temp;
        } else if (!strcmp(current_word, "alpha")) {
            fscanf(fp, "%g", &temp);
            *alpha = (double)temp;
        } else if (!strcmp(current_word, "delta_time")) {
            fscanf(fp, "%g", &temp);
            *delta_time = (double)temp;
        }
    }
    *delta_y = fabs((*y_0 - *y_N) / (*N-1));

    if (*alpha <= 0.5) {
        check_delta_time(*delta_time, *delta_y, *mu);
    }

    fclose(fp);
}

/* initializing the solution vector with ones and the boundry conditions
argument list:
u   - pointer to the solution array
u_0 - boundry condition at zero
u_N - boundry condition at N
N   - number of grid points */
void init(double *u, int u_0, int u_N, int N)
{
    for (int i = 1; i < N-1; i++) {
        u[i] = 1;
    }
    u[0] = u_0;
    u[N-1] = u_N;
}

/* printing a double array with length 'len' to 'stdin'
argument list:
array - pointer to the array
len   - number of elements in the array */
void print_array(double *array, int len)
{
    printf("--------------------\n");
    for (int i = 0; i < len; i++) {
        printf("%g ", array[i]);
    }
    printf("\n");
    printf("--------------------\n");
}

/* checking if the inputed delta t will result in converging solution
arugemnt list: 
delta_time - time step between iterations
delta_y    - the distance between grid points
mu         - viscosity */
void check_delta_time(double delta_time, double delta_y, double mu)
{
    if (delta_time >= delta_y * delta_y / 2 / mu) {
        printf("%s:%d:[Warning] delta_time too big\n", __FILE__, __LINE__);
    }
}

/* calculate the RHS (vector D) 
argumet list:
D - pointer to the array
u - pointer to the flow solution
N - number of grid points */
void RHS(double *D, double *u, int N)
{
    for (int i = 1; i < N-1; i++) {
        D[i] = u[i+1] - 2*u[i] + u[i-1];
    }
}

/* calculate the LHS (vectors A, B, C)
A          - pointer to the A array
B          - pointer to the B array
C          - pointer to the C array 
delta_time - time step between iterations
delta_y    - distance between to grid points
mu         - viscosity
alpha      - parameter that control the scheme
N          - number of grid points */
void LHS(double *A, double *B, double *C, double delta_time, double delta_y, double mu, double alpha, int N)
{
    for (int i = 1; i < N-1; i++) {
        A[i] = -alpha;
        B[i] = delta_y * delta_y / mu / delta_time + 2*alpha;
        C[i] = -alpha;
    }
}

/* setting the boundry conditions
A   - pointer to the A array
B   - pointer to the B array
C   - pointer to the C array 
D   - pointer to the D array 
u_0 - boundry condition at zero
u_N - boundry condition at N
N   - number of grid points */
void BC(double *A, double *C, double *D, double u_0, double u_N, int N)
{
    D[1] = D[1] - A[1]*u_0;
    A[1] = 0;

    D[N-1] = D[N-1] - C[N-1]*u_N;
    C[N-1] = 0;
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

/* calculating the second norma of the vector 'delta_u'
argument list:
delta_u - pointer to the vector elements array
N       - number of gird points */
double calculate_norm(double *delta_u, int N)
{
    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += delta_u[i] * delta_u[i];
    }
    return sqrt(sum);
}

/* preforming the step of the shceme 
argument list:
A          - pointer to the A array
B          - pointer to the B array
C          - pointer to the C array 
D          - pointer to the D array 
u          - pointer to the u array 
delta_u    - pointer to the delta_u array 
delta_time - time step between iterations
delta_y    - distance between to grid points
mu         - viscosity
alpha      - parameter that control the scheme
u_0        - boundry condition at zero
u_N        - boundry condition at N
N          - number of grid points */
double step(double *A, double *B, double *C, double *D, double *u, double *delta_u, double delta_time, double delta_y, double mu, double alpha, double u_0, double u_N, int N)
{
    /* zero A, B, C, D*/
    for (int i = 0; i < N; i++) {
        A[i] = 0; 
        B[i] = 0; 
        C[i] = 0; 
        D[i] = 0; 
    }

    RHS(D, u, N);
    if (alpha == 0) {
        for (int i = 0; i < N; i++) {
            delta_u[i] = mu * delta_time / delta_y / delta_y * D[i];
        }
        for (int i = 0; i < N; i++) {
            u[i] = u[i] + delta_u[i];
        }
        return calculate_norm(delta_u, N);
    }

    LHS(A, B, C, delta_time, delta_y, mu, alpha, N);
    BC(A, C, D, u_0, u_N, N);

    double norm = calculate_norm(D, N);

    tridiag(A, B, C, D, delta_u, 1, N-2);


    for (int i = 0; i < N; i++) {
        u[i] = u[i] + delta_u[i];
    }

    return norm;
}

/* printing a double array with length 'len' to file 'fp' 
argument list:
fp    - file pointer
array - pointer to the array
len   - number of elements in the array */
void print_array_to_file(FILE *fp, double *array, int len)
{
    for (int i = 0; i < len; i++) {
        fprintf(fp, "%g ", array[i]);
    }
    fprintf(fp, "\n");
}

/* ouputing metadata of the solution 
argumetn list:
output_dir - name of the output directory
N          - number of grid points
y_0, u_0   - boundry condtions at zero
y_N, u_N   - boundry condtions at N
delta_time - time step between iterations
alpha      - parameter that control the scheme
delta_y    - distance between to grid points
mu         - viscosity */
void output(char *output_dir, int N, double y_0, double u_0, double y_N, double u_N, double delta_time, double alpha, double delta_y, double mu)
{
    char temp_word[MAXWORD];
    FILE *meta_data_file;

    strcpy(temp_word, output_dir);
    strcat(temp_word, "/mata_data.txt");
    meta_data_file = fopen(temp_word, "wt");

    fprintf(meta_data_file, "%s, %s, %s, %s\n", "N", "y_0", "u_0", "y_N");
    fprintf(meta_data_file, "%d, %g, %g, %g", N, y_0, u_0, y_N);


    fclose(meta_data_file);
}
