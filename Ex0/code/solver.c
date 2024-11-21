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
void read_input(char *input_file, int *N, int *len, double *alpha, double *y_0, double *u_0, double *y_N, double *u_N, double *delta_time, double *delta_y, double *mu);
void init(double *u, int len);
void print_array(double *array, int len);
void check_delta_time(double delta_time, double delta_y, double mu);
void RHS(double *D, double *u, int N);
void LHS(double *A, double *B, double *C, double delta_time, double delta_y, double mu, double alpha, int N);
void BC(double *A, double *C, double *D, double u_0, double u_1, int N);

int main(int argc, char const *argv[])
{
/* declerations */
    char input_fill[MAXDIR], output_dir[MAXDIR];

    double *A, *B, *C, *D, *u,
           y_0, u_0, y_N, u_N, delta_time, alpha, delta_y, mu;

    int N, len;

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

/*------------------------------------------------------------*/
/* reading the input */
    read_input(input_fill, &N, &len, &alpha, &y_0, &u_0, &y_N, &u_N, &delta_time, &delta_y, &mu);

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

    A = (double *)malloc(sizeof(double) * len);
    for (int i = 0; i < len; i++) {
        A[i] = 0;
    }
    B = (double *)malloc(sizeof(double) * len);
    for (int i = 0; i < len; i++) {
        B[i] = 0;
    }
    C = (double *)malloc(sizeof(double) * len);
    for (int i = 0; i < len; i++) {
        C[i] = 0;
    }
    D = (double *)malloc(sizeof(double) * len);
    for (int i = 0; i < len; i++) {
        D[i] = 0;
    }
    u = (double *)malloc(sizeof(double) * len);
    for (int i = 0; i < len; i++) {
        u[i] = 0;
    }

/*------------------------------------------------------------*/
/* initializtion */
    init(u, len);

    print_array(u, len);

/*------------------------------------------------------------*/
/* the loop */

/*------------------------------------------------------------*/
/* output */

/*------------------------------------------------------------*/
/* freeing the memory */

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
fp - file pointer to input file */
void read_input(char *input_file, int *N, int *len, double *alpha, double *y_0, double *u_0, double *y_N, double *u_N, double *delta_time, double *delta_y, double *mu)
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
            *len = *N+1;
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

void init(double *u, int len)
{
    for (int i = 0; i < len; i++) {
        u[i] = 1;
    }
}

void print_array(double *array, int len)
{
    printf("--------------------\n");
    for (int i = 0; i < len; i++) {
        printf("%g ", array[i]);
    }
    printf("\n");
    printf("--------------------\n");
}

void check_delta_time(double delta_time, double delta_y, double mu)
{
    if (delta_time >= delta_y * delta_y / 2 / mu) {
        printf("%s:%d:[Warning] delta_time too big\n", __FILE__, __LINE__);
    }
}

void RHS(double *D, double *u, int N)
{
    for (int i = 1; i < N; i++) {
        D[i] = u[i+1] - 2*u[i] + u[i-1];
    }
}

void LHS(double *A, double *B, double *C, double delta_time, double delta_y, double mu, double alpha, int N)
{
    for (int i = 1; i < N; i++) {
        A[i] = -alpha;
        B[i] = delta_y * delta_y / mu / delta_time + 2*alpha;
        C[i] = -alpha;
    }
}

void BC(double *A, double *C, double *D, double u_0, double u_N, int N)
{
    D[1] = D[1] - A[1]*u_0;
    A[1] = 0;

    D[N-1] = D[N-1] - C[N-1]*u_N;
}
