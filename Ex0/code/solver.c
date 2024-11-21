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
void read_input(char *input_file, int *N, int *alpha, double *y_0, double *u_0, double *y_1, double *u_1, double *delta_time);
void init(double *u);

int main(int argc, char const *argv[])
{
/* declerations */
    char input_fill[MAXDIR], output_dir[MAXDIR];

    double *A, *B, *C, *D, *u,
           y_0, u_0, y_1, u_1, delta_time;

    int N, alpha;

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
    read_input(input_fill, &N, &alpha, &y_0, &u_0, &y_1, &u_1, &delta_time);

/* Checking the input */
    printf("--------------------\n");
    dprintINT(N);
    dprintINT(alpha);
    dprintD(y_0);
    dprintD(u_0);
    dprintD(y_1);
    dprintD(u_1);
    dprintD(delta_time);
    printf("--------------------\n");

/*------------------------------------------------------------*/
/* allocating the matrices */

    A = (double *)malloc(sizeof(double) * N+1);
    for (int i = 0; i < N+1; i++) {
        A[i] = 0;
    }
    B = (double *)malloc(sizeof(double) * N+1);
    for (int i = 0; i < N+1; i++) {
        B[i] = 0;
    }
    C = (double *)malloc(sizeof(double) * N+1);
    for (int i = 0; i < N+1; i++) {
        C[i] = 0;
    }
    D = (double *)malloc(sizeof(double) * N+1);
    for (int i = 0; i < N+1; i++) {
        D[i] = 0;
    }
    u = (double *)malloc(sizeof(double) * N+1);
    for (int i = 0; i < N+1; i++) {
        u[i] = 0;
    }

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
void read_input(char *input_file, int *N, int *alpha, double *y_0, double *u_0, double *y_1, double *u_1, double *delta_time)
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
        } else if (!strcmp(current_word, "alpha")) {
            fscanf(fp, "%d", alpha);
        } else if (!strcmp(current_word, "y_0")) {
            fscanf(fp, "%g", &temp);
            *y_0 = (double)temp;
        } else if (!strcmp(current_word, "u_0")) {
            fscanf(fp, "%g", &temp);
            *u_0 = (double)temp;
        } else if (!strcmp(current_word, "y_1")) {
            fscanf(fp, "%g", &temp);
            *y_1 = (double)temp;
        } else if (!strcmp(current_word, "u_1")) {
            fscanf(fp, "%g", &temp);
            *u_1 = (double)temp;
        } else if (!strcmp(current_word, "delta_time")) {
            fscanf(fp, "%g", &temp);
            *delta_time = (double)temp;
        }
    }
}

void init(double *u)
{

}
