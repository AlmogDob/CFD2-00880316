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

#ifdef __linux__
    #include <sys/stat.h>
    #include <dirent.h>
    #define ON_LINUX 1
    int create_empty_dir(char *parent_directory);
#endif

void read_input(char *input_file, int *N, double *u0, double *u1, double *x_max, double *x_min, double *k, double *b, double *c, double *mu, double *CFL, double *w, double *theta, double *delta_x, double *delta_time, int *iterations, double *final_time);

int main(int argc, char const *argv[])
{
/* declerations */
    char input_file[MAXDIR], output_dir[MAXDIR], temp_word[MAXWORD];

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

/* Checking the input */
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

/*------------------------------------------------------------*/
/* initializtion */

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
u0         - double pointer
u1         - double pointer
x_max      - double pointer
x_min      - double pointer
k          - double pointer
b          - double pointer
c          - double pointer
mu         - double pointer
CFL        - double pointer
w          - double pointer
theta      - double pointer
delta_x    - double pointer 
delta_time - double pointer 
iterations - int pointer max number of desierd iterations 
final_time - int pointer to the final time of the program */
void read_input(char *input_file, int *N, double *u0, double *u1, double *x_max, double *x_min, double *k, double *b, double *c, double *mu, double *CFL, double *w, double *theta, double *delta_x, double *delta_time, int *iterations, double *final_time)
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
        }
    }

    fclose(fp);
}
