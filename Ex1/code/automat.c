#include <stdio.h>
#include <sys/stat.h>
#include <errno.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdarg.h>

#define dfprintINT(fp, expr) do{fprintf(fp, #expr "\n%d\n\n", expr);} while(0)     /* macro for easy debuging*/
#define dfprintD(fp, expr) do{fprintf(fp, #expr "\n%g\n\n", expr);} while(0)     /* macro for easy debuging*/
#define dfprintS(fp, expr) do{fprintf(fp, #expr "\n%s\n\n", expr);} while(0)     /* macro for easy debuging*/

int create_empty_dir(char *parent_directory);
void print_command_to_file(FILE *fp, char *program, ...);
void create_input_file(char *file_name, int N, double u0, double u1, double x_max, double x_min, double delta_x, double delta_time, double k, double b, double c, double mu, char *method, char *limiter, double CFL, double w, double theta, int iterations, double final_time);

int main()
{
    char temp_dir[BUFSIZ], temp1[BUFSIZ], temp_input[BUFSIZ];
    int output_counter = 0;

    char parent_dir[] = "./auto";
    if (create_empty_dir(parent_dir) != 0) {
        return 1;
    }

    strncpy(temp_dir, parent_dir, BUFSIZ);
    strncat(temp_dir, "/command_to_run.txt", BUFSIZ/2);
    FILE *fp = fopen(temp_dir, "wt");

    fprintf(fp, "make build_solver\n");

    char *methods[] = {"Beam_and_Warming", "MacCormack", "Roe_first"};
    double ws[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
    double thetas[] = {0.5, 1};
    double mus[] = {0.25, 0.001};
    double delta_times[] = {1, 0.75, 0.5};

    for (int index5 = 0; index5 < 3; index5++) {
        for (int index4 = 0; index4 < 11; index4++) {
            for (int index3 = 0; index3 < 2; index3++) {
                for (int index2 = 0; index2 < 2; index2++) {
                    for (int index1 = 0; index1 < 3; index1++) {
                        int N = 41;
                        double u0 = 0;
                        double u1 = 1;
                        double x_max = 40;
                        double x_min = 0;
                        double delta_x = 1;
                        double delta_time = delta_times[index1];
                        double k = -1;
                        double b = -1;
                        double c = 0.5;
                        double mu = mus[index2];
                        char *method = methods[index5];
                        char *limiter = 0;
                        double CFL = 0;
                        double w = ws[index4];
                        double theta = thetas[index3];
                        int iterations = 1e4;
                        double final_time = 0;


                        strncpy(temp_input, parent_dir, BUFSIZ);
                        strncat(temp_input, "/input", BUFSIZ/2);
                        sprintf(temp1, "%d.txt", output_counter++);
                        strncat(temp_input, temp1, BUFSIZ/2);
                        create_input_file(temp_input, N, u0, u1, x_max, x_min, delta_x, delta_time, k, b, c, mu, method, limiter, CFL, w, theta, iterations, final_time);

                        strncpy(temp_dir, parent_dir, BUFSIZ);
                        strncat(temp_dir, "/results", BUFSIZ/2);
                        
                        print_command_to_file(fp,
                                            "./solver",
                                            temp_input,
                                            temp_dir,
                                            NULL);

                    }
                }
            }
        }
    }

    fprintf(fp, "make clean_solver\n");

    return 0;
}

/* if allready exisest, delet all the files inside 
returns 0 on success
this functin handls the errors so on fail just quit */
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
            printf("\n");
            while (entity != NULL) {   /* directory is not empty */
                strncpy(path_to_remove, parent_directory, BUFSIZ);
                strncat(path_to_remove, "/", BUFSIZ/2);
                strncat(path_to_remove, entity->d_name, BUFSIZ);
                printf("%hhd: %s\n", entity->d_type, path_to_remove);
                if (entity->d_type == DT_REG) {
                    if (remove(path_to_remove) != 0) {
                        fprintf(stderr, "%s:%d: [Error] problem removing '%s': %s\n", __FILE__, __LINE__, path_to_remove, strerror(errno));
                        return 1;
                    }
                    printf("remove %s\n", path_to_remove);
                }
                entity = readdir(dir);
            }


            printf("\ndirectory already exist\n\n");

            closedir(dir);

            return 0;
        }

        fprintf(stderr, "%s:%d: [Error] problem making '%s': %s\n", __FILE__, __LINE__, parent_directory, strerror(errno));
        return 1;
    }
    return 0;
}

void print_command_to_file(FILE *fp, char *program, ...)
{
    fprintf(fp, "%s ", program);
    va_list args;
    va_start(args, program);
    char *temp_string;
    while ((temp_string = va_arg(args, char *)) != NULL) {
        fprintf(fp, "%s ", temp_string);
    }
    va_end(args);
    fprintf(fp, "\n");
}

void create_input_file(char *file_name, int N, double u0, double u1, double x_max, double x_min, double delta_x, double delta_time, double k, double b, double c, double mu, char *method, char *limiter, double CFL, double w, double theta, int iterations, double final_time)
{
    FILE *input_fp = fopen(file_name, "wt");

    dfprintINT(input_fp, N);
    dfprintD(input_fp, u0);
    dfprintD(input_fp, u1);
    dfprintD(input_fp, x_max);
    dfprintD(input_fp, x_min);
    dfprintD(input_fp, delta_x);
    dfprintD(input_fp, delta_time);
    dfprintD(input_fp, k);
    dfprintD(input_fp, b);
    dfprintD(input_fp, c);
    dfprintD(input_fp, mu);
    dfprintS(input_fp, method);
    dfprintS(input_fp, limiter);
    dfprintD(input_fp, CFL);
    dfprintD(input_fp, w);
    dfprintD(input_fp, theta);
    dfprintINT(input_fp, iterations);
    dfprintD(input_fp, final_time);

    fclose(input_fp);
}
