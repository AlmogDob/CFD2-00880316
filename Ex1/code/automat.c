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

int create_empty_dir(char *parent_directory);
void print_command_to_file(FILE *fp, char *program, ...);
void create_input_file(char *file_name, int i_TEL, int i_LE, int i_TEU, 
                       int j_TEL, int j_LE, int j_TEU, double Mach_inf, 
                       double angle_of_attack_deg, double density, 
                       double environment_pressure, double delta_t,
                       double Gamma, double epse, double max_iteration);
void create_mesh_input_file(char *file_name, double t, int i_max,
                            int j_max, int i_TEL, int i_LE, int i_TEU,
                            double delta_y, double XSF, double YSF, 
                            double x_int, double r, double omega,
                            int phi, int psi);

int main()
{
    char parent_dir[] = "./auto_results";
    if (create_empty_dir("./auto_results") != 0) {
        return 1;
    }

    char temp_dir[BUFSIZ], temp1[BUFSIZ], temp_input[BUFSIZ], temp_num[BUFSIZ];

    strncpy(temp_dir, parent_dir, BUFSIZ);
    strncat(temp_dir, "/command_to_run.txt", BUFSIZ/2);
    FILE *fp = fopen(temp_dir, "wt");

    fprintf(fp, "make build_solver\n");

    double times[] = {1.00000000000000e-05,
                      1.35387618002254e-05,
                      1.83298071083244e-05,
                      2.48162892283682e-05,
                      3.35981828628378e-05,
                      4.54877794700378e-05,
                      6.15848211066027e-05,
                      8.33782223471788e-05,
                      0.000112883789168469,
                      0.000152830673265877,
                      0.000206913808111479,
                      0.000280135676119887,
                      0.000379269019073225,
                      0.000513483290743755,
                      0.000695192796177561,
                      0.000941204967268067,
                      0.00127427498570313,
                      0.00172521054994204,
                      0.00233572146909012,
                      0.00316227766016838};

    for (int i = 1; i <= 2; i++) {
        
        // strncpy(temp_dir, parent_dir, BUFSIZ);
        // strncat(temp_dir, "/mesh_input", BUFSIZ/2);
        // sprintf(temp1, "%d.txt", i);
        // strncat(temp_dir, temp1, BUFSIZ/2);
        // create_mesh_input_file(temp_dir, 
        //                        0.12,            /* t       */
        //                        50+4*(i-1),      /* i_max   */
        //                        25+2*(i-1),      /* j_max   */
        //                        11+4*(i-1),      /* i_TEL   */
        //                        (50+4*(i-1))/2,  /* i_LE    */
        //                        39+4*(i-1)*2,    /* i_TEU   */
        //                        0.02,            /* delta_y */
        //                        1.15,            /* XSF     */
        //                        1.15,            /* YSF     */
        //                        1.008930411365,  /* x_int   */
        //                        0.001,           /* r       */
        //                        1,               /* omega   */
        //                        -1,              /* phi     */
        //                        -1);             /* psi     */

        strncpy(temp_dir, parent_dir, BUFSIZ);
        strncat(temp_dir, "/input", BUFSIZ/2);
        sprintf(temp1, "%d.txt", i);
        strncat(temp_dir, temp1, BUFSIZ/2);
        create_input_file(temp_dir,
                          11,           /* i_TEL                */
                          25,           /* i_LE                 */
                          39,           /* i_TEU                */
                          0,            /* j_TEL                */
                          0,            /* j_LE                 */
                          0,            /* j_TEU                */
                          0.9,          /* Mach_inf             */
                          0,            /* angle_of_attack_deg  */
                          1.225,        /* density              */
                          101325,       /* environment_pressure */
                          times[i-1],   /* delta_t              */
                          1.4,          /* Gamma                */
                          0.06,         /* epse                 */
                          1e6);         /* max_iteration        */


        // strncpy(temp_input, parent_dir, BUFSIZ);
        // strncat(temp_input, "/mesh_input", BUFSIZ/2);
        // sprintf(temp1, "%d.txt", i);
        // strncat(temp_input, temp1, BUFSIZ/2);

        // print_command_to_file(fp,
        //                     "./mesh_generate",
        //                     temp_input,
        //                     "./mesh_output.txt",
        //                     "./auto_results",
        //                     temp_num,
        //                     NULL);


        strncpy(temp_input, parent_dir, BUFSIZ);
        strncat(temp_input, "/input", BUFSIZ/2);
        sprintf(temp1, "%d.txt", i);
        strncat(temp_input, temp1, BUFSIZ/2);

        sprintf(temp_num, "%d", i);

        print_command_to_file(fp,
                            "./solver",
                            temp_input,
                            "./mesh_output.txt",
                            "./auto_results",
                            temp_num,
                            NULL);

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

void create_input_file(char *file_name, int i_TEL, int i_LE, int i_TEU, 
                       int j_TEL, int j_LE, int j_TEU, double Mach_inf, 
                       double angle_of_attack_deg, double density, 
                       double environment_pressure, double delta_t,
                       double Gamma, double epse, double max_iteration)
{
    FILE *input_fp = fopen(file_name, "wt");

    dfprintINT(input_fp, i_TEL);
    dfprintINT(input_fp, i_LE);
    dfprintINT(input_fp, i_TEU);
    dfprintINT(input_fp, j_TEL);
    dfprintINT(input_fp, j_LE);
    dfprintINT(input_fp, j_TEU);
    dfprintD(input_fp, Mach_inf);
    dfprintD(input_fp, angle_of_attack_deg);
    dfprintD(input_fp, density);
    dfprintD(input_fp, environment_pressure);
    dfprintD(input_fp, delta_t);
    dfprintD(input_fp, Gamma);
    dfprintD(input_fp, epse);
    dfprintD(input_fp, max_iteration);

    fclose(input_fp);
}
