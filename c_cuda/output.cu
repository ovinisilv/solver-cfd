#include "comum.h"
#include <sys/stat.h>
#include <sys/types.h>
#define idx i*(jmax+1)+j

static void ensure_directory(const char *path){
    struct stat st;
    if(stat(path, &st) != 0){
        mkdir(path, 0777);
    }
}

static FILE *open_output_file(const char *path, const char *mode){
    FILE *arquivo = fopen(path, mode);
    if(!arquivo){
        fprintf(stderr, "ERROR: Could not open '%s' for writing\n", path);
    }
    return arquivo;
}

void output(double *dev_um, double *dev_vm, double *dev_u, double *dev_v, double *dev_p, double *dev_t, double *dev_c, int k){
    FILE *arquivo;
    
    ensure_directory("data");
    ensure_directory("data/restart");

    arquivo = open_output_file("data/output_variables.dat", "w");
    if(arquivo){
        for(int i = 1; i <= imax; i++){
            for(int j = 1; j <= jmax; j++)
                fprintf(arquivo, "%lf %lf %lf %lf %lf %lf %lf\n", dev_x[i], dev_y[j], dev_u[idx], dev_v[idx], dev_p[idx], dev_t[idx], dev_c[idx]);
        }
        fflush(arquivo);
        fclose(arquivo);
    }

    //--- RESTART/RESTART.dat ---
    arquivo = fopen("data/restart/restartU.dat", "w");
    for(int i = 1; i <= (imax+1); i++){
        for(int j = 1; j <= jmax; j++)
            fprintf(arquivo, "%lf\n", dev_um[idx]);
    }
    fclose(arquivo);

    arquivo = fopen("data/restart/restartV.dat","w");
    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= (jmax+1); j++)
            fprintf(arquivo, "%lf\n", dev_vm[i*(jmax+2)+j]);
    }
    fclose(arquivo);

    arquivo = fopen("data/restart/restartPTC.dat", "w");
    for(int i = 1; i <= imax; i++){
        for(int j = 1; j <= jmax; j++)
            fprintf(arquivo, "%lf %lf %lf\n", dev_p[idx], dev_t[idx], dev_c[idx]);       
    }
    fclose(arquivo);
}