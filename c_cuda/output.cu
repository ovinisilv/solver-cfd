#include "comum.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <stdarg.h>
#include <limits.h>
#ifdef _WIN32
#include <direct.h>
#define getcwd _getcwd
#else
#include <unistd.h>
#endif
#define idx i*(jmax+1)+j

static void ensure_directory(const char *path){
    struct stat st;
    if(stat(path, &st) != 0){
        mkdir(path, 0777);
    }
}

static FILE *open_log_file(void){
    FILE *log = fopen("output_log.txt", "a");
    return log;
}

static void log_message(const char *fmt, ...){
    FILE *log = open_log_file();
    if(!log) return;
    va_list args;
    va_start(args, fmt);
    vfprintf(log, fmt, args);
    va_end(args);
    fprintf(log, "\n");
    fflush(log);
    fclose(log);
}

static FILE *open_output_file(const char *path, const char *mode){
    FILE *arquivo = fopen(path, mode);
    if(!arquivo){
        log_message("ERROR: Could not open '%s' for writing", path);
    }
    return arquivo;
}

void output(double *dev_um, double *dev_vm, double *dev_u, double *dev_v, double *dev_p, double *dev_t, double *dev_c, int k){
    FILE *arquivo;
    char cwd[PATH_MAX];

    if(getcwd(cwd, sizeof(cwd))){
        log_message("output() cwd=%s", cwd);
    } else {
        log_message("output() cwd=UNKNOWN");
    }
    log_message("output() called with imax=%d jmax=%d", imax, jmax);
    
    ensure_directory("data");
    ensure_directory("data/restart");

    arquivo = open_output_file("data/output_variables.dat", "w");
    if(arquivo){
        for(int i = 1; i <= imax; i++){
            for(int j = 1; j <= jmax; j++)
                fprintf(arquivo, "%lf %lf %lf %lf %lf %lf %lf\n", dev_x[i], dev_y[j], dev_u[idx], dev_v[idx], dev_p[idx], dev_t[idx], dev_c[idx]);
        }
        fflush(arquivo);
        if(fclose(arquivo) != 0){
            log_message("ERROR: fclose failed for data/output_variables.dat");
        } else {
            log_message("WROTE: data/output_variables.dat successfully");
        }
    } else {
        log_message("ERROR: Could not open data/output_variables.dat");
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