#ifndef COMUM_H
#define COMUM_H

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

// Estruturas de iterações e referências
typedef struct {
    int nc;
    int n_tr;
    int n_out;
    int n_vort;
    int itc_max;
    int start_mode;
    double beta;
    double b_art;
    double dtau_f;
    double final_time;
    double eps;
    double eps_mass;
} Iteracoes;

typedef struct {
    double tnu;
    double yf_b;
    double yo_oo;
    double ts;
    double tn_too;
} Referencia;

// Variáveis globais declaradas no comum.cu (acessíveis via extern)
extern __managed__ int imax;
extern __managed__ int jmax;
extern __managed__ int c_f;
extern __managed__ double dt;
extern __managed__ double dtau;
extern __managed__ double beta;
extern __managed__ double v_i;
extern __managed__ double tinf;
extern __managed__ double temp_cylinder;
extern __managed__ double concentracao_inicial;
extern __managed__ double s;
extern __managed__ double lf;
extern __managed__ double q;
extern __managed__ double too;
extern __managed__ double tsup;
extern __managed__ double tempo;

extern __managed__ double *dx_c;
extern __managed__ double *dy_c;
extern __managed__ double *dev_x;
extern __managed__ double *dev_y;
extern __managed__ int *dev_flag;

extern Iteracoes iterations;
extern Referencia ref;

// Funções de configuração e malha
void input_parameters();
void alocar_globais();
void desalocar_globais();
void mesh();
void IC(double *dev_um, double *dev_vm, double *dev_u, double *dev_v, double *dev_p, double *dev_t, double *dev_c, int *dev_flag);
void velocidade_centro_celula(double *dev_um, double *dev_vm, double *dev_u, double *dev_v);
double erro(double *dev_rc, double *dev_res_u, double *dev_res_v);
void output(double *dev_um, double *dev_vm, double *dev_u, double *dev_v, double *dev_p, double *dev_t, double *dev_c, int k);

// PROTÓTIPOS CORRIGIDOS DOS SOLVERS
void solve_U(double *dev_um, double *dev_vm, double *dev_um_n, double *dev_um_tau, double *dev_um_n_tau, double *dev_p, double *dev_t, int *dev_flag, double *dev_ru, double *dev_ui, double *dev_res_u);
void solve_V(double *dev_um, double *dev_vm, double *dev_vm_n, double *dev_vm_tau, double *dev_vm_n_tau, double *dev_p, double *dev_t, int *dev_flag, double *dev_rv, double *dev_vi, double *dev_res_v);
void solve_P(double *dev_um_tau, double *dev_vm_tau, double *dev_p, double *dev_pn, int *dev_flag, double *dev_rc);
void solve_Z(double *dev_um_tau, double *dev_vm_tau, double *dev_t, double *dev_t_tau, double *dev_t_n_tau, double *dev_z, double *dev_h, int *dev_flag, double *dev_rz);
void solve_C(double *dev_um_tau, double *dev_vm_tau, double *dev_c, double *dev_c_tau, double *dev_c_n_tau, int *dev_flag, double *dev_rc);

#endif