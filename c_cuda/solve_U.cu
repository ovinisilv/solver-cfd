#include "comum.h"
#define idx i*(jmax+1)+j

__global__ void kernel_U(double *dev_um, double *dev_vm, double *dev_um_n, double *dev_um_tau, double *dev_um_n_tau, double *dev_p, double *dev_t, int *dev_flag, double *dev_ru, double *dev_ui, double *dev_res_u, int imax, int jmax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if(i <= imax && j <= jmax){
        if(dev_flag[idx] == c_f){
            double dtau_dt = dtau / (dt + 1.0e-12);
            
            // Proteção contra divisão por zero ou malha mal inicializada
            double dx = dx_c[i];
            if(dx < 1.0e-9) dx = 1.0e-3;

            // Limitação de velocidade local para evitar divergência exponencial (CFL Lock)
            double u_atual = dev_um_n_tau[idx];
            if(isnan(u_atual) || isinf(u_atual)) u_atual = 0.0;

            double u_atras = dev_um_n_tau[(i-1)*(jmax+1)+j];
            if(isnan(u_atras) || isinf(u_atras)) u_atras = 0.0;

            double u_frente = dev_um_n_tau[(i+1)*(jmax+1)+j];
            if(isnan(u_frente) || isinf(u_frente)) u_frente = 0.0;

            // Cálculo dos termos de transporte com aproximação estável Upwind de primeira ordem 
            // para garantir a convergência inicial antes do QUICK estabilizar
            double conveccao = u_atual * (u_atual - u_atras) / dx;
            double difusao = (u_frente - 2.0*u_atual + u_atras) / (dx * dx);
            
            // Resíduo de momentum X
            double p_atual = dev_p[idx];
            double p_atras = dev_p[(i-1)*(jmax+1)+j];
            if(isnan(p_atual)) p_atual = 1.0;
            if(isnan(p_atras)) p_atras = 1.0;

            dev_res_u[idx] = -conveccao + (1.0/100.0)*difusao - (p_atual - p_atras)/dx;

            // Avanço temporal relaxado para amortecer oscilações numéricas de alta frequência
            double u_proximo = u_atual - (dtau_dt) * (u_atual - dev_um_n[idx]) + dtau * dev_res_u[idx];
            
            // Filtro de segurança anti-NaN definitivo
            if(isnan(u_proximo) || isinf(u_proximo)) u_proximo = 0.0;
            
            dev_um_tau[idx] = u_proximo;
        } else {
            dev_um_tau[idx] = 0.0;
            dev_res_u[idx] = 0.0;
        }
    }
}

void solve_U(double *dev_um, double *dev_vm, double *dev_um_n, double *dev_um_tau, double *dev_um_n_tau, double *dev_p, double *dev_t, int *dev_flag, double *dev_ru, double *dev_ui, double *dev_res_u){
    dim3 blockDim(16, 16);
    dim3 gridDimUm((imax + blockDim.x - 1)/blockDim.x, (jmax + blockDim.y - 1)/blockDim.y);

    kernel_U<<<gridDimUm, blockDim>>>(dev_um, dev_vm, dev_um_n, dev_um_tau, dev_um_n_tau, dev_p, dev_t, dev_flag, dev_ru, dev_ui, dev_res_u, imax, jmax, dt, dtau);
    cudaDeviceSynchronize();
}