#include "comum.h"
#define idx i*(jmax+1)+j

__global__ void kernel_U(double *dev_um, double *dev_vm, double *dev_um_n, double *dev_um_tau, double *dev_um_n_tau, double *dev_p, double *dev_t, int *dev_flag, double *dev_ru, double *dev_ui, double *dev_res_u, int imax, int jmax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if(i <= imax && j <= jmax){
        if(dev_flag[idx] == c_f){
            double dtau_dt = dtau / dt;

            // Coeficientes e diferenciação espacial do escoamento (QUICK / Central)
            double u_atual = dev_um_n_tau[idx];
            double u_atras = dev_um_n_tau[(i-1)*(jmax+1)+j];
            double u_frente = dev_um_n_tau[(i+1)*(jmax+1)+j];
            double u_baixo = dev_um_n_tau[i*(jmax+1)+j-1];
            double u_cima = dev_um_n_tau[i*(jmax+1)+j+1];

            double de_dx = (u_frente - u_atras) / (2.0 * dev_dx[i]);
            double de_dy = (u_cima - u_baixo) / (2.0 * dev_dy[j]);
            double d2u_dx2 = (u_frente - 2.0*u_atual + u_atras) / (dev_dx[i]*dev_dx[i]);
            double d2u_dy2 = (u_cima - 2.0*u_atual + u_baixo) / (dev_dy[j]*dev_dy[j]);

            // Viscosidade cinemática efetiva aproximada obtida dos parâmetros de referência (Re ~ 100)
            double nu = 0.01; 

            dev_res_u[idx] = -(u_atual * de_dx) + nu * (d2u_dx2 + d2u_dy2) - (dev_p[idx] - dev_p[(i-1)*(jmax+1)+j]) / dev_dx[i];

            // Atualização explícita padrão do passo pseudo-temporal transiente
            dev_um_tau[idx] = dev_um_n_tau[idx] - dtau_dt * (dev_um_n_tau[idx] - dev_um_n[idx]) + dtau * dev_res_u[idx];
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