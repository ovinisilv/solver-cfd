#include "comum.h"
#define idx i*(jmax+1)+j

__global__ void kernel_V(double *dev_um, double *dev_vm, double *dev_vm_n, double *dev_vm_tau, double *dev_vm_n_tau, double *dev_p, double *dev_t, int *dev_flag, double *dev_rv, double *dev_vi, double *dev_res_v, int imax, int jmax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if(i <= imax && j <= jmax + 1){
        int idx_v = i*(jmax+2)+j;
        if(dev_flag[idx] == c_f){
            double dtau_dt = dtau / dt;

            double v_atual = dev_vm_n_tau[idx_v];
            double v_atras = dev_vm_n_tau[idx_v-1];
            double v_frente = dev_vm_n_tau[idx_v+1];
            double v_esquerda = dev_vm_n_tau[(i-1)*(jmax+2)+j];
            double v_direita = dev_vm_n_tau[(i+1)*(jmax+2)+j];

            double dv_dx = (v_direita - v_esquerda) / (2.0 * dev_dx[i]);
            double dv_dy = (v_frente - v_atras) / (2.0 * dev_dy[j]);
            double d2v_dx2 = (v_direita - 2.0*v_atual + v_esquerda) / (dev_dx[i]*dev_dx[i]);
            double d2v_dy2 = (v_frente - 2.0*v_atual + v_atras) / (dev_dy[j]*dev_dy[j]);

            double nu = 0.01;

            dev_res_v[idx_v] = -(v_atual * dv_dy) + nu * (d2v_dx2 + d2v_dy2) - (dev_p[idx] - dev_p[i*(jmax+1)+j-1]) / dev_dy[j];

            dev_vm_tau[idx_v] = dev_vm_n_tau[idx_v] - dtau_dt * (dev_vm_n_tau[idx_v] - dev_vm_n[idx_v]) + dtau * dev_res_v[idx_v];
        } else {
            dev_vm_tau[idx_v] = (1.0e-3) * v_i;
            dev_res_v[idx_v] = 0.0;
        }
    }
}

void solve_V(double *dev_um, double *dev_vm, double *dev_vm_n, double *dev_vm_tau, double *dev_vm_n_tau, double *dev_p, double *dev_t, int *dev_flag, double *dev_rv, double *dev_vi, double *dev_res_v){
    dim3 blockDim(16, 16);
    dim3 gridDimVm((imax + blockDim.x - 1)/blockDim.x, (jmax + 1 + blockDim.y - 1)/blockDim.y);

    kernel_V<<<gridDimVm, blockDim>>>(dev_um, dev_vm, dev_vm_n, dev_vm_tau, dev_vm_n_tau, dev_p, dev_t, dev_flag, dev_rv, dev_vi, dev_res_v, imax, jmax, dt, dtau);
    cudaDeviceSynchronize();
}