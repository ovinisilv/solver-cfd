#include "comum.h"

__global__ void col_p(double *dev_um_tau, double *dev_vm_tau, double *dev_p, double *dev_pn, int *dev_flag, double *dev_rc, int imax, int jmax, double beta_param, double dt_dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if(i <= imax && j <= jmax){
        int idx_local = i * (jmax + 1) + j;
        if(dev_flag[idx_local] == c_f){
            
            // Garantindo os divisores e propriedades locais de malha via variáveis explicitadas
            double dx = dx_c[i];
            double dy = dy_c[j];

            dev_rc[idx_local] = ((dev_um_tau[(i)*(jmax+1)+j] - dev_um_tau[(i-1)*(jmax+1)+j]) / dx + 
                                 (dev_vm_tau[i*(jmax+2)+j] - dev_vm_tau[i*(jmax+2)+j-1]) / dy);

            dev_p[idx_local] = dev_pn[idx_local] - beta_param * dt_dtau * dev_rc[idx_local];
        } else {
            dev_p[idx_local] = 1.0;
            dev_rc[idx_local] = 0.0;
        }
    }
}

void solve_P(double *dev_um_tau, double *dev_vm_tau, double *dev_p, double *dev_pn, int *dev_flag, double *dev_rc){
    dim3 blockDim(16, 16);
    dim3 gridDimPtc((imax + blockDim.x - 1)/blockDim.x, (jmax + blockDim.y - 1)/blockDim.y);

    double dt_dtau = dt / (dt + dtau);

    // Passando explicitamente a variável global beta para o kernel
    col_p<<<gridDimPtc, blockDim>>>(dev_um_tau, dev_vm_tau, dev_p, dev_pn, dev_flag, dev_rc, imax, jmax, beta, dt_dtau);
    cudaDeviceSynchronize();
}