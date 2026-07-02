#include "comum.h"

// Garantindo a declaração explícita do escopo global das matrizes e parâmetros para o NVCC
extern __managed__ double *dx_c;
extern __managed__ double *dy_c;
extern __managed__ double beta;

__global__ void col_p(double *dev_um_tau, double *dev_vm_tau, double *dev_p, double *dev_pn, int *dev_flag, double *dev_rc, double *dev_dx, double *dev_dy, int imax, int jmax, double beta_param, double dt_dtau) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if (i <= imax && j <= jmax) {
        int idx_local = i * (jmax + 1) + j;
        if (dev_flag[idx_local] == c_f) {
            
            // Lendo dos ponteiros de malha passados explicitamente para o Kernel
            double dx = dev_dx[i];
            double dy = dev_dy[j];

            dev_rc[idx_local] = ((dev_um_tau[(i) * (jmax + 1) + j] - dev_um_tau[(i - 1) * (jmax + 1) + j]) / dx + 
                                 (dev_vm_tau[i * (jmax + 2) + j] - dev_vm_tau[i * (jmax + 2) + j - 1]) / dy);

            dev_p[idx_local] = dev_pn[idx_local] - beta_param * dt_dtau * dev_rc[idx_local];
        } else {
            dev_p[idx_local] = 1.0;
            dev_rc[idx_local] = 0.0;
        }
    }
}

void solve_P(double *dev_um_tau, double *dev_vm_tau, double *dev_p, double *dev_pn, int *dev_flag, double *dev_rc) {
    dim3 blockDim(16, 16);
    dim3 gridDimPtc((imax + blockDim.x - 1) / blockDim.x, (jmax + blockDim.y - 1) / blockDim.y);

    double dt_dtau = dt / (dt + dtau);

    // Passando os vetores globais dx_c e dy_c diretamente como argumentos para evitar problemas de escopo oculto na GPU
    col_p<<<gridDimPtc, blockDim>>>(dev_um_tau, dev_vm_tau, dev_p, dev_pn, dev_flag, dev_rc, dx_c, dy_c, imax, jmax, beta, dt_dtau);
    cudaDeviceSynchronize();
}