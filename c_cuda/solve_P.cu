#include "comum.h"
#define idx i*(jmax+1)+j

__global__ void col_p(double *dev_um_tau, double *dev_vm_tau, double *dev_p, double *dev_pn, int *dev_flag, double *dev_rc, int imax, int jmax, double beta, double dt_dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if(i <= imax && j <= jmax){
        if(dev_flag[idx] == c_f){
            // Correção aqui: dev_rc mapeado de forma consistente com a alocação da main (jmax+1)
            dev_rc[idx] = ((dev_um_tau[(i)*(jmax+1)+j] - dev_um_tau[(i-1)*(jmax+1)+j])/dx_c[i] + 
                           (dev_vm_tau[i*(jmax+2)+j] - dev_vm_tau[i*(jmax+2)+j-1])/dy_c[j]);

            dev_p[idx] = dev_pn[idx] - beta * dt_dtau * dev_rc[idx];
        } else {
            dev_p[idx] = 1.0;
            dev_rc[idx] = 0.0;
        }
    }
}

void solve_P(double *dev_um_tau, double *dev_vm_tau, double *dev_p, double *dev_pn, int *dev_flag, double *dev_rc){
    dim3 blockDim(16, 16);
    dim3 gridDimPtc((imax + blockDim.x - 1)/blockDim.x, (jmax + blockDim.y - 1)/blockDim.y);

    double dt_dtau = dt / (dt + dtau);

    col_p<<<gridDimPtc, blockDim>>>(dev_um_tau, dev_vm_tau, dev_p, dev_pn, dev_flag, dev_rc, imax, jmax, beta, dt_dtau);
    cudaDeviceSynchronize();
}