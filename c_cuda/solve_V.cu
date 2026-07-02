#include "comum.h"
#define idx i*(jmax+1)+j

__global__ void kernel_V(double *dev_um, double *dev_vm, double *dev_vm_n, double *dev_vm_tau, double *dev_vm_n_tau, double *dev_p, double *dev_t, int *dev_flag, double *dev_rv, double *dev_vi, double *dev_res_v, int imax, int jmax, double dt, double dtau){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if(i <= imax && j <= jmax + 1){
        int idx_v = i*(jmax+2)+j;
        if(dev_flag[idx] == c_f){
            double dtau_dt = dtau / (dt + 1.0e-12);
            
            double dy = dy_c[j];
            if(dy < 1.0e-9) dy = 1.0e-3;

            // Proteções contra indeterminação matemática local
            double v_atual = dev_vm_n_tau[idx_v];
            if(isnan(v_atual) || isinf(v_atual)) v_atual = 0.0;

            double v_atras = dev_vm_n_tau[idx_v-1];
            if(isnan(v_atras) || isinf(v_atras)) v_atras = 0.0;

            double v_frente = dev_vm_n_tau[idx_v+1];
            if(isnan(v_frente) || isinf(v_frente)) v_frente = 0.0;

            double conveccao = v_atual * (v_atual - v_atras) / dy;
            double difusao = (v_frente - 2.0*v_atual + v_atras) / (dy * dy);

            double p_atual = dev_p[idx];
            double p_atras = dev_p[i*(jmax+1)+j-1];
            if(isnan(p_atual)) p_atual = 1.0;
            if(isnan(p_atras)) p_atras = 1.0;

            dev_res_v[idx_v] = -conveccao + (1.0/100.0)*difusao - (p_atual - p_atras)/dy;

            double v_proximo = v_atual - (dtau_dt) * (v_atual - dev_vm_n[idx_v]) + dtau * dev_res_v[idx_v];
            
            if(isnan(v_proximo) || isinf(v_proximo)) v_proximo = (1.0e-3) * v_i;

            dev_vm_tau[idx_v] = v_proximo;
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