//Artificial Compressibility Methods
//Solve Momentum Equation with QUICK Scheme
//BASED ON:
//Versteeg, H. K., and W. Malalasekera. 
//"An introduction to computational Fluid Dynamics, The finite volume control, ed." (1995).
#include "comum.h"
#include <math.h>
#define N_IMAX (51*2)
#define N_ITC 800000
#define idx i*(jmax+1)+j

__global__ void atualizar_matrizes_linearizadas(double *origem, double *destino, int tamanhoLinha, int tamanhoColuna, int inicio, int coluna){
    int i = blockIdx.x * blockDim.x + threadIdx.x + inicio;
    int j = blockIdx.y * blockDim.y + threadIdx.y + inicio;
    
    if(i <= tamanhoLinha && j <= tamanhoColuna){
        destino[i*(coluna)+j] = origem[i*(coluna)+j];
    }
}

__global__ void atualiza_tc(double *dev_t, double *dev_c, int *dev_flag, double temp_cylinder, double concentracao_inicial, int c_f, int imax, int jmax){
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    if(i > imax || j > jmax) return;
    
    if(dev_flag[idx] != c_f){
        dev_t[idx] = temp_cylinder;
        dev_c[idx] = concentracao_inicial;
    }
}

int main(int argc, char const *argv[]){
    //--- iterations structured type ---
    Iteracoes iterations;
    iterations.itc_max = N_ITC;
    iterations.final_time = 0.50; 

    printf("\n==============================================");
    printf("\n       CFD POROUS MEDIA - CUDA SOLVER        ");
    printf("\n==============================================\n");

    //lendo parametros de entrada
    input_parameters();

    //alocando variaveis globais
    alocar_globais();

    //alocando ponteiros locais por causa do escopo
    double *dev_rc, *dev_ru, *dev_ui, *dev_res_u, *dev_rv, *dev_vi, *dev_res_v, *dev_rz;
    double *dev_um, *dev_vm, *dev_um_n, *dev_um_tau, *dev_um_n_tau, *dev_vm_n, *dev_vm_tau, *dev_vm_n_tau;
    double *dev_pn, *dev_h, *dev_z, *dev_t_n_tau, *dev_t_tau, *dev_c_n_tau, *dev_c_tau;
    double *dev_u, *dev_v, *dev_p, *dev_t, *dev_c;

    cudaMalloc((void**)&dev_rc, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_ru, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_ui, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_res_u, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_rv, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_vi, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_res_v, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_rz, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_um, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_vm, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_um_n, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_um_tau, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_um_n_tau, sizeof(double)*(imax+2)*(jmax+1));
    cudaMalloc((void**)&dev_vm_n, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_vm_tau, sizeof(double)*(imax+1)*(jmax+2));
    cudaMalloc((void**)&dev_vm_n_tau, sizeof(double)*(imax+1)*(jmax+2));

    cudaMalloc((void**)&dev_pn, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_h, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_z, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_t_n_tau, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_t_tau, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_c_n_tau, sizeof(double)*(imax+1)*(jmax+1));
    cudaMalloc((void**)&dev_c_tau, sizeof(double)*(imax+1)*(jmax+1));

    cudaMallocManaged((void**)&dev_u, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_v, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_p, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_t, sizeof(double)*(imax+1)*(jmax+1));
    cudaMallocManaged((void**)&dev_c, sizeof(double)*(imax+1)*(jmax+1));

    // ================================================================
    // INICIALIZAÇÃO DE BUFFERS CONTRA LIXO DE MEMÓRIA (NaN)
    // ================================================================
    cudaMemset(dev_rc, 0, sizeof(double)*(imax+1)*(jmax+1));
    cudaMemset(dev_ru, 0, sizeof(double)*(imax+2)*(jmax+1));
    cudaMemset(dev_ui, 0, sizeof(double)*(imax+2)*(jmax+1));
    cudaMemset(dev_res_u, 0, sizeof(double)*(imax+2)*(jmax+1));
    cudaMemset(dev_rv, 0, sizeof(double)*(imax+1)*(jmax+2));
    cudaMemset(dev_vi, 0, sizeof(double)*(imax+1)*(jmax+2));
    cudaMemset(dev_res_v, 0, sizeof(double)*(imax+1)*(jmax+2));
    cudaMemset(dev_rz, 0, sizeof(double)*(imax+1)*(jmax+1));
    cudaMemset(dev_pn, 0, sizeof(double)*(imax+1)*(jmax+1));
    cudaMemset(dev_h, 0, sizeof(double)*(imax+1)*(jmax+1));
    cudaMemset(dev_z, 0, sizeof(double)*(imax+1)*(jmax+1));
    cudaMemset(dev_t_n_tau, 0, sizeof(double)*(imax+1)*(jmax+1));
    cudaMemset(dev_t_tau, 0, sizeof(double)*(imax+1)*(jmax+1));
    cudaMemset(dev_c_n_tau, 0, sizeof(double)*(imax+1)*(jmax+1));
    cudaMemset(dev_c_tau, 0, sizeof(double)*(imax+1)*(jmax+1));

    //--- CUDA Events para profiling de tempo ---
    cudaEvent_t ev_ts_start, ev_ts_stop;
    cudaEvent_t ev_solveU_start, ev_solveU_stop;
    cudaEvent_t ev_solveV_start, ev_solveV_stop;
    cudaEvent_t ev_solveP_start, ev_solveP_stop;
    cudaEvent_t ev_solveZ_start, ev_solveZ_stop;
    cudaEvent_t ev_solveC_start, ev_solveC_stop;

    cudaEventCreate(&ev_ts_start);    cudaEventCreate(&ev_ts_stop);
    cudaEventCreate(&ev_solveU_start); cudaEventCreate(&ev_solveU_stop);
    cudaEventCreate(&ev_solveV_start); cudaEventCreate(&ev_solveV_stop);
    cudaEventCreate(&ev_solveP_start); cudaEventCreate(&ev_solveP_stop);
    cudaEventCreate(&ev_solveZ_start); cudaEventCreate(&ev_solveZ_stop);
    cudaEventCreate(&ev_solveC_start); cudaEventCreate(&ev_solveC_stop);

    float time_ts, time_solveU, time_solveV, time_solveP, time_solveZ, time_solveC;
    float acc_ts=0.0f, acc_solveU=0.0f, acc_solveV=0.0f, acc_solveP=0.0f, acc_solveZ=0.0f, acc_solveC=0.0f;

    //gerando malha
    mesh();

    //condicoes iniciais
    IC(dev_um, dev_vm, dev_u, dev_v, dev_p, dev_t, dev_c, dev_flag);

    //--- Execution configuration ---
    dim3 blockDim(16, 16);
    dim3 gridDimUm( (imax+1 + blockDim.x - 1) / blockDim.x, (jmax + blockDim.y - 1) / blockDim.y );
    dim3 gridDimVm( (imax + blockDim.x - 1) / blockDim.x, (jmax+1 + blockDim.y - 1) / blockDim.y );
    dim3 gridDimPtc( (imax + blockDim.x - 1) / blockDim.x, (jmax + blockDim.y - 1) / blockDim.y );

    //--- FORÇANDO INICIALIZAÇÃO CORRETA DE EMISSÃO DE TEMPO ---
    dt = 0.005;
    dtau = 0.05;
    beta = 1.0;

    int k = 0;
    double tempo = 0.0;

    printf("\nIniciando simulação transiente...\n");
    
    //--- Physical time step ---
    while(tempo < iterations.final_time){
        tempo = tempo + dt;

        // Salva o passo físico anterior antes das sub-iterações
        atualizar_matrizes_linearizadas<<<gridDimUm, blockDim>>>(dev_um, dev_um_n, imax+1, jmax, 1, jmax+1);
        atualizar_matrizes_linearizadas<<<gridDimVm, blockDim>>>(dev_vm, dev_vm_n, imax, jmax+1, 1, jmax+2);
        
        // Inicializa estimativas iniciais do pseudo-tempo (tau)
        atualizar_matrizes_linearizadas<<<gridDimUm, blockDim>>>(dev_um, dev_um_tau, imax+1, jmax, 1, jmax+1);
        atualizar_matrizes_linearizadas<<<gridDimVm, blockDim>>>(dev_vm, dev_vm_tau, imax, jmax+1, 1, jmax+2);
        atualizar_matrizes_linearizadas<<<gridDimPtc, blockDim>>>(dev_t, dev_t_tau, imax, jmax, 1, jmax+1);
        atualizar_matrizes_linearizadas<<<gridDimPtc, blockDim>>>(dev_c, dev_c_tau, imax, jmax, 1, jmax+1);

        cudaEventRecord(ev_ts_start);
        acc_solveU = 0.0f; acc_solveV = 0.0f; acc_solveP = 0.0f;
        acc_solveZ = 0.0f; acc_solveC = 0.0f;

        int itc = 0;
        double erro_max = 1.0;

        //--- Pseudo time step (Dual-Time Stepping Loop) ---
        while(itc < iterations.itc_max){
            itc = itc + 1;

            // Sincroniza passo n_tau anterior do pseudo-tempo para cálculo correto do resíduo temporal
            atualizar_matrizes_linearizadas<<<gridDimUm, blockDim>>>(dev_um_tau, dev_um_n_tau, imax+1, jmax, 1, jmax+1);
            atualizar_matrizes_linearizadas<<<gridDimVm, blockDim>>>(dev_vm_tau, dev_vm_n_tau, imax, jmax+1, 1, jmax+2);
            atualizar_matrizes_linearizadas<<<gridDimPtc, blockDim>>>(dev_t_tau, dev_t_n_tau, imax, jmax, 1, jmax+1);
            atualizar_matrizes_linearizadas<<<gridDimPtc, blockDim>>>(dev_c_tau, dev_c_n_tau, imax, jmax, 1, jmax+1);

            //--- SOLVE MOMENTUM X ---
            cudaEventRecord(ev_solveU_start);
            solve_U(dev_um, dev_vm, dev_um_n, dev_um_tau, dev_um_n_tau, dev_p, dev_t, dev_flag, dev_ru, dev_ui, dev_res_u);
            cudaEventRecord(ev_solveU_stop);
            cudaEventSynchronize(ev_solveU_stop);
            cudaEventElapsedTime(&time_solveU, ev_solveU_start, ev_solveU_stop);
            acc_solveU += time_solveU;

            //--- SOLVE MOMENTUM Y ---
            cudaEventRecord(ev_solveV_start);
            solve_V(dev_um, dev_vm, dev_vm_n, dev_vm_tau, dev_vm_n_tau, dev_p, dev_t, dev_flag, dev_rv, dev_vi, dev_res_v);
            cudaEventRecord(ev_solveV_stop);
            cudaEventSynchronize(ev_solveV_stop);
            cudaEventElapsedTime(&time_solveV, ev_solveV_start, ev_solveV_stop);
            acc_solveV += time_solveV;

            //--- SOLVE PRESSURE CORRECTION ---
            cudaEventRecord(ev_solveP_start);
            solve_P(dev_um_tau, dev_vm_tau, dev_p, dev_pn, dev_flag, dev_rc);
            cudaEventRecord(ev_solveP_stop);
            cudaEventSynchronize(ev_solveP_stop);
            cudaEventElapsedTime(&time_solveP, ev_solveP_start, ev_solveP_stop);
            acc_solveP += time_solveP;

            //--- SOLVE ENERGY (CORRIGIDO: Passando os buffers pseudo-temporais certos) ---
            cudaEventRecord(ev_solveZ_start);
            solve_Z(dev_um_tau, dev_vm_tau, dev_t_tau, dev_t, dev_t_n_tau, dev_z, dev_h, dev_flag, dev_rz);
            cudaEventRecord(ev_solveZ_stop);
            cudaEventSynchronize(ev_solveZ_stop);
            cudaEventElapsedTime(&time_solveZ, ev_solveZ_start, ev_solveZ_stop);
            acc_solveZ += time_solveZ;

            //--- SOLVE SPECIES (CORRIGIDO) ---
            cudaEventRecord(ev_solveC_start);
            solve_C(dev_um_tau, dev_vm_tau, dev_c_tau, dev_c, dev_c_n_tau, dev_flag, dev_rc);
            cudaEventRecord(ev_solveC_stop);
            cudaEventSynchronize(ev_solveC_stop);
            cudaEventElapsedTime(&time_solveC, ev_solveC_start, ev_solveC_stop);
            acc_solveC += time_solveC;

            erro_max = erro(dev_rc, dev_res_u, dev_res_v);

            if(erro_max < 1.0e-5){
                break;
            }
        } 

        // Atualiza os arrays principais com o resultado convergido da sub-iteração
        atualizar_matrizes_linearizadas<<<gridDimUm, blockDim>>>(dev_um_tau, dev_um, imax+1, jmax, 1, jmax+1);
        atualizar_matrizes_linearizadas<<<gridDimVm, blockDim>>>(dev_vm_tau, dev_vm, imax, jmax+1, 1, jmax+2);
        atualizar_matrizes_linearizadas<<<gridDimPtc, blockDim>>>(dev_t_tau, dev_t, imax, jmax, 1, jmax+1);
        atualizar_matrizes_linearizadas<<<gridDimPtc, blockDim>>>(dev_c_tau, dev_c, imax, jmax, 1, jmax+1);

        atualiza_tc<<<gridDimPtc, blockDim>>>(dev_t, dev_c, dev_flag, temp_cylinder, concentracao_inicial, c_f, imax, jmax);

        velocidade_centro_celula(dev_um, dev_vm, dev_u, dev_v);

        k = k + 1;
        printf("Passo: %d | Tempo Físico: %.4f s | Sub-iterações: %d | Resíduo: %.4e\n", k, tempo, itc, erro_max);

        if(k % 1 == 0){
            output(dev_um, dev_vm, dev_u, dev_v, dev_p, dev_t, dev_c, k);
        }

        cudaEventRecord(ev_ts_stop);
        cudaEventSynchronize(ev_ts_stop);
        cudaEventElapsedTime(&time_ts, ev_ts_start, ev_ts_stop);
        acc_ts += time_ts;
    } 

    cudaEventDestroy(ev_ts_start);    cudaEventDestroy(ev_ts_stop);
    cudaEventDestroy(ev_solveU_start); cudaEventDestroy(ev_solveU_stop);
    cudaEventDestroy(ev_solveV_start); cudaEventDestroy(ev_solveV_stop);
    cudaEventDestroy(ev_solveP_start); cudaEventDestroy(ev_solveP_stop);
    cudaEventDestroy(ev_solveZ_start); cudaEventDestroy(ev_solveZ_stop);
    cudaEventDestroy(ev_solveC_start); cudaEventDestroy(ev_solveC_stop);

    cudaFree(dev_um);         cudaFree(dev_um_n);
    cudaFree(dev_um_tau);     cudaFree(dev_um_n_tau);
    cudaFree(dev_vm);         cudaFree(dev_vm_n);
    cudaFree(dev_vm_tau);     cudaFree(dev_vm_n_tau);
    cudaFree(dev_u);          cudaFree(dev_v);
    cudaFree(dev_rc);         cudaFree(dev_ru);
    cudaFree(dev_ui);         cudaFree(dev_res_u);
    cudaFree(dev_rv);         cudaFree(dev_vi);
    cudaFree(dev_res_v);      cudaFree(dev_rz);
    cudaFree(dev_pn);         cudaFree(dev_h);
    cudaFree(dev_z);          cudaFree(dev_t_n_tau);
    cudaFree(dev_t_tau);      cudaFree(dev_c_n_tau);
    cudaFree(dev_c_tau);
    cudaFree(dev_p);          cudaFree(dev_t);          cudaFree(dev_c);

    desalocar_globais();

    printf("\nSimulação concluída com sucesso!\n");
    return 0;
}