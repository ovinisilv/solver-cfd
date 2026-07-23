# Relatório Técnico — Solver CFD CUDA

## Resumo
Este documento descreve a implementação, arquitetura, métodos numéricos, entrada/saída, e instruções de execução do projeto `solver-cfd` presente neste diretório. O solver resolve escoamentos incompressíveis usando o método de compressibilidade artificial, com aceleração por GPU via CUDA.

## 1. Objetivo do projeto
- Simular escoamento incompressível em malhas 2D (possibilidade de malha simétrica ou não-simétrica).
- Resolver equações de momento, continuidade, energia e transporte de espécie (mistura/concentração).
- Acelerar computação usando kernels CUDA e Unified Memory (`cudaMallocManaged`).

## 2. Arquitetura do código
- Entrada principal: `main.cu` — inicialização, loop temporal (físico e pseudo), chamadas para os solvers e geração de saída.
- Configuração global: `comum.cu`, `comum.h` — parâmetros físicos, malha, alocação/desalocação de arrays e estruturas globais.
- Solvers modulares (cada um em arquivo separado):
  - `solve_U.cu` — componente U do momentum (QUICK / esquema implícito-pseudo-tempo).
  - `solve_V.cu` — componente V do momentum.
  - `solve_P.cu` — pressão/continuidade por compressibilidade artificial.
  - `solve_Z.cu` — transporte da fração de mistura (Z) usando Ralston RK2.
  - `solve_C.cu` — transporte de concentração (C) usando Ralston RK2.
- Condições de contorno: `bcUV.cu`, `bcP.cu`, `bcC.cu` — aplicam contornos por kernels.
- Utilitários de malha/volumes/áreas: `area_das_faces.cu`, `grids.cu`, `vol.cu`, `x_y.cu`, `xm_ym.cu`, `dx_dy.cu`.
- Protótipos e interfaces: `functions.h`.

## 3. Principais parâmetros e entradas
- Malha: controlada por `imax` (número de pontos em x) e `jmax` calculado a partir de `lhori`, `y_up`, `y_down` e `dx_c`.
- Parâmetros numéricos: definidos na estrutura `Iterations` (`itc_max`, `eps`, `dtau_f`, `final_time`, etc.). Esses podem ser lidos de arquivos em `input/` (ex.: `iteration.dat`) ou ajustados no código.
- Propriedades de combustível/fluido: em `properties/` (ex.: `properties.f90`) e possivelmente em `input/reference.dat`.
- Restart e dados: `data/restart*` contém arquivos binários de restart; `data/results/` armazena saída por variáveis (U, V, PTZH).

## 4. Métodos numéricos
- Formulação: método de compressibilidade artificial para acoplar velocidade e pressão, conforme referência clássica (Versteeg & Malalasekera).
- Discretização espacial: esquemas QUICK para termos convectivos e opções de upwind.
- Integração temporal:
  - Pseudo-tempo (iteração interna) para convergência das equações de momento/pressão.
  - Passo físico `dt` para evolução temporal real.
  - Integrador de Ralston (RK2) para equações escalares de transporte (`Z` e `C`).
- Critério de convergência: máximo residual entre `residual_u`, `residual_v` e `residual_p` comparado a `iterations.eps`.

## 5. Implementação CUDA e paralelismo
- Uso intensivo de kernels CUDA para cálculos locais e aplicação de BCs.
- Uso de `cudaMallocManaged` para Unified Memory: simplifica movimentação H2D/D2H, porém pode implicar page faults e overhead de migração.
- Estratégia de blocos/threads: kernels usam blocos 2D (`dim3 blockDim(16,16)`), e grades dimensionadas com `imax`/`jmax`.
- Sincronizações explícitas via `cudaEvent` para medir tempos de kernels e seccionalmente via `cudaDeviceSynchronize` quando necessário.

## 6. Compilação e execução
- Recomendado usar `nvcc`. Exemplo mínimo para compilar (ajuste `-arch=sm_XX` conforme GPU):

```bash
cd /home/pep/Documents/solver-cfd/c_cuda
nvcc -O3 -arch=sm_70 -o cylinder_solver main.cu *.cu
```

- Executar com argumentos opcionais `imax` e `itc_max`:

```bash
./cylinder_solver [imax] [itc_max]
```

- Scripts existentes: `run.sh`, `executar.sh`, `executar2.sh` — verifique e adapte conforme seu ambiente.

## 7. Saída e pós-processamento
- Saída principal: arquivos em `data/results/` e `transient/` para frames de animação.
- Plot/anim: scripts `pos*.sh`, `animation.sh` e `transient.c` (usam `convert` e ferramentas de plotagem). Atenção que `convert` (ImageMagick) pode consumir muita RAM.

## 8. Observações e recomendações
- Criar um `Makefile` ou `CMakeLists.txt` para simplificar compilação com flags corretas de arquitetura e debug.
- Incluir README atualizado com exemplos completos de execução e descrição dos arquivos em `input/`.
- Considerar perfilamento com `nsys` para diagnosticar migrações de página da Unified Memory.
- Validar e documentar formatos binários em `data/restart*` para interoperabilidade.

## 9. Referências
- H. K. Versteeg, W. Malalasekera — "An Introduction to Computational Fluid Dynamics" (1995).
- Código fonte no diretório atual.

## 10. Apêndice — arquivos principais
- `main.cu` — controlador principal.
- `comum.cu`, `comum.h` — parâmetros globais, alocação.
- `functions.h` — protótipos das rotinas.
- `solve_*.cu` — `solve_U.cu`, `solve_V.cu`, `solve_P.cu`, `solve_Z.cu`, `solve_C.cu`.
- `bcUV.cu`, `bcP.cu`, `bcC.cu` — condições de contorno.
- Scripts: `run.sh`, `executar.sh`, `pos*.sh`, `animation.sh`, `cleaning.sh`.

---

*Gerado automaticamente por análise do código no repositório `solver-cfd`.*
