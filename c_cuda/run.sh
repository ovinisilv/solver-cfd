#!/bin/bash
set -euo pipefail

rm -f *.mod  # remove arquivos .mod, sem erro se não existir
rm -f *.out  # remove arquivos .out
sh cleaning.sh  # executa script para limpeza extra

NVCC=nvcc
ARCH=sm_86
OPT="-O2"
STD="-std=c++11"

NVCC_PATH=$(command -v "$NVCC" || true)
if [ -z "$NVCC_PATH" ]; then
    echo "ERROR: nvcc not found in PATH" >&2
    exit 1
fi

NVCC_BINDIR="$(dirname "$NVCC_PATH")"
if [ -z "${CUDA_HOME:-}" ]; then
    CUDA_HOME="$(dirname "$(dirname "$NVCC_BINDIR")")"
    echo "Setting CUDA_HOME from nvcc path: $CUDA_HOME"
fi

echo "Using nvcc from: $NVCC_PATH"
echo "Using CUDA_HOME: $CUDA_HOME"

target_include_candidates=(
    "/usr/local/cuda/include"
    "$CUDA_HOME/include"
    "$CUDA_HOME/cuda/12.4/targets/x86_64-linux/include"
    "$NVCC_BINDIR/../include"
    "$NVCC_BINDIR/../targets/x86_64-linux/include"
)

CUDA_INC=""
for candidate in "${target_include_candidates[@]}"; do
    if [ -n "$candidate" ] && [ -d "$candidate" ] && [ -f "$candidate/cuda_runtime.h" ]; then
        CUDA_INC="$candidate"
        break
    fi
done

if [ -z "$CUDA_INC" ] && [ -n "${CUDA_HOME:-}" ]; then
    found_runtime=$(find "$CUDA_HOME" -path '*/cuda_runtime.h' 2>/dev/null | head -n 1 || true)
    if [ -n "$found_runtime" ]; then
        CUDA_INC="$(dirname "$found_runtime")"
    fi
fi

if [ -z "$CUDA_INC" ] && [ -d "$NVCC_BINDIR/../.." ]; then
    found_runtime=$(find "$NVCC_BINDIR/../.." -path '*/cuda_runtime.h' 2>/dev/null | head -n 1 || true)
    if [ -n "$found_runtime" ]; then
        CUDA_INC="$(dirname "$found_runtime")"
    fi
fi

CFLAGS=""
if [ -n "$CUDA_INC" ] && [ -d "$CUDA_INC" ]; then
    echo "Using CUDA include: $CUDA_INC"
    CFLAGS="-I$CUDA_INC"
else
    echo "WARNING: CUDA include path not found; gcc may fail to compile C files" >&2
fi

gcc $CFLAGS -c probe.c -o probe.o
gcc $CFLAGS -c convergence.c -o convergence.o
gcc $CFLAGS -c flametip.c -o flametip.o

$NVCC -arch=$ARCH $OPT $STD -lineinfo \
  xm_ym.cu \
  x_y.cu \
  vol.cu \
  solve_C.cu \
  solve_P.cu \
  solve_U.cu \
  solve_V.cu \
  solve_Z.cu \
  upwind_U_pair.cu \
  upwind_V_pair.cu \
  transient.c \
  resz.cu \
  resv.cu \
  resu.cu \
  resc.cu \
  nonsymetric_mesh.cu \
  max_reduce.cu \
  main.cu \
  initial.cu \
  grids.cu \
  dx_dy.cu \
  comum.cu \
  comp_mean.cu \
  bcZ.cu \
  bcUV.cu \
  bcP.cu \
  bcC.cu \
  area_das_faces.cu \
  output.cu \
  probe.o \
  convergence.o \
  flametip.o \
  -o cylinder_solver.out -lm

#export OMP_NUM_THREADS=8
#nohup time ./cylinder_solver.out | tee archivelog.log

#sh posv.sh

