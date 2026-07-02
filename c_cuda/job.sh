#!/bin/bash
#$ -S /bin/bash
#$ -q gpu.q
#$ -N cylinder_solver
#$ -cwd
#$ -o saida_job.txt
#$ -e erro_job.txt
set -euo pipefail

echo "PWD: $(pwd)"
echo "HOSTNAME: $(hostname)"
echo "USER: $(whoami)"
echo "DATE: $(date)"
echo "PATH=$PATH"
command -v nvcc || true
command -v module || true
if command -v module >/dev/null 2>&1; then
    module list 2>&1 | cat
    module avail 2>&1 | grep -E 'cuda|nvhpc|nvidia' || true
fi

echo "--- Running build ---"
bash ./run.sh

echo "--- Build finished, making binary executable ---"
chmod +x ./cylinder_solver.out

echo "--- Running executable ---"
bash ./executar.sh