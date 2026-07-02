#!/bin/bash
#$ -S /bin/bash
#$ -q gpu.q
#$ -N cylinder_solver
#$ -cwd
#$ -o saida_job.txt
#$ -e erro_job.txt

# Usando o bash explicitamente ignora travas de permissão do arquivo
bash run.sh
bash executar.sh