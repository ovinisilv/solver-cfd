#!/bin/bash
#$ -S /bin/bash
#$ -q gpu.q
#$ -N cylinder_solver
#$ -cwd
#$ -o saida_job.txt
#$ -e erro_job.txt

bash ./run.sh
bash ./executar.sh