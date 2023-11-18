#!/bin/bash
#SBATCH -p checkpt
#SBATCH -t 3-00:00:00
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -o Build_ref_run.out
#SBATCH -e Build_ref_run.err
source Build_ref.sh