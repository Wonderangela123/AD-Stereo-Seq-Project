#!/bin/bash
#SBATCH --job-name=D02175A4    ### Job Name
#SBATCH --output=D02175A4.out       ### File in which to store job output
#SBATCH --error=D02175A4.err        ### File in which to store job error messages
#SBATCH --qos=normal
#SBATCH --time=1-00:00:00     ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1             ### Node count required for the job
#SBATCH --ntasks-per-node=20   ### Number of tasks to be launched per Node
#SBATCH --partition=centos7
#SBATCH --mem=128000
module load anaconda3/2023.07
export CONDA_ENVS_PATH=/lustre/project/hdeng2/ygong/env
source activate st
python segmentation.py
