#!/bin/bash
#SBATCH --job-name=annotation_single_cell    ### Job Name
#SBATCH --output=annotation_single_cell.log       ### File in which to store job output
#SBATCH --error=annotation_single_cell.err        ### File in which to store job error messages
#SBATCH --qos=normal        ### Quality of Service (like a queue in PBS)
#SBATCH --partition=centos7  ### Partition to run on (not needed with normal and long queues)
#SBATCH --time=1-00:00:00     ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1             ### Node count required for the job
#SBATCH --ntasks-per-node=20   ### Number of tasks to be launched per Node
#SBATCH --mem=256000

source ~/.bashrc
source activate st1
Rscript 1_Annotation_cell_level.r
