#!/bin/bash
#SBATCH --job-name=tn33
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=12:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --array=0-99
#SBATCH --mail-user=yz4281@princeton.edu

julia run_var_exp.jl 3 3 500 2 100 data/230403/230403_d1_$SLURM_ARRAY_TASK_ID