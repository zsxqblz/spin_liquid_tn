#!/bin/bash
#SBATCH --job-name=tn
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=2:00:00           # total run time limit (HH:MM:SS)
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --array=0-9
#SBATCH --mail-user=yz4281@princeton.edu

julia run_fix_exp.jl 10 2 100 data/230403/230403_d3_$SLURM_ARRAY_TASK_ID