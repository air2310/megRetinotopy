#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time 20:00:00 # time (D-HH:MM)
#SBATCH --mem=64GB # memory pool for all cores
#SBATCH --job-name=megRet_wlsubj068
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ek99@nyu.edu
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

module load matlab/2016b

cd /scratch/ek99/megRetData/wlsubj068

# If the files you are running are not in the same folder as this script,
# you can insert "addpath(genpath('/PATH/TO/FILES/'));" before the command
# you want to run.
matlab -nodisplay -r "addpath(genpath('/scratch/ek99/megRetinotopy/')); addpath(genpath('/scratch/ek99/megRetData/')); addpath(genpath('/scratch/ek99/vistasoft')); addpath('/scratch/ek99/fieldtrip'); ft_defaults; cd('/scratch/ek99/megRetData/wlsubj068'); run('runModel$SLURM_ARRAY_TASK_ID.m'); exit()"

exit

