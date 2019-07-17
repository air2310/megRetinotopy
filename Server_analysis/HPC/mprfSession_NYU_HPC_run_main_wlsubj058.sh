#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time 20:00:00 # time (D-HH:MM)
#SBATCH --mem=64GB # memory pool for all cores
#SBATCH --job-name=mprf_main
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ek99@nyu.edu
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

module load matlab/2016b

cd /scratch/ek99

# If the files you are running are not in the same folder as this script,
# you can insert "addpath(genpath('/PATH/TO/FILES/'));" before the command
# you want to run.
matlab -nodisplay -r "addpath(genpath(brainstorm3)); addpath(genpath(meg_utils)); addpath(genpath(noisepoolPCADenoise)); addpath(genpath(vistasoft)); addpath(fieldtrip); ft_defaults; cd('megRetinotopy'); mprf_addPaths; run('runModelMain$SLURM_ARRAY_TASK_ID('wlsubj058')); exit()"

exit

