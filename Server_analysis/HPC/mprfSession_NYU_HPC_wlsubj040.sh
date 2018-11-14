#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time 20:00:00 # time (D-HH:MM)
#SBATCH --mem=64GB # memory pool for all cores
#SBATCH --job-name=megRetRunModel_wlsubj040
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ek99@nyu.edu
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

module load matlab/2016b

cd /scratch/ek99/megRetinotopy/

# If the files you are running are not in the same folder as this script,
# you can insert "addpath(genpath('/PATH/TO/FILES/'));" before the command
# you want to run.
matlab -nodisplay -r "addpath(genpath('/scratch/ek99/megRetinotopy')); cd('/scratch/ek99/'); sub_sess_dir = '/scratch/ek99/megRetData/Subjects_sessions/wl_subj040'; save_path = fullfile(sub_sess_dir, 'modeling', 'results'); meg_data_file = fullfile(sub_sess_dir, 'data', 'meg', 'raw'); mprfSession_model_server(save_path,meg_data_file,sub_sess_dir,$SLURM_ARRAY_TASK_ID); exit()"

exit

