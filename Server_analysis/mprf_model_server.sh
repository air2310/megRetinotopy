#!/bin/bash
 #
 #$ -N creating parameters for the MEG model to run
 #$ -S /bin/bash
 #$ -j Y
 #$ -q long.q
 #$ -o /home/edadan/meg_test.o
 #$ -u edadan
 #$ -M a.edadan@uu.nl 
 
 if [ -z "$1" ]
 then
  echo 'Run_mprf_model. Inputs:'
  echo 'save_path_model_params=$1, path to save model parameters' 
  echo 'MEG_data_file=$2, Preprocessed MEG data'
  echo 'sub_sess_dir=$3, path to mprfSESSION'
  echo 'model_type=$4, subject session directory'
  echo 'variable parameters'
  exit 1
 fi

 # Load environment module 
 # module load matlab/R2018b
 echo "making params file for running MEG model"
 echo "Running $4"
 
 echo "$1" 
 echo "$2" 
 echo "$3"
 echo "$4"
 echo "$5"
 echo "$6"
 echo "$7" 

 # Run matlab job
 matlab -nodesktop -nosplash -nodisplay  -r "addpath(genpath('/home/akhi/Documents/MATLAB/toolboxes/')); mprfSession_model_server('$1','$2','$3','$4','n_iteration_scrambled','$5','n_cores','$6','phase_fit_loo','$7'); quit"


