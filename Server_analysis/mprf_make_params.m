function [mprf_params,params_file] = mprf_make_params(save_path,model_type,varargin)
% mprf_make_params - Creates parameters that are required to run the
% mprfSession_run_model. 
%
% Input - save_path: path to save the output
%       - model_type: Which model to run (Eg: "original model phase ref amplitude" for running the model with the phase referenced amplitude fits,
%                                             "original phase ref amplitude model fit leave one out" for model with phase reference amplitude fit and leave one out calculation of the reference phase)
% Output - params_file: Contains all the parameters for running the model     

global mprfSESSION

if ~exist('save_path','var') || ~exist(save_path,'file')
    save_path = pwd;
end

% Adding the variable arguments to model
if nargin > 2
    addArg = varargin;
    if numel(addArg) == 1
       addArg = addArg{1};
    end
else
    addArg = [];
end

model = mprf_generate_model_params(model_type,addArg);


%tmp_stimulus = load((fullfile(pwd,'/stimuli/meg/imported_stimulus/meg_stimulus.mat')));
tmp_stimulus = load(fullfile(mprfSESSION.meg.imported_stim,'/meg_stimulus.mat'));
stimulus  = tmp_stimulus.meg_stim;

channels.data = 0:156;
channels.trigger = 160:167;
channels.diode = 191;

syn.do = 0; 
syn.pp = 0;

mprf_params.model = model;
mprf_params.stimulus = stimulus;
mprf_params.channels = channels;
mprf_params.syn = syn;

%cur_time = datestr(now);
%cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
params_file = fullfile(save_path,strcat('mprf_params','.mat'));
save(params_file,'model','stimulus','syn','channels');

end