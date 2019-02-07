function mprfSessionWrapper(subject, model_type)
% Function that loads subject paths and runs the requested model

% NOTE THIS FUNCTION IS STILL UNDER CONSTRUCTIOIN -- DO NOT USE

% mprfSessionWrapper(subject, modelType)

% INPUTS
% subject       : subjectID (string)
% model_type    : string with an integer (default = '1');
%                    '1' - standard / original model (comparing data with prediction)
%                    '2' - prf size range model (comparing data with model prediction for different prf sizes)
%                    '3' - prf position range model (comparing data with model prediction for different prf positions)
%                    '4' - splif half reliability check using standard / original model
%                    '5' - scramble prf parameters across vertices (keeping x,y,sigma together)

% Load model parameters
model = getModelParams(model_type);

% Load relevant data structures: meg data, prf parameters per vertex, lead field
data  = getData(subject);

% Get MEG stimulus
stimulus = getStimulus(subject);

% Make prediction for predicted vertex amplitudes
prediction.MRI = predictPRFResponse(model, data, stimulus);
    
% Make prediction for predicted MEG phase references amplitudes
prediction.MEG = predictedMEGResponse(bs, predicted_MRI_resp);

if savePredictions
    % Get current time and day for file name
    curTime = datestr(now); curTime(curTime == ' ' | curTime == ':' | curTime == '-') = '_';
    
    % Get directory to save predictions
    save_dir = mprf__get_directory('model_predictions');
    
    % Save it!
    save(fullfile(main_dir, save_dir, ['model_predictions_' curTime]),...
        'prf','bs','roi','model','stimulus','pred_resp','syn','meg_resp','channels');
end


% Run the model:
switch model_type
    case '1'
        mprfSession_run_original(prediction);
    case '2'
        mprfSession_run_prf_size_range;
    case '3'
        mprfSession_run_position_range;
    case '4'
        mprfSession_run_splithalf_reliability;
    case '5'
        mprfSession_run_scramble;
end