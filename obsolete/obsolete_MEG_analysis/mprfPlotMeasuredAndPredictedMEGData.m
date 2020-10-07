%function out = mprfPlotMeasuredAndPredictedMEGData(data, pred, periods, pred_stim);


% TODO:
% 1. make GUI showing scalp activity, per position for both prediciton
% and data -> takes long to create these images, beter save them
% 2. make GUI that allows users to select relevant channels and plot
% prediction with data for these channels.
% 3. Make sure both GUI can exist simultaneously




if strcmpi(measure,'amp_ratio')
    if exist(fullfile(pwd, 'mprfSESSION.mat'),'file')
        load(fullfile(pwd, 'mprfSESSION.mat'))
    end
    
    data_scalp_im_02 = mprfCreateScalpImages(data_pp_02); 
    pred_scalp_im_02 = mprfCreateScalpImages(pred_pp_02);
    stim_im = pred_stim.full_im;
    
    mprfDataPredScalp_gui(data_scalp_im.snr_median, pred_scalp_im.snr_median, stim_im);
    
else
    error('Unrecognized measure type')
    
    
    
end



%end






























