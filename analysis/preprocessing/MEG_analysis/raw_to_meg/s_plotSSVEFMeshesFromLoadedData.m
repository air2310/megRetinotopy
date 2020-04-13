% s_plotSSVEFMeshesFromLoadedData.m

% Define subject ID
subjs = {'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};
for ii=1:length(subjs)
close all;

subjID = subjs{ii};
%%
% Load paths with data files for this subject
dirPth = loadPaths(subjID);

load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'epoched_data_hp_preproc_denoised.mat'), 'data');
load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'megStimConditions.mat'), 'triggers');
opt = getOpts('saveFig',1,'verbose',1);

data = data.data;

dataPermuted = permute(data, [4 1 2 3]);

[nSensors, nTimePoints, nEpochs, nRuns] = size(dataPermuted);

dataReshaped =reshape(dataPermuted, [nSensors, nTimePoints, nEpochs*nRuns]);

plotSSVEFmesh(dataReshaped, triggers.stimConditions, subjID, dirPth, opt)
clear data
end