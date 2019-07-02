function G_constrained = loadGainMtx(subjID, dirPth, opt)
%
% Function to load gain matrix produced by Brainstorm

if opt.verbose; sprintf('(%s) Load MEG Gain matrix for subject %s...\n', mfilename, subjID); end

if opt.fullSizeGainMtx
    d = dir(fullfile(dirPth.bs.dataPth, '*', 'headmodel*meg_02.mat'));
else
    d = dir(fullfile(dirPth.bs.dataPth, '*', 'headmodel*meg.mat'));
end

headmodel = load(fullfile(d.folder,d.name));

% Keep all sensors in Gain matrix
keep_sensors = logical([ones(1,157), zeros(1,size(headmodel.Gain,1)-157)]);

% Get Gain matrix and truncate to not-nan sensors
G = headmodel.Gain(keep_sensors,:); % [Nsensors x 3*Nvertices]

% Contrained gain matrix
G_constrained = bst_gain_orient(G, headmodel.GridOrient); % [Nsensors x Nsources], equivalent to size BS pial cortex [1x15002]

end