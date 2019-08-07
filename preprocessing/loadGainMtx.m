function G_constrained = loadGainMtx(subjID, dirPth, opt)
% Function to load gain matrix produced by Brainstorm
%
%   G_constrained = loadGainMtx(subjID, dirPth, opt)
%
% INPUTS:
%   subjID          : subject name (string)
%   dirPth          : paths to files for given subject (struct)
%   opt             : pipeline options (struct with boolean flags).
%
% OUTPUT:
%   G_constrained   : constrained Gain matrix, i.e. only dipoles with
%                       perpendicular dipoles (relative to vertex)
%                       (sensors x vertices)
%
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019


if opt.verbose; sprintf('(%s) Load MEG Gain matrix for subject %s...\n', mfilename, subjID); end

if opt.fullSizeMesh
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