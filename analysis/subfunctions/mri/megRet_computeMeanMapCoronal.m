function [] = megRet_computeMeanMapCoronal(subject, sessionDir)

% Computes mean maps and averaged scans in the Inplane. Transforms 
% averaged mean map and coranal to gray using trilinear interpolation.

% Open the session
cd(sessionDir)

% Open Vista
vw = mrVista;

%% Compute mean maps for original scans

vw = computeMeanMap(vw, 0);
updateGlobal(vw)

%view meap map
vw = setDisplayMode(vw, 'map');
vw = viewSet(vw, 'cothresh', 0);
vw = refreshScreen(vw);


%% Average 

% Make a new dataTYPE of the averages  
vw = averageTSeries(vw, 1:6, 'Averages', 'Average of Original 1:6');

vw = viewSet(vw, 'current data type', 'Averages'); 


%% Compute mean map on averaged scan
vw = computeMeanMap(vw, 0);
updateGlobal(vw)

%view meap map
vw = setDisplayMode(vw, 'map');
vw = viewSet(vw, 'cothresh', 0);
vw = refreshScreen(vw);

%% Xform maps to gray
gr = mrVista('3');
gr = viewSet(gr, 'cur dt', 'Averages');
vw = viewSet(vw, 'cur dt', 'Averages');

% Xform mean map
gr = ip2volParMap(vw, gr, 0, [], 'linear');
gr = viewSet(gr, 'cothresh', 0);
gr = setDisplayMode(gr, 'map');
gr = refreshScreen(gr);

gr = ip2volTSeries(vw, gr, [], 'linear');


% Save parameters
saveSession;
 











