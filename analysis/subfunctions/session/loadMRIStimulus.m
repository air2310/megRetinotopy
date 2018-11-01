function retStim = loadMRIStimulus(stimPths)

% subfunction of loadStim, to specifically prepare the MRI stimulus

% Options for vistasoft's rmStimulusMatrix function
useFinalImages = false;
scans          = false;
xRange         = [];
yRange         = [];

% Get MRI prf scan parameters from file
prfParams = load(stimPths.PRFParams.pth);

% [EK]: not sure what the difference is between unconvolved and all
% stim images. And why we even need this?
%         retStim.im = prfParams.params.stim.images_unconvolved;
%         retStim.im_conv = prfParams.params.analysis.allstimimages';

retStim.window  = prfParams.params.stim.stimwindow;
X = prfParams.params.analysis.X;
Y = prfParams.params.analysis.Y;

% Create X and Y axis of stimulus grid
[retStim.X, retStim.Y] = meshgrid(unique(X), unique(Y));

% Get final retinotopy stimulus images (previously called full_im)
[~, retStim.im] = rmStimulusMatrix(prfParams.params,xRange,yRange,useFinalImages,scans);

end