function phaseRefMEGResponse = mprf_MEGPhaseReferenceData(megData, predMEGResponse)
% Function to computing phase referenced amplitude from preprocessed MEG data 
% and predicted MEG responses from cortical surface
%   phaseRefMEGResponse = mprf_MEGPhaseReferenceData(megData, predMEGResponse)
%   	
% INPUTS:
%   megData         : preprocessed MEG data (time x epochs x run x sensors)
%   predMEGResponse : predicted MEG responses (epochs x sensors)
%  
% OUTPUT:
%   phaseRefMEGResponse : Phase referenced MEG time series (sensors x epochs)



phaseRange = linspace(0,pi,20); % range of values to search for the reference phase

% to do:
% - FFT of MEG data
% - Get phase & amplitude at 10 Hz
% - Pick reference phase that results in positive amplitudes, and thus
% positive beta values for stimulus peaks in predicted MEG response





return


