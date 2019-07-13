function amps_out = rescaleAmpsWithRefPhase(amps_in, phase_in, ref_phase)
% Function to scale amplitude data with the cosine of a particular angle
%  amps_out = rescaleAmpsWithRefPhase(amps_in, phase_in, ref_phase);
%
% INPUTS:
%   amps_in     : Amplitudes of spectral data
%   phase_in    : Phase  of spectral data (radians)
%   ref_phase   : Phase reference, to take compare with phase (radians)
%
% OUTPUTS:
%   amps_out    : Scaled amplitudes 

% Compute difference between phase and reference phase
phase_diff  = phase_in - ref_phase;

% Compute vector length of phase difference
scale_factor = cos(phase_diff);

% Rescale amplitudes with phase vector
amps_out    = amps_in.*scale_factor;

end