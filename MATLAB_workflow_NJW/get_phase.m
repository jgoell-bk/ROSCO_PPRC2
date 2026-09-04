%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nathan Wei
% Princeton University
% FAST Lab: Dynamic Induction Control project
%
% Function to compute spectrum of a quantity
%
% Inputs: time vector, data vector
% Output: spectrum and frequency vectors
% Dependencies: none
% Created: 06 August 2024
% Updated: 06 August 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phase, amp] = get_phase(t, x)

% Compute amplitude and phase using fft
% Code from https://www.mathworks.com/matlabcentral/answers/275358-...
fts = fft(x)/length(t);
Fv = linspace(0, 1, fix(length(t)/2)+1)*(1/mean(diff(t)))/2; % freq. vector
Iv = 1:length(Fv); % index vector
amp_fts = abs(fts(Iv))*2; % note: first element is doubled
phs_fts = angle(fts(Iv)); % in rad
fft_amp = amp_fts(2); % only works cuz we're FFTing a single period
fft_phs = wrapTo2Pi(phs_fts(2) + pi/2); % sine phase, in nondimensional time units
phase = wrapTo180(rad2deg(fft_phs)); % phase in degrees
amp = fft_amp;

end