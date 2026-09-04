%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nathan Wei and Adina Fleisher
% FAST Lab, Princeton University
% Dynamic Induction Control project
%
% Function to compute the phase average and standard deviation of a given
% input signal.
% Assumes t is uniformly sampled
%
% Inputs: data vector, number of samples in phase average, 
%         fraction of offset (e.g. 30 degrees / 360)
% Output: phase-averaged vector, standard deviation of phase-averaged
%         vector
% Dependencies: none
% Created: 06 June 2024
% Updated: 06 June 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = phase_average(y, N, offset)

y_pa = zeros(size(y));
y_pa = y_pa(1:N); % preserve row/column vector
y_pa_std = zeros(size(y_pa));

% Phase-averaging FUN TIMES!!!
for ii = 1 : N
    inds = ii - round(offset*N) : N : length(y); % all instances at given time step
    inds(inds < 1) = []; % erase negative indices
    y_pa(ii) = mean(y(inds), 'omitnan');
    y_pa_std(ii) = std(y(inds), [], 'omitnan');
end

varargout{1} = y_pa;
if nargout > 1
    varargout{2} = y_pa_std;
end