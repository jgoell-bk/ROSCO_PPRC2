% Nathan Wei
% Proportional-resonant control for FOWTs
% Script calculates gains based on desired input parameters
% Then rewrites those gains in the UMaineSemi_DISCON.IN file
% Created: 18 Aug 2026 (modified from PR_control_nonlinear_sim.m)

function params = calculateGains(TSR0, beta0, Uinf, freq, zeta, kp_kr_ratio)

if nargin == 0
    TSR0    = 8.5; % optimal TSR is 8.5, design TSR is 9.0
    beta0   = 0; % (deg) optimal blade pitch is -1, design pitch is 0
    Uinf    = 10; % wind speed (m/s) - rated is 10.75 m/s
    freq    = 0.1; % desired forcing frequency (rad/s)
    zeta    = 0.7; % desired damping ratio [-]
end

fz = freq/10; % desired PI controller zero frequency (rad/s)
R  = 120; % turbine radius (m)

%% Turbine parameters
Ng = 1; % gearbox ratio
Jr = 3.525e8; % turbine moment of inertia (kg m^2) - IEA 15 MW
rho = 1.2; % air density (kg/m^3)
interpMethod = 'linear'; % 2D interpolation method for quasi-steady aerodynamics

%% Sensitivity coefficients from IEA 15 MW steady data
data = load('IEA15MW_Cp_Ct_Cq.mat');
omega0 = TSR0*Uinf/R;
[dCp_dbetas, dCp_dTSRs] = gradient(data.Cp, data.angles, data.TSRs);
dCp_dTSR = interp2(data.betas, data.lambdas, dCp_dTSRs, beta0, TSR0, interpMethod); % value at setpoint
dCp_dbeta = interp2(data.betas, data.lambdas, dCp_dbetas, beta0, TSR0, interpMethod) * 180/pi; % value at setpoint (data are in deg, convert to rad -- matches dTa_dbeta/dFa_dbeta convention below)
Cp0 = interp2(data.betas, data.lambdas, data.Cp, beta0, TSR0, interpMethod);

%% Calculate gains using modified formulae from Abbas et al. (2022)
A       = 1/2*rho*pi*R^4*Uinf^2 / (Jr*TSR0^2*Uinf) * (dCp_dTSR*TSR0 - Cp0);
B_beta  = Ng/(2*Jr*TSR0^2) * rho*pi*R^3*Uinf^2 * (dCp_dbeta*TSR0);
B_Tg    = -Ng^2 / Jr;
kp      = -1/(2*pi*B_beta)*(2*zeta*freq + A); % for beta-phi PR controller
kr      = kp/kp_kr_ratio; % for beta-phi PR controller
kp_Tg   = -1/B_Tg*(2*zeta*freq + A);
kr_Tg   = kp_Tg/kp_kr_ratio;

params = struct('kp', kp, 'kr', kr, 'kp_Tg', kp_Tg, 'kr_Tg', kr_Tg, 'omega_z', ...
    fz, 'omega0', omega0, 'freq', freq);

%% Rewrite the calculated gains, omega_z, omega0, and freq into the UMaineSemi DISCON.IN file
disconFile = fullfile(fileparts(mfilename('fullpath')), 'Examples', 'Test_Cases', ...
    'IEA-15-240-RWT', 'IEA-15-240-RWT-UMaineSemi', 'IEA-15-240-RWT-UMaineSemi_DISCON.IN');
writeDisconGains(disconFile, params);

end

function writeDisconGains(disconFile, params)
% Rewrite the PPPR_freq_phi, PPPR_freq_omega, PPPR_offset_omega,
% PPPR_CntrGains_phi, and PPPR_CntrGains_omega lines in disconFile with the
% newly calculated values, leaving everything else in the file -- including
% the trailing "! ParamName - description" comments -- untouched.

lines = strsplit(fileread(disconFile), '\n', 'CollapseDelimiters', false);

lines = setDisconValue(lines, 'PPPR_freq_phi', sprintf('%.6g', params.freq));
lines = setDisconValue(lines, 'PPPR_freq_omega', sprintf('%.6g', params.freq));
lines = setDisconValue(lines, 'PPPR_offset_omega', sprintf('%.6g', params.omega0));
lines = setDisconValue(lines, 'PPPR_CntrGains_phi', ...
    sprintf('%.6g   %.6g  %.6g', params.kp, params.kr, params.omega_z));
lines = setDisconValue(lines, 'PPPR_CntrGains_omega', ...
    sprintf('%.6g  %.6g  %.6g', params.kp_Tg, params.kr_Tg, params.omega_z));

fid = fopen(disconFile, 'w');
fwrite(fid, strjoin(lines, '\n'));
fclose(fid);

end

function lines = setDisconValue(lines, paramName, newValueStr)
% Replace the value portion (before "!") of the line whose trailing
% comment names paramName (e.g. "! PPPR_offset_omega - ..."), preserving
% that comment text as-is.
marker = ['! ', paramName];
idx = find(contains(lines, marker), 1);
if isempty(idx)
    error('setDisconValue:notFound', 'Could not find "%s" in %s', paramName, 'DISCON.IN');
end
commentStart = strfind(lines{idx}, '!');
lines{idx} = [newValueStr, '    ', lines{idx}(commentStart(1):end)];
end
