% Nathan Wei
% Proportional-resonant control for FOWTs
% Script calculates gains based on desired input parameters
% Then rewrites those gains in the UMaineSemi_DISCON.IN file
% Created: 18 Aug 2026 (modified from PR_control_nonlinear_sim.m)

function params = calculateGains(TSR0, beta0, Uinf, freq, freqd, zeta, ...
    freqd_Tg, zeta_Tg, kp_kr_ratio)

if nargin == 0
    TSR0     = 8.5; % optimal TSR is 8.5, design TSR is 9.0
    beta0    = 0; % (deg) optimal blade pitch is -1, design pitch is 0
    Uinf     = 10; % wind speed (m/s) - rated is 10.75 m/s
    freq     = 0.1; % desired forcing frequency (rad/s)
    freqd    = 0.1; % desired natural frequency, beta-phi controller (rad/s)
    zeta     = 0.7; % desired damping ratio, beta-phi controller [-]
    freqd_Tg = 0.1; % desired natural frequency, Tau_gen-omega controller (rad/s)
    zeta_Tg  = 1; % desired damping ratio, Tau_gen-omega controller [-]
    kp_kr_ratio = 10; % somewhat arbitrary but can be theoretically justified
end

% fz = freq/10; % desired PI controller zero frequency (rad/s)
R  = 120; % turbine radius (m)

%% Turbine parameters
Ng = 1; % gearbox ratio
Jr = 3.525e8; % turbine moment of inertia (kg m^2) - IEA 15 MW
ht = 150; % hub height of turbine (m)
Jt = 1.251e10; % tower moment of inertia (kg m^2)
rho = 1.2; % air density (kg/m^3)
mt = 249718; % tower (OC3 spar) mass (kg)
Dt = 4.7e7; % platform damping coefficient (e.g. hydrodynamics)
Kt = 5.98e8; % platform restoring coefficient (e.g. from mooring lines)
% Dt and Kt fit from single-DOF (PtfmPDOF only) OpenFAST free-decay tests,
% no wind, all other DOFs off, fit to
% phi(t) = C + exp(-sigma*t)*(a*cos(wd*t)+b*sin(wd*t)).
% Kt=5.98e8 confirmed to 4 sig figs across two independent tests (phi0=10
% deg with waves, phi0=2 deg no waves) -- wn=0.2187 rad/s both times,
% matching DISCON.IN's stated platform pitch natural freq (0.213 rad/s)
% to ~3%. Dt is amplitude-dependent (nonlinear/quadratic damping): the
% phi0=10deg test gave Dt=9.0e7, phi0=2deg gave Dt=4.7e7 (both tskip=20s
% fit window). 4.7e7 is used here since PPPR_amp_phi=2deg is close to
% the actual operating oscillation amplitude.
interpMethod = 'linear'; % 2D interpolation method for quasi-steady aerodynamics

%% Sensitivity coefficients from IEA 15 MW steady data
data = load('IEA15MW_Cp_Ct_Cq.mat');
omega0 = TSR0*Uinf/R;
[dCp_dbetas, dCp_dTSRs] = gradient(data.Cp, data.angles, data.TSRs);
% [dCq_dbetas, dCq_dTSRs] = gradient(data.Cq, data.angles, data.TSRs);
[dCt_dbetas, dCt_dTSRs] = gradient(data.Ct, data.angles, data.TSRs);
% dTa_dUs = 1/2*rho*pi*R^3*(2*Uinf*data.Cq - R*omega0*dCq_dTSRs); % sensitivity of aero torque to wind speed (kg m/s)
% dTa_dU0 = interp2(data.betas, data.lambdas, dTa_dUs, beta0, TSR0, interpMethod); % value at setpoint
dCp_dTSR0 = interp2(data.betas, data.lambdas, dCp_dTSRs, beta0, TSR0, interpMethod); % value at setpoint
% dCp_dbeta0 = interp2(data.betas, data.lambdas, dCp_dbetas, beta0, TSR0, interpMethod) * 180/pi; % value at setpoint (data are in deg, convert to rad -- matches dTa_dbeta/dFa_dbeta convention below)
Cp0 = interp2(data.betas, data.lambdas, data.Cp, beta0, TSR0, interpMethod);

dCt_dTSR0 = interp2(data.betas, data.lambdas, dCt_dTSRs, beta0, TSR0, interpMethod);
dCt_dbeta0 = interp2(data.betas, data.lambdas, dCt_dbetas, beta0, TSR0, interpMethod) * 180/pi; % value at setpoint (data are in deg, convert to rad -- matches dTa_dbeta/dFa_dbeta convention below)
Ct0 = interp2(data.betas, data.lambdas, data.Ct, beta0, TSR0, interpMethod);

%% Calculate gains using modified formulae from Abbas et al. (2022)
A       = 1/2*rho*pi*R^4*Uinf^2 / (Jr*TSR0^2*Uinf) * (dCp_dTSR0*TSR0 - Cp0);
% B_beta  = Ng/(2*Jr*TSR0^2) * rho*pi*R^3*Uinf^2 * (dCp_dbeta0*TSR0);
B_Tg    = -Ng^2 / Jr;
% kp      = -1/(B_beta)*(2*zeta*freqd + A); % for beta-phi PR controller
kp      = -(Jt*freqd^2*(1-4*zeta^2) + 2*zeta*freqd*(Dt + rho*pi*R^2*ht^2*Uinf ...
          * (Ct0 - (TSR0/2)*dCt_dTSR0)) - Kt) / (1/2*rho*pi*R^2*ht*Uinf^2*dCt_dbeta0);
ki      = -abs((Dt + rho*pi*R^2*ht^2*Uinf*(Ct0 - (TSR0/2)*dCt_dTSR0))*freqd^2 ...
          - 2*Jt*zeta*freqd^3) / (1/2*rho*pi*R^2*ht*Uinf^2*dCt_dbeta0);
fz      = ki/kp;
kr      = kp/kp_kr_ratio; % for beta-phi PR controller
% k_Tg  = m_Tg*ht/Ng*dTa_dU0; % Eqn. 29, from Fischer (2013) and Stockhouse et al. (2021)
kp_Tg   = -1/B_Tg*(2*zeta_Tg*freqd_Tg + A);
kr_Tg   = kp_Tg/kp_kr_ratio;
fz_Tg   = -freqd_Tg^2 / B_Tg / kp_Tg;

params = struct('kp', kp, 'kr', kr, 'kp_Tg', kp_Tg, 'kr_Tg', kr_Tg, 'freqz', ...
    fz, 'freqz_Tg', fz_Tg, 'omega0', omega0, 'freq', freq);

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
    sprintf('%.6g   %.6g  %.6g', params.kp, params.kr, params.freqz));
lines = setDisconValue(lines, 'PPPR_CntrGains_omega', ...
    sprintf('%.6g  %.6g  %.6g', params.kp_Tg, params.kr_Tg, params.freqz_Tg));

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
