% Nathan Wei
% Proportional-resonant control for FOWTs
% New version of calculateGains.m: the beta->phi loop gains (kp, kr, fz)
% are now found by NUMERICAL SEARCH -- maximizing phase margin of the real,
% fully-coupled beta->phi open loop (built from the actual A0/Bc matrices,
% not the isolated single-DOF reduction), subject to hard constraints on
% closed-loop stability (checked against ALL poles -- restricting this
% check to a subset, e.g. only the largest-magnitude pair, was verified by
% direct sweep to miss real instability in a large fraction of cases) and
% on |kp/kr| -- instead of the old analytical (Abbas et al. 2022-style)
% formula. The Tau_g->omega loop gains (kp_Tg, kr_Tg, fz_Tg) are UNCHANGED,
% still computed with the old analytical formula.
% Then rewrites those gains in the UMaineSemi DISCON.IN file.
% Created: 29 August 2026 (modified from calculateGains.m)

function params = optimizeGains(TSR0, beta0, Uinf, freq, freqd, zeta, ...
    stability_margin, kp_kr_ratio_max)

if nargin == 0
    TSR0     = 8.5; % optimal TSR is 8.5, design TSR is 9.0
    beta0    = 0; % (deg) optimal blade pitch is -1, design pitch is 0
    Uinf     = 8; % wind speed (m/s) - rated is 10.75 m/s
    freq     = 0.1; % desired forcing frequency (rad/s)
    freqd    = 0.1; % desired natural frequency, Tau_gen-omega controller ONLY (rad/s)
    zeta     = 1; % desired damping ratio, Tau_gen-omega controller ONLY [-]
    stability_margin = 0.02; % beta->phi search: require max(Re(ALL poles)) < -stability_margin
    kp_kr_ratio_max  = 10;   % beta->phi search: hard ceiling on |kp/kr|. Fixed kp/kr for the Tau_gen-omega controller.
end

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
[~, dCp_dTSRs] = gradient(data.Cp, data.angles, data.TSRs);
[dCq_dbetas, dCq_dTSRs] = gradient(data.Cq, data.angles, data.TSRs);
[dCt_dbetas, dCt_dTSRs] = gradient(data.Ct, data.angles, data.TSRs);
dCp_dTSR0 = interp2(data.betas, data.lambdas, dCp_dTSRs, beta0, TSR0, interpMethod);
Cp0 = interp2(data.betas, data.lambdas, data.Cp, beta0, TSR0, interpMethod);

dCq_dTSR = interp2(data.betas, data.lambdas, dCq_dTSRs, beta0, TSR0, interpMethod);
dCt_dTSR = interp2(data.betas, data.lambdas, dCt_dTSRs, beta0, TSR0, interpMethod);
dCq_dbeta = interp2(data.betas, data.lambdas, dCq_dbetas, beta0, TSR0, interpMethod);
dCt_dbeta = interp2(data.betas, data.lambdas, dCt_dbetas, beta0, TSR0, interpMethod);
Ct0 = interp2(data.betas, data.lambdas, data.Ct, beta0, TSR0, interpMethod);
dCt_dTSR0 = dCt_dTSR; dCt_dbeta0 = dCt_dbeta*180/pi; % dCt_dbeta0 needs the
    % deg->rad conversion (matches dFa_dbeta below); dCt_dTSR0 does not

%% Tau_g->omega loop gains: UNCHANGED analytical formula (Abbas et al. 2022)
A       = 1/2*rho*pi*R^4*Uinf^2 / (Jr*TSR0^2*Uinf) * (dCp_dTSR0*TSR0 - Cp0);
B_Tg    = -Ng^2 / Jr;
kp_Tg   = -1/B_Tg*(2*zeta*freqd + A);
kr_Tg   = kp_Tg/kp_kr_ratio_max;
fz_Tg   = -freqd^2 / B_Tg / kp_Tg;

%% Build the real, fully-coupled beta->phi open-loop transfer function
%% (from the actual A0/Bc matrices) -- the plant the search below checks
%% stability and phase margin against.
dTa_domega = 1/2*rho*Uinf*pi*R^4*dCq_dTSR;
dFa_domega = 1/2*rho*Uinf*pi*R^3*dCt_dTSR;
dTa_dbeta = 1/2*rho*Uinf^2*pi*R^3*dCq_dbeta * 180/pi;
dFa_dbeta = 1/2*rho*Uinf^2*pi*R^2*dCt_dbeta * 180/pi;
dTa_dU = 1/2*rho*pi*R^3*(2*Uinf ...
    * interp2(data.betas, data.lambdas, data.Cq, beta0, TSR0) ...
    - R*omega0*dCq_dTSR);
dFa_dU = 1/2*rho*pi*R^2*(2*Uinf ...
    * interp2(data.betas, data.lambdas, data.Ct, beta0, TSR0) ...
    - R*omega0*dCt_dTSR);

s = tf('s');
sys_beta_to_phi_OL = (ht*(Ng*dFa_domega*dTa_dbeta - Ng*dFa_dbeta*dTa_domega ...
    + Jr*dFa_dbeta*s))/(Jr*Kt*s - Kt*Ng*dTa_domega + Dt*Jr*s^2 + Jr*Jt*s^3 ...
    - Jt*Ng*dTa_domega*s^2 + Jr*dFa_dU*ht^2*s^2 - Dt*Ng*dTa_domega*s ...
    - Ng*dFa_dU*dTa_domega*ht^2*s + Ng*dFa_domega*dTa_dU*ht^2*s);

%% beta->phi loop gains: numerical search over (kp, fz, kr), ALL free,
%% maximizing phase margin subject to (1) closed-loop stability (ALL
%% poles) with margin, (2) 0 < fz < freq (ROSCO's own CheckInputs
%% requirement, buffered), and (3) |kp/kr| <= kp_kr_ratio_max.

% Internal-only seed constants (NOT design targets -- phase margin is the
% actual objective now -- just used to center the search and set its
% range, exactly as validated in PR_control_gain_optimization.m).
zeta_seed = 0.7;
omega_des_seed = 0.1156;
kp0 = (Jt*omega_des_seed^2*(1-4*zeta_seed^2) + 2*zeta_seed*omega_des_seed*(Dt + rho*pi*R^2*ht^2*Uinf ...
      * (Ct0 - (TSR0/2)*dCt_dTSR0)) - Kt) / (1/2*rho*pi*R^2*ht*Uinf^2*dCt_dbeta0);
ki0 = -((Dt + rho*pi*R^2*ht^2*Uinf*(Ct0 - (TSR0/2)*dCt_dTSR0))*omega_des_seed^2 ...
      - 2*Jt*zeta_seed*omega_des_seed^3) / (1/2*rho*pi*R^2*ht*Uinf^2*dCt_dbeta0);
fz0 = ki0/kp0;

fz_lo = 0.1*freq;
fz_hi = 0.9*freq;
kr0 = kp0/kp_kr_ratio_max;
fz0 = min(max(fz0, fz_lo), fz_hi); % clip the analytical seed into the valid band

gainBound = struct('kp', 5*abs(kp0), 'kr', 5*abs(kr0));

evalGains = @(x) evaluate_gains(x(1), x(2), x(3), kp_kr_ratio_max, fz_lo, fz_hi, ...
    sys_beta_to_phi_OL, freq, stability_margin, gainBound);

N = 13; % grid points per dimension for the coarse search: 13^3=2197 evaluations.
kp_range = linspace(-gainBound.kp, gainBound.kp, N);
fz_range = linspace(fz_lo, fz_hi, N);
kr_range = linspace(-gainBound.kr, gainBound.kr, N);

bestCost = inf; bestX = [kp0, fz0, kr0];
for ikp = 1:N
    for ifz = 1:N
        for ikr = 1:N
            cost = evalGains([kp_range(ikp), fz_range(ifz), kr_range(ikr)]);
            if cost < bestCost
                bestCost = cost;
                bestX = [kp_range(ikp), fz_range(ifz), kr_range(ikr)];
            end
        end
    end
end

opts = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-10, 'MaxFunEvals', 8000, 'MaxIter', 8000);
xopt = fminsearch(evalGains, bestX, opts);

[~, feasibleOpt, diagOpt] = evaluate_gains(xopt(1), xopt(2), xopt(3), ...
    kp_kr_ratio_max, fz_lo, fz_hi, sys_beta_to_phi_OL, freq, stability_margin, gainBound);
kp = xopt(1); fz = xopt(2); kr = xopt(3);

fprintf('=== beta->phi gains (numerical search) ===\n');
fprintf('kp=%.6f, fz=%.6f, kr=%.6f  (|kp/kr|=%.3f, feasible=%d)\n', kp, fz, kr, abs(kp/kr), feasibleOpt);
fprintf('max(Re(ALL poles))=%.6f, Phase margin=%.2f deg, Gain margin=%.2f dB\n', ...
    diagOpt.maxRe, diagOpt.PM, diagOpt.GM_dB);
if ~feasibleOpt
    warning('optimizeGains:infeasible', ...
        'Numerical search did not converge to a feasible (stable) beta->phi design -- widen the search range or relax stability_margin.');
end

params = struct('kp', kp, 'kr', kr, 'kp_Tg', kp_Tg, 'kr_Tg', kr_Tg, 'freqz', ...
    fz, 'freqz_Tg', fz_Tg, 'omega0', omega0, 'freq', freq);

%% Rewrite the calculated gains, omega_z, omega0, and freq into the UMaineSemi DISCON.IN file
disconFile = fullfile(fileparts(mfilename('fullpath')), 'Examples', 'Test_Cases', ...
    'IEA-15-240-RWT', 'IEA-15-240-RWT-UMaineSemi', 'IEA-15-240-RWT-UMaineSemi_DISCON.IN');
writeDisconGains(disconFile, params);

end

%% ------------------------------------------------------------------
function [cost, feasible, diagnostics] = evaluate_gains(kp, fz, kr, ...
    kp_kr_ratio_max, fz_lo, fz_hi, sys_beta_to_phi_OL, freq, stability_margin, gainBound)
% Evaluate a candidate (kp,fz,kr) -- ALL free -- for the beta->phi loop
% against the real, fully-coupled plant. Returns a scalar cost (lower is
% better) = -PhaseMargin(loop_bp) for feasible candidates, plus a soft box
% penalty keeping |kp|,|kr| within gainBound.kp/.kr of the analytical seed.
% Hard constraints, each penalized like instability if violated:
%   (1) closed-loop stability, checked against ALL poles, with at least
%       stability_margin of clearance from the imaginary axis
%   (2) fz_lo < fz < fz_hi
%   (3) |kp/kr| <= kp_kr_ratio_max
s = tf('s');
diagnostics = struct('maxRe', NaN, 'PM', NaN, 'GM_dB', NaN, 'Wcp', NaN, 'boxPenalty', 0);
if abs(kp) < 1e-8
    cost = 1e8; feasible = false; return;
end
boxPenalty = max(0, abs(kp)-gainBound.kp)^2 + max(0, abs(kr)-gainBound.kr)^2;
diagnostics.boxPenalty = boxPenalty;

if fz <= fz_lo || fz >= fz_hi
    distToValid = min(abs(fz-fz_lo), abs(fz-fz_hi));
    cost = 1e6 + distToValid*1e2 + boxPenalty;
    feasible = false;
    return;
end

ratio = abs(kp) / max(abs(kr), 1e-10);
if ratio > kp_kr_ratio_max
    cost = 1e6 + (ratio - kp_kr_ratio_max)*1e2 + boxPenalty;
    feasible = false;
    return;
end

try
    controller_beta = -(kp*(1+fz/s) + kr*s/(s^2+freq^2));
    loop_bp = minreal(sys_beta_to_phi_OL * controller_beta);
    sys_bp  = minreal(loop_bp / (1 + loop_bp));
    p = pole(sys_bp); % ALL closed-loop poles -- checked in full, not a subset
    maxRe = max(real(p));
    diagnostics.maxRe = maxRe;
    feasible = maxRe < -stability_margin;
    if ~feasible
        cost = 1e6 + (maxRe + stability_margin)*1e4 + boxPenalty;
        return;
    end

    [Gm, Pm, ~, Wcp] = margin(loop_bp);
    if isempty(Pm) || isnan(Pm) || isinf(Pm)
        Pm = 180;
    end
    if isempty(Gm) || isnan(Gm) || isinf(Gm)
        GM_dB = 40;
    else
        GM_dB = 20*log10(Gm);
    end
    diagnostics.PM = Pm;
    diagnostics.GM_dB = GM_dB;
    diagnostics.Wcp = Wcp;

    cost = -Pm + boxPenalty;
catch
    cost = 1e8;
    feasible = false;
end
end

%% ------------------------------------------------------------------
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
