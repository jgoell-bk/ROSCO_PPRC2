% Nathan Wei / Claude Code
% Proportional-resonant control for FOWTs
% beta->phi loop gains (kp, kr, fz): numerical search maximizing phase
% margin of the real, coupled beta->phi open loop (built from the actual
% A0/Bc matrices), subject to (1) closed-loop stability (ALL poles) with
% margin, (2) fz < freq (ROSCO's CheckInputs requirement), and (3)
% |kp/kr| <= kp_kr_ratio_max. Tau_g->omega loop gains (kp_Tg, kr_Tg,
% fz_Tg): unchanged old analytical formula (Abbas et al. 2022).
% Returns the gains as a struct.
% Created: 18 Aug 2026 (modified from PR_control_nonlinear_sim.m)
% Updated: 3 Sep 2026 -- beta->phi gains now from numerical search
% Updated: 3 Sep 2026 -- DISCON.IN writing moved out to setupSim.m

function params = calculateGains(TSR0, beta0, Uinf, freq, zeta, ...
    stability_margin, kp_kr_ratio_max)

if nargin == 0
    TSR0     = 8.5; % optimal TSR is 8.5, design TSR is 9.0
    beta0    = 0; % (deg) optimal blade pitch is -1, design pitch is 0
    Uinf     = 8; % wind speed (m/s) - rated is 10.75 m/s
    freq     = 0.1; % desired forcing frequency (rad/s)
    zeta     = 1; % desired damping ratio, Tau_gen-omega controller ONLY [-]
    stability_margin = 0.02; % beta->phi search: require max(Re(ALL poles)) < -stability_margin
    kp_kr_ratio_max  = 10;   % beta->phi search ceiling on |kp/kr|; also fixed kp/kr for Tau_gen-omega
end
freqd = freq; % desired natural frequency, Tau_gen-omega controller ONLY (rad/s)
R  = 120; % turbine radius (m)

%% Turbine parameters
Ng = 1; % gearbox ratio
Jr = 3.525e8; % turbine moment of inertia (kg m^2) - IEA 15 MW
ht = 150; % hub height of turbine (m)
Jt = 1.251e10; % tower moment of inertia (kg m^2)
rho = 1.2; % air density (kg/m^3)
% mt = 249718; % tower (OC3 spar) mass (kg)
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

%% Tau_g->omega loop gains: unchanged analytical formula (Abbas et al. 2022)
A     = 1/2*rho*pi*R^4*Uinf^2 / (Jr*TSR0^2*Uinf) * (dCp_dTSR0*TSR0 - Cp0);
B_Tg  = -Ng^2 / Jr;
kp_Tg = -1/B_Tg*(2*zeta*freqd + A);
kr_Tg = kp_Tg/kp_kr_ratio_max;
fz_Tg = -freqd^2 / B_Tg / kp_Tg;

%% Real, coupled beta->phi open-loop transfer function (from A0/Bc matrices)
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

%% beta->phi loop gains: numerical search over (kp,fz,kr), ALL free,
%% maximizing phase margin subject to (1) closed-loop stability (ALL
%% poles) with margin, (2) fz < freq, and (3) |kp/kr| <= kp_kr_ratio_max.
kp0 = 0.5; % fixed seed, not the final gains
fz0 = 0.1*freq;
fz_hi = freq*(1 - 1e-6); % HARD CONSTRAINT: fz < freq, tiny fp safety margin
fz_grid_lo = 1e-3*freq;  % grid seeding range only, not a constraint
kr0 = kp0/kp_kr_ratio_max;

gainBound = struct('kp', 5*abs(kp0), 'kr', 5*abs(kr0));

evalGains = @(x) evaluate_gains(x(1), x(2), x(3), kp_kr_ratio_max, fz_hi, ...
    sys_beta_to_phi_OL, freq, stability_margin, gainBound);

N = 13; % grid points per dimension: 13^3=2197 evaluations
kp_range = linspace(-gainBound.kp, gainBound.kp, N);
fz_range = linspace(fz_grid_lo, fz_hi, N);
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
    kp_kr_ratio_max, fz_hi, sys_beta_to_phi_OL, freq, stability_margin, gainBound);
kp = xopt(1); fz = xopt(2); kr = xopt(3);

fprintf('=== beta->phi gains (numerical search) ===\n');
fprintf('kp=%.6f, fz=%.6f, kr=%.6f  (|kp/kr|=%.3f, feasible=%d)\n', kp, fz, kr, abs(kp/kr), feasibleOpt);
fprintf('max(Re(ALL poles))=%.6f, Phase margin=%.2f deg, Gain margin=%.2f dB\n', ...
    diagOpt.maxRe, diagOpt.PM, diagOpt.GM_dB);
if ~feasibleOpt
    warning('calculateGains:infeasible', ...
        'Numerical search did not converge to a feasible (stable) beta->phi design -- widen the search range or relax stability_margin.');
end

params = struct('kp', kp, 'kr', kr, 'kp_Tg', kp_Tg, 'kr_Tg', kr_Tg, 'freqz', ...
    fz, 'freqz_Tg', fz_Tg, 'omega0', omega0, 'freq', freq);

end

%% ------------------------------------------------------------------
function [cost, feasible, diagnostics] = evaluate_gains(kp, fz, kr, ...
    kp_kr_ratio_max, fz_hi, sys_beta_to_phi_OL, freq, stability_margin, gainBound)
% Cost = -PhaseMargin for feasible (kp,fz,kr), plus a soft box penalty on
% kp,kr. Hard constraints, each penalized like instability if violated:
%   (1) closed-loop stability, checked against ALL poles, with at least
%       stability_margin of clearance from the imaginary axis
%   (2) fz < fz_hi
%   (3) |kp/kr| <= kp_kr_ratio_max
s = tf('s');
diagnostics = struct('maxRe', NaN, 'PM', NaN, 'GM_dB', NaN, 'Wcp', NaN);
if abs(kp) < 1e-8
    cost = 1e8; feasible = false; return; % avoid degenerate near-zero-gain controllers
end
boxPenalty = max(0, abs(kp)-gainBound.kp)^2 + max(0, abs(kr)-gainBound.kr)^2;

if fz >= fz_hi
    cost = 1e6 + abs(fz-fz_hi)*1e2 + boxPenalty;
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
    p = pole(sys_bp); % ALL closed-loop poles, not a subset
    maxRe = max(real(p));
    diagnostics.maxRe = maxRe;
    feasible = maxRe < -stability_margin;
    if ~feasible
        cost = 1e6 + (maxRe + stability_margin)*1e4 + boxPenalty;
        return;
    end

    [Gm, Pm, ~, Wcp] = margin(loop_bp);
    if isempty(Pm) || isnan(Pm) || isinf(Pm)
        Pm = 180; % no gain crossover -> excellent margin
    end
    if isempty(Gm) || isnan(Gm) || isinf(Gm)
        GM_dB = 40; % sentinel "very good" gain margin
    else
        GM_dB = 20*log10(Gm);
    end
    diagnostics.PM = Pm;
    diagnostics.GM_dB = GM_dB;
    diagnostics.Wcp = Wcp;

    cost = -Pm + boxPenalty; % minimize -PM = maximize PM
catch
    cost = 1e8;
    feasible = false;
end
end
