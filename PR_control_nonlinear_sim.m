% Nathan Wei
% Proportional-resonant control for FOWTs - waveform optimization
% Specify waveform (phi and omega amplitude/freq/phase)
% Outputs system dynamics, torque/thrust/power mean/amplitude
% Created: 8 April 2026 (modified from PR_control_test.m)

%% Inputs

% Operating setpoint
TSR0    = 8.5; % optimal TSR is 8.5, design TSR is 9.0
beta0   = 0; % (deg) optimal blade pitch is -1, design pitch is 0
Uinf    = 8; % wind speed (m/s) - rated is 10.75 m/s
freq    = 0.1; % desired forcing frequency (rad/s)
% fz      = freq/10; % desired PI controller zero frequency (rad/s)
zeta    = 0.7; % desired damping ratio [-]
zeta_Tg = 1;
R       = 120; % turbine radius (m)

% Reference waveforms
phi_amp         = 1; % deg
phi_offset      = 3; % deg
omega_amp       = 0.01; % rad/s
omega_phase     = 0; % deg
omega_offset    = 1*TSR0*Uinf/R;

pds     = 10; % number of periods
uStar   = 0.1; % surge-velocity amplitude / Uinf

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

% Set up simulation parameters
dt = 1e-2;
% w0 = sqrt(Kt/Jt); % natural frequency of platform (rad/s)
t = (0 : dt : pds/freq*2*pi)';
interpMethod = 'linear'; % 2D interpolation method for quasi-steady aerodynamics

kp_kr_ratio = 5; % kp/kr for PR controller gains
m_Tg = 1; % in range [0, 1]

%% Sensitivity coefficients from IEA 15 MW steady data
data = load('IEA15MW_Cp_Ct_Cq.mat');
omega0 = TSR0*Uinf/R;
[dCp_dbetas, dCp_dTSRs] = gradient(data.Cp, data.angles, data.TSRs);
[dCq_dbetas, dCq_dTSRs] = gradient(data.Cq, data.angles, data.TSRs);
[dCt_dbetas, dCt_dTSRs] = gradient(data.Ct, data.angles, data.TSRs);
dTa_dUs = 1/2*rho*pi*R^3*(2*Uinf*data.Cq - R*omega0*dCq_dTSRs); % sensitivity of aero torque to wind speed (kg m/s)
dTa_dU0 = interp2(data.betas, data.lambdas, dTa_dUs, beta0, TSR0, interpMethod); % value at setpoint
dCp_dTSR0 = interp2(data.betas, data.lambdas, dCp_dTSRs, beta0, TSR0, interpMethod); % value at setpoint
dCp_dbeta0 = interp2(data.betas, data.lambdas, dCp_dbetas, beta0, TSR0, interpMethod) * 180/pi; % value at setpoint (data are in deg, convert to rad -- matches dTa_dbeta/dFa_dbeta convention below)
Cp0 = interp2(data.betas, data.lambdas, data.Cp, beta0, TSR0, interpMethod);

dCt_dTSR0 = interp2(data.betas, data.lambdas, dCt_dTSRs, beta0, TSR0, interpMethod);
dCt_dbeta0 = interp2(data.betas, data.lambdas, dCt_dbetas, beta0, TSR0, interpMethod) * 180/pi; % value at setpoint (data are in deg, convert to rad -- matches dTa_dbeta/dFa_dbeta convention below)
Ct0 = interp2(data.betas, data.lambdas, data.Ct, beta0, TSR0, interpMethod);

% Characterization
phi_eq  = ht*(1/2*rho*pi*R^2*Uinf^2*Ct0)/Kt;
freqd_max = (Dt + ht^2*rho*pi*R^2*Uinf*(Ct0-(TSR0/2)*dCt_dTSR0)) / (2*Jt*zeta);
freqd = freqd_max/2;
freqd_Tg = freqd_max/2;

%% Calculate gains using modified formulae from Abbas et al. (2022)
A       = 1/2*rho*pi*R^4*Uinf^2 / (Jr*TSR0^2*Uinf) * (dCp_dTSR0*TSR0 - Cp0);
B_beta  = Ng/(2*Jr*TSR0^2) * rho*pi*R^3*Uinf^2 * (dCp_dbeta0*TSR0);
B_Tg    = -Ng^2 / Jr;
% kp      = -1/(B_beta)*(2*zeta*omega_des + A); % for beta-phi PR controller
kp      = (Jt*freqd^2*(1-4*zeta^2) + 2*zeta*freqd*(Dt + rho*pi*R^2*ht^2*Uinf ...
          * (Ct0 - (TSR0/2)*dCt_dTSR0)) - Kt) / (1/2*rho*pi*R^2*ht*Uinf^2*dCt_dbeta0); % dCt_dbeta < 0, dCt_dTSR > 0
ki      = -((Dt + rho*pi*R^2*ht^2*Uinf*(Ct0 - (TSR0/2)*dCt_dTSR0))*freqd^2 ...
          - 2*Jt*zeta*freqd^3) / (1/2*rho*pi*R^2*ht*Uinf^2*dCt_dbeta0);
fz      = ki/kp;
kr      = kp/kp_kr_ratio; % for beta-phi PR controller
% k_Tg  = m_Tg*ht/Ng*dTa_dU0; % Eqn. 29, from Fischer (2013) and Stockhouse et al. (2021)
kp_Tg   = -1/B_Tg*(2*zeta_Tg*freqd_Tg + A);
kr_Tg   = kp_Tg/kp_kr_ratio;
fz_Tg   = -freqd_Tg^2 / B_Tg / kp_Tg;

%% Simulate controller in discrete time
nSkip = 3; % 3rd order difference equation = pad first 3 time points
t = [zeros([nSkip, 1]); t]; % first two entries are padding
omega = ones(size(t))*omega0;
omega_ref = omega_amp*sin(freq*t - omega_phase*pi/180) + omega_offset;
omega_error = zeros(size(t));
phi_ref = deg2rad(phi_amp)*sin(freq*t) + deg2rad(phi_offset);
phi = ones(size(t))*deg2rad(phi_offset); % in rad
phi_error = zeros(size(t));
phi_dot = zeros(size(t));
beta = ones(size(t))*beta0; % now in deg
TSR = ones(size(t))*TSR0;
Tau_g = zeros(size(t));
Urel = zeros(size(t)); % relative velocity (including phi_dot contribution)
Uprime = zeros(size(t)); % not implemented yet

% PR controller parameters (copied from Joeri Frederik's pull request)
b0 = 4 + freq^2*dt^2; % b coeffs are attached to output vars
b1 = -8 + 2*freq^2*dt^2;
b2 = 4 + freq^2*dt^2;
% a0 = b0*kp2 + 2*dt*kr; % a coeffs are attached to input vars
% a1 = b1*kp2;
% a2 = b2*kp2 - 2*dt*kr;
n0 = kp*(1 + fz*dt/2)*b0 + 2*kr*dt;
n1 = kp*(b1-b0) + kp*fz*dt/2*(b0+b1) - 2*kr*dt;
n2 = kp*(b0-b1) + kp*fz*dt/2*(b0+b1) - 2*kr*dt;
n3 = kp*(fz*dt/2-1)*b0 + 2*kr*dt;

% PR controller parameters (copied from Joeri Frederik's pull request)
% fz_Tg = freq / 10;
% a0g = b0*k_Tg + 2*dt*kr_Tg; % a coeffs are attached to input vars
% a1g = b1*k_Tg;
% a2g = b2*k_Tg - 2*dt*kr_Tg;
n0g = kp_Tg*(1 + fz_Tg*dt/2)*b0 + 2*kr_Tg*dt;
n1g = kp_Tg*(b1-b0) + kp_Tg*fz_Tg*dt/2*(b0+b1) - 2*kr_Tg*dt;
n2g = kp_Tg*(b0-b1) + kp_Tg*fz_Tg*dt/2*(b0+b1) - 2*kr_Tg*dt;
n3g = kp_Tg*(fz_Tg*dt/2-1)*b0 + 2*kr_Tg*dt;

% Run simulation (discrete time)
for ti = nSkip+1 : length(t)
    % System dynamics: calculate current state based on previous state and
    % previous control action
    x = [0; omega(ti-1); phi(ti-1); phi_dot(ti-1)]; % previous state
    u = [beta(ti-1); Tau_g(ti-1); Uprime(ti-1)]; % previous control action
    Urel(ti-1) = Uinf - ht*cos(phi(ti-1))*phi_dot(ti-1) + Uprime(ti-1);
    % TSR(ti-1) = (omega0 + omega(ti-1))*R/Urel(ti-1);
    TSR(ti-1) = omega(ti-1)*R/Urel(ti-1);

    % Clamp the aero-table query point to the lookup table's domain. The
    % startup transient (cold-start phi_dot spike feeding into Urel, hence
    % TSR) briefly pushes TSR outside the table range within the first few
    % timesteps; interp2 returns NaN outside its domain with no
    % extrapolation method specified, which then propagates through the
    % rest of the simulation. Clamping here only affects the aero
    % sensitivity lookup, not the physical beta/TSR state themselves.
    beta_q = min(max(beta(ti-1), min(data.betas(:))), max(data.betas(:)));
    TSR_q  = min(max(TSR(ti-1),  min(data.lambdas(:))), max(data.lambdas(:)));

    % Compute local sensitivity derivatives
    dCq_dTSR = interp2(data.betas, data.lambdas, dCq_dTSRs, beta_q, TSR_q, interpMethod);
    dCt_dTSR = interp2(data.betas, data.lambdas, dCt_dTSRs, beta_q, TSR_q, interpMethod);
    dCq_dbeta = interp2(data.betas, data.lambdas, dCq_dbetas, beta_q, TSR_q, interpMethod);
    dCt_dbeta = interp2(data.betas, data.lambdas, dCt_dbetas, beta_q, TSR_q, interpMethod);
    % Dimensional versions
    dTa_domega = 1/2*rho*Urel(ti-1)*pi*R^4*dCq_dTSR; % sensitivity of aero torque to rotation rate (N s/rad)
    dFa_domega = 1/2*rho*Urel(ti-1)*pi*R^3*dCt_dTSR; % sensitivity of thrust force to rotation rate (N s / rad*m)
    dTa_dbeta = 1/2*rho*Urel(ti-1)^2*pi*R^3*dCq_dbeta * 180/pi; % sensitivity of aero torque to blade pitch angle (data are in deg, convert to rad)
    dFa_dbeta = 1/2*rho*Urel(ti-1)^2*pi*R^2*dCt_dbeta * 180/pi; % sensitivity of thrust force to blade pitch angle (data are in deg, convert to rad)
    dTa_dU = 1/2*rho*pi*R^3*(2*Urel(ti-1) ...
        * interp2(data.betas, data.lambdas, data.Cq, beta_q, TSR_q) ...
        - R*omega(ti-1)*dCq_dTSR); % sensitivity of aero torque to wind speed (kg m/s)
    dFa_dU = 1/2*rho*pi*R^2*(2*Urel(ti-1) ...
        * interp2(data.betas, data.lambdas, data.Ct, beta_q, TSR_q) ...
        - R*omega(ti-1)*dCt_dTSR); % sensitivity of thrust force to wind speed (kg/s)

    % System definition
    A0 = [0 1 0 0;
        0 Ng/Jr*dTa_domega 0 -ht*Ng/Jr*dTa_dU;
        0 0 0 1;
        0 ht/Jt*dFa_domega -Kt/Jt -1/Jt*(Dt+ht^2*dFa_dU)]; % dynamics
    Bc = [0 0 0;
        Ng/Jr*dTa_dbeta -Ng^2/Jr Ng/Jr*dTa_dU;
        0 0 0;
    ht/Jt*dFa_dbeta 0 ht/Jt*dFa_dU]; % control contributions to dynamics
    x_dot = A0*x + Bc*u; % dynamics
    x_next = x + x_dot*dt; % forward Euler (can use different scheme if needed)
    omega(ti) = x_next(2); % current state
    phi(ti) = x_next(3);
    phi_dot(ti) = x_next(4);
    phi_error(ti) = phi(ti) - phi_ref(ti);
    omega_error(ti) = omega(ti) - omega_ref(ti);
    % Feedback control: calculate current control action based on current state
    % PR control (using biquad transform of PR controller)
    % http://cn.imperix.com/doc/implementation/proportional-resonant-controller.html
    % https://en.wikipedia.org/wiki/Bilinear_transform
    % Here, input is phi-phi_ref and output is beta
    % beta(ti) = 1/b0 * (a0*phi_error(ti) + a1*phi_error(ti-1) + a2*phi_error(ti-2) ...
    %     - b1*beta(ti-1) - b2*beta(ti-2)); %  - (beta(ti-1) + freq/100*(dt/2) * (phi_error(ti) - phi_error(ti-1)));
    % beta(ti) = beta(ti); % can add OL beta if desired here
    % Tau_g(ti) = 1/b0 * (a0g*omega_error(ti) + a1g*omega_error(ti-1) + a2g*omega_error(ti-2) ...
    %     - b1*Tau_g(ti-1) - b2*Tau_g(ti-2)); %  - (Tau_g(ti-1) + freq/100*(dt/2) * (omega_error(ti) - omega_error(ti-1)));
    beta(ti) = 1/b0 * (n0*phi_error(ti) + n1*phi_error(ti-1) + n2*phi_error(ti-2) ...
        + n3*phi_error(ti-3) - (b1-b0)*beta(ti-1) + (b1-b0)*beta(ti-2) + b0*beta(ti-3));
    Tau_g(ti) = 1/b0 * (n0g*omega_error(ti) + n1g*omega_error(ti-1) + n2g*omega_error(ti-2) ...
        + n3g*omega_error(ti-3) - (b1-b0)*Tau_g(ti-1) + (b1-b0)*Tau_g(ti-2) + b0*Tau_g(ti-3));
end

% Remove padding (first Nskip elements)
t(1:nSkip) = [];
omega(1:nSkip) = [];
omega_ref(1:nSkip) = [];
omega_error(1:nSkip) = [];
phi(1:nSkip) = [];
phi_dot(1:nSkip) = [];
beta(1:nSkip) = [];
phi_ref(1:nSkip) = [];
phi_error(1:nSkip) = [];
Tau_g(1:nSkip) = [];
Urel(1:nSkip) = [];
Uprime(1:nSkip) = [];
TSR(1:nSkip) = [];

%% Plot results
colors = orderedcolors('gem'); % matlab default colors
figure;
hold on;
plot(t*freq/(2*pi), rad2deg(phi), 'c', 'LineWidth', 2);
plot(t*freq/(2*pi), rad2deg(phi_ref), 'b--', 'LineWidth', 2);
plot(t*freq/(2*pi), beta, 'Color', colors(5,:), 'LineWidth', 2);
xlabel('Time, $tf$', 'interpreter', 'latex');
ylabel('Angle (degrees)', 'interpreter', 'latex');
yyaxis right;
plot(t*freq/(2*pi), omega, 'Color', colors(7,:), 'LineWidth', 2);
ax = gca;
ax.YAxis(2).Color = colors(7,:);
hold on;
plot(t*freq/(2*pi), omega_ref, 'r--', 'LineWidth', 2);
ylabel('Rotation Rate, $\omega$ (rad/s)', 'interpreter', 'latex')
legend({'$\phi$', '$\phi_{ref}$', '$\beta$', '$\omega$', '$\omega_{ref}$'}, 'interpreter', 'latex');

%% Lastly, analyze transfer functions (from PR_control_test.m)
% Compute local sensitivity derivatives
dCq_dTSR = interp2(data.betas, data.lambdas, dCq_dTSRs, beta0, TSR0, interpMethod);
dCt_dTSR = interp2(data.betas, data.lambdas, dCt_dTSRs, beta0, TSR0, interpMethod);
dCq_dbeta = interp2(data.betas, data.lambdas, dCq_dbetas, beta0, TSR0, interpMethod);
dCt_dbeta = interp2(data.betas, data.lambdas, dCt_dbetas, beta0, TSR0, interpMethod);
% Dimensional versions
dTa_domega = 1/2*rho*Uinf*pi*R^4*dCq_dTSR; % sensitivity of aero torque to rotation rate (N s/rad)
dFa_domega = 1/2*rho*Uinf*pi*R^3*dCt_dTSR; % sensitivity of thrust force to rotation rate (N s / rad*m)
dTa_dbeta = 1/2*rho*Uinf^2*pi*R^3*dCq_dbeta * 180/pi; % sensitivity of aero torque to blade pitch angle (data are in deg, convert to rad)
dFa_dbeta = 1/2*rho*Uinf^2*pi*R^2*dCt_dbeta * 180/pi; % sensitivity of thrust force to blade pitch angle (data are in deg, convert to rad)
dTa_dU = 1/2*rho*pi*R^3*(2*Uinf ...
    * interp2(data.betas, data.lambdas, data.Cq, beta0, TSR0) ...
    - R*omega0*dCq_dTSR); % sensitivity of aero torque to wind speed (kg m/s)
dFa_dU = 1/2*rho*pi*R^2*(2*Uinf ...
    * interp2(data.betas, data.lambdas, data.Ct, beta0, TSR0) ...
    - R*omega0*dCt_dTSR); % sensitivity of thrust force to wind speed (kg/s)

% NEW: TFs manually copied from PR_control_analysis.m (symbolic solve)
s = tf('s');
sys_beta_to_phi_OL = (ht*(Ng*dFa_domega*dTa_dbeta - Ng*dFa_dbeta*dTa_domega ...
    + Jr*dFa_dbeta*s))/(Jr*Kt*s - Kt*Ng*dTa_domega + Dt*Jr*s^2 + Jr*Jt*s^3 ...
    - Jt*Ng*dTa_domega*s^2 + Jr*dFa_dU*ht^2*s^2 - Dt*Ng*dTa_domega*s ...
    - Ng*dFa_dU*dTa_domega*ht^2*s + Ng*dFa_domega*dTa_dU*ht^2*s); % no control
sys_Tg_to_phi_OL = -(Ng^2*dFa_domega*ht)/(Jr*Kt*s - Kt*Ng*dTa_domega + Dt*Jr*s^2 ...
    + Jr*Jt*s^3 - Jt*Ng*dTa_domega*s^2 + Jr*dFa_dU*ht^2*s^2 - Dt*Ng*dTa_domega*s ...
    - Ng*dFa_dU*dTa_domega*ht^2*s + Ng*dFa_domega*dTa_dU*ht^2*s); % no control
sys_beta_to_omega_OL = (Jt*Ng*dTa_dbeta*s^2 + Ng*(Dt*dTa_dbeta + dFa_dU*dTa_dbeta*ht^2 ...
    - dFa_dbeta*dTa_dU*ht^2)*s + Kt*Ng*dTa_dbeta) ...
    / (Jr*Jt*s^3 + (Jr*dFa_dU*ht^2 + Dt*Jr - Jt*Ng*dTa_domega)*s^2 ...
    + (Jr*Kt - Dt*Ng*dTa_domega - Ng*dFa_dU*dTa_domega*ht^2 ...
    + Ng*dFa_domega*dTa_dU*ht^2)*s - Kt*Ng*dTa_domega);
sys_Tg_to_omega_OL = (- Jt*Ng^2*s^2 + (- dFa_dU*ht^2 - Dt)*Ng^2*s - Kt*Ng^2) ...
    / (Jr*Jt*s^3 + (Jr*dFa_dU*ht^2 + Dt*Jr - Jt*Ng*dTa_domega)*s^2 ...
    + (Jr*Kt - Dt*Ng*dTa_domega - Ng*dFa_dU*dTa_domega*ht^2 ...
    + Ng*dFa_domega*dTa_dU*ht^2)*s - Kt*Ng*dTa_domega);

%% Define controllers, open-loop systems, and closed-loop systems
controller_beta = -(kp*(1+fz/s) + kr*s/(s^2+freq^2));
controller_Tg = -(kp_Tg*(1+fz_Tg/s) + kr_Tg*s/(s^2+freq^2));
loop_bp = minreal(sys_beta_to_phi_OL * controller_beta);
sys_bp = minreal(loop_bp / (1 + loop_bp)); % y/r = PC / (1+PC), Rowley 2.16a
loop_to = minreal(sys_Tg_to_omega_OL * controller_Tg);
sys_to = minreal(loop_to / (1 + loop_to)); % y/r = PC / (1+PC), Rowley 2.16a
loop_bo = minreal(sys_beta_to_omega_OL * controller_beta);
sys_bo = minreal(loop_bo / (1 + loop_bo)); % y/r = PC / (1+PC), Rowley 2.16a
loop_tp = minreal(sys_Tg_to_phi_OL * controller_Tg);
sys_tp = minreal(loop_tp / (1 + loop_tp)); % y/r = PC / (1+PC), Rowley 2.16a

% margin(loop_bp)

%% Root locus analysis: kp (fz, kp/kr fixed)
% With fz and kp_kr_ratio both held fixed, controller_beta is EXACTLY
% proportional to kp:
%   controller_beta = -kp*[(1+fz/s) + (1/kp_kr_ratio)*s/(s^2+freq^2)]
%                    = kp * C_shape_kp(s)
% so sweeping kp is a standard root-locus gain sweep of
% sys_beta_to_phi_OL*C_shape_kp(s) closed in unity negative feedback.
C_shape_kp = -(1 + fz/s + (1/kp_kr_ratio)*s/(s^2+freq^2));
G_rlocus_kp = minreal(sys_beta_to_phi_OL * C_shape_kp);
figure;
rlocus(G_rlocus_kp);
title(sprintf('Root locus vs. k_p  (f_z=%.4g, k_p/k_r=%.4g fixed)', fz, kp_kr_ratio));
hold on;
plot(real(pole(sys_bp)), imag(pole(sys_bp)), 'kd', 'MarkerSize', 9, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('current design point k_p=%.4g', kp));
legend('show', 'Location', 'best');

% Numerical cross-check: the closed-loop poles at the CURRENT kp should be
% identical whether computed directly (controller_beta/loop_bp/sys_bp) or
% as the root locus evaluated at gain=kp -- confirms the K*G(s)
% reformulation above is algebraically correct, not just plausible-looking.
poles_direct_kp = pole(sys_bp);
poles_rlocus_kp = rlocus(G_rlocus_kp, kp);
fprintf('\n=== Root locus (k_p) cross-check at k_p=%.5f ===\n', kp);
fprintf('  max(Re), direct (sys_bp)      = %.6f\n', max(real(poles_direct_kp)));
fprintf('  max(Re), via rlocus(G,k_p)    = %.6f\n', max(real(poles_rlocus_kp)));

%% Root locus analysis: fz (kp, kp/kr fixed)
% With kp and kp/kr (hence kr) held fixed, fz enters the closed loop as an
% OUTER loop closed around the already kp+kr-compensated system. Writing
% controller_beta = A_fz(s) + fz*B_fz(s), with
%   A_fz(s) = -(kp + kr*s/(s^2+freq^2))   [proportional+resonant part only, fz=0]
%   B_fz(s) = -kp/s                       [coefficient multiplying fz]
% the closed-loop characteristic equation 1+sys_beta_to_phi_OL*controller_beta=0,
% after dividing through by (1+sys_beta_to_phi_OL*A_fz(s)) -- i.e. closing
% the proportional+resonant loop first -- becomes the standard root-locus
% form 1 + fz*S0_fz(s)*B_fz(s) = 0, a gain sweep in fz.
A_fz = -(kp + kr*s/(s^2+freq^2));
loop_fz0 = minreal(sys_beta_to_phi_OL * A_fz);
S0_fz = minreal(sys_beta_to_phi_OL / (1 + loop_fz0));
B_fz = -kp/s;
G_rlocus_fz = minreal(S0_fz * B_fz);
figure;
rlocus(G_rlocus_fz);
title(sprintf('Root locus vs. f_z  (k_p=%.4g, k_p/k_r=%.4g fixed)', kp, kp_kr_ratio));
hold on;
plot(real(pole(sys_bp)), imag(pole(sys_bp)), 'kd', 'MarkerSize', 9, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('current design point f_z=%.4g', fz));
legend('show', 'Location', 'best');

% Same cross-check as above, this time sweeping fz through the inner-loop
% reformulation.
poles_direct_fz = pole(sys_bp);
poles_rlocus_fz = rlocus(G_rlocus_fz, fz);
fprintf('\n=== Root locus (f_z) cross-check at f_z=%.5f ===\n', fz);
fprintf('  max(Re), direct (sys_bp)      = %.6f\n', max(real(poles_direct_fz)));
fprintf('  max(Re), via rlocus(G,f_z)    = %.6f\n', max(real(poles_rlocus_fz)));

%% Critical kr / kr_Tg via numerical pole sweep, then set gains to 10% of critical
% For each loop, C(s) = Cpi(s) + kr*Cres(s), with Cpi(s) = kp*(1+fz/s) and
% Cres(s) = s/(s^2+freq^2) (both carrying the same leading minus sign as
% controller_beta/controller_Tg above, to match this script's controller
% sign convention). find_critical_kr sweeps kr and finds where the
% closed-loop poles first cross into the right half-plane.

Cpi_phi   = -kp*(1 + fz/s);
Cres_phi  = -s/(s^2 + freq^2);
kr_crit    = find_critical_kr(Cpi_phi,   Cres_phi,   sys_beta_to_phi_OL);

Cpi_omega  = -kp_Tg*(1 + fz_Tg/s);
Cres_omega = -s/(s^2 + freq^2);
% kr_Tg_crit = find_critical_kr(Cpi_omega, Cres_omega, sys_Tg_to_omega_OL);

kr_2    = 0.1*kr_crit;
% kr_Tg_2 = 0.1*kr_Tg_crit;

fprintf('\n=== Critical kr (numerical pole sweep) ===\n');
fprintf('phi loop:   kr_crit    = %.6g   ->   kr    = 0.1*kr_crit    = %.6g\n', kr_crit, kr_2);
% fprintf('omega loop: kr_Tg_crit = %.6g   ->   kr_Tg = 0.1*kr_Tg_crit = %.6g\n', kr_Tg_crit, kr_Tg_2);

%% ------------------------------------------------------------------
function kr_crit = find_critical_kr(Cpi, Cres, G)
% Numerically find the smallest kr>0 at which the closed-loop poles of
% G*(Cpi + kr*Cres), closed in unity negative feedback, first cross into
% the right half-plane (max(Re(pole)) crosses zero).
%
% NOTE: near-cancelling pole/zero pairs (e.g. at very small kr, where
% Cres's contribution to the loop is tiny) can leave minreal with leftover
% poles whose real part is pure floating-point noise (~1e-10 to 1e-12),
% which a bare ">0" check would misread as instability. tol filters that
% out; genuine crossings here have real parts many orders larger (~1e-4+).
tol = 1e-6;

    % Exponential search to bracket the order of magnitude of kr_crit
    kr_hi = 1e-6;
    maxRe = -1;
    while maxRe < tol
        kr_hi = kr_hi*10;
        loop_i = minreal(G*(Cpi + kr_hi*Cres));
        maxRe = max(real(pole(minreal(loop_i/(1+loop_i)))));
        if kr_hi > 1e15
            error('find_critical_kr: no imaginary-axis crossing found up to kr=1e15');
        end
    end

    % Linear sweep within the bracket, then interpolate to max(Re)=0
    kr_sweep = linspace(0, kr_hi, 500);
    maxReVec = zeros(size(kr_sweep));
    for i = 1:length(kr_sweep)
        loop_i = minreal(G*(Cpi + kr_sweep(i)*Cres));
        maxReVec(i) = max(real(pole(minreal(loop_i/(1+loop_i)))));
    end
    idx = find(maxReVec > tol, 1, 'first');
    kr_crit = interp1(maxReVec(idx-1:idx), kr_sweep(idx-1:idx), 0);
end