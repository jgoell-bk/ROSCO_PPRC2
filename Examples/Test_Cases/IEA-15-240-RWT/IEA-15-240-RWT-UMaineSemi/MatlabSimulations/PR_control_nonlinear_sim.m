% Nathan Wei
% Proportional-resonant control for FOWTs - waveform optimization
% Specify waveform (phi and omega amplitude/freq/phase)
% Outputs system dynamics, torque/thrust/power mean/amplitude
% Created: 8 April 2026 (modified from PR_control_test.m)

%% Inputs

% Operating setpoint
TSR0    = 8.5; % optimal TSR is 8.5, design TSR is 9.0
beta0   = 0; % (deg) optimal blade pitch is -1, design pitch is 0
Uinf    = 10; % wind speed (m/s) - rated is 10.75 m/s
freq    = 0.2; % desired forcing frequency (rad/s)
fz      = freq/10; % desired PI controller zero frequency (rad/s)
zeta    = 0.7; % desired damping ratio [-]
R       = 120; % turbine radius (m)

% Reference waveforms
phi_amp         = 1; % deg
phi_offset      = -2; % deg
omega_amp       = 0.01; % rad/s
omega_phase     = 90; % deg
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
Dt = 1e7; % platform damping coefficient (e.g. hydrodynamics)
Kt = 1e8; % platform restoring coefficient (e.g. from mooring lines)
% Dt and Kt are set based on platform pitch simulation data; analysis code
% can be found in IEA15MW_steady.m.

% Set up simulation parameters
dt = 1e-2;
% w0 = sqrt(Kt/Jt); % natural frequency of platform (rad/s)
t = (0 : dt : pds/freq*2*pi)';
interpMethod = 'linear'; % 2D interpolation method for quasi-steady aerodynamics

kp_kr_ratio = 10; % kp/kr for PR controller gains
m_Tg = 1; % in range [0, 1]

%% Sensitivity coefficients from IEA 15 MW steady data
data = load('IEA15MW_Cp_Ct_Cq.mat');
omega0 = TSR0*Uinf/R;
[dCp_dbetas, dCp_dTSRs] = gradient(data.Cp, data.angles, data.TSRs);
[dCq_dbetas, dCq_dTSRs] = gradient(data.Cq, data.angles, data.TSRs);
[dCt_dbetas, dCt_dTSRs] = gradient(data.Ct, data.angles, data.TSRs);
dTa_dUs = 1/2*rho*pi*R^3*(2*Uinf*data.Cq - R*omega0*dCq_dTSRs); % sensitivity of aero torque to wind speed (kg m/s)
dTa_dU0 = interp2(data.betas, data.lambdas, dTa_dUs, beta0, TSR0, interpMethod); % value at setpoint
dCp_dTSR = interp2(data.betas, data.lambdas, dCp_dTSRs, beta0, TSR0, interpMethod); % value at setpoint
dCp_dbeta = interp2(data.betas, data.lambdas, dCp_dbetas, beta0, TSR0, interpMethod); % value at setpoint
Cp0 = interp2(data.betas, data.lambdas, data.Cp, beta0, TSR0, interpMethod);

%% Calculate gains using modified formulae from Abbas et al. (2022)
A       = 1/2*rho*pi*R^4*Uinf^2 / (Jr*TSR0^2*Uinf) * (dCp_dTSR*TSR0 - Cp0);
B_beta  = Ng/(2*Jr*TSR0^2) * rho*pi*R^3*Uinf^2 * (dCp_dbeta*TSR0);
B_Tg    = -Ng^2 / Jr;
kp      = -1/(2*pi*B_beta)*(2*zeta*freq + A); % for beta-phi PR controller
kr      = kp/kp_kr_ratio; % for beta-phi PR controller
% k_Tg  = m_Tg*ht/Ng*dTa_dU0; % Eqn. 29, from Fischer (2013) and Stockhouse et al. (2021)
kp_Tg   = -1/B_Tg*(2*zeta*freq + A);
kr_Tg   = kp_Tg/kp_kr_ratio;

%% Simulate controller in discrete time
nSkip = 3; % 3rd order difference equation = pad first 3 time points
t = [zeros([nSkip, 1]); t]; % first two entries are padding
omega = ones(size(t))*omega0;
omega_ref = omega_amp*sin(freq*t - omega_phase*pi/180) + omega_offset;
omega_error = zeros(size(t));
phi_ref = deg2rad(phi_amp)*sin(freq*t) + deg2rad(phi_offset);
phi = zeros(size(t)); % in rad
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
fz_Tg = freq / 10;
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

    % Compute local sensitivity derivatives
    dCq_dTSR = interp2(data.betas, data.lambdas, dCq_dTSRs, beta(ti-1), TSR(ti-1), interpMethod);
    dCt_dTSR = interp2(data.betas, data.lambdas, dCt_dTSRs, beta(ti-1), TSR(ti-1), interpMethod);
    dCq_dbeta = interp2(data.betas, data.lambdas, dCq_dbetas, beta(ti-1), TSR(ti-1), interpMethod);
    dCt_dbeta = interp2(data.betas, data.lambdas, dCt_dbetas, beta(ti-1), TSR(ti-1), interpMethod);
    % Dimensional versions
    dTa_domega = 1/2*rho*Urel(ti-1)*pi*R^4*dCq_dTSR; % sensitivity of aero torque to rotation rate (N s/rad)
    dFa_domega = 1/2*rho*Urel(ti-1)*pi*R^3*dCt_dTSR; % sensitivity of thrust force to rotation rate (N s / rad*m)
    dTa_dbeta = 1/2*rho*Urel(ti-1)^2*pi*R^3*dCq_dbeta * 180/pi; % sensitivity of aero torque to blade pitch angle (data are in deg, convert to rad)
    dFa_dbeta = 1/2*rho*Urel(ti-1)^2*pi*R^2*dCt_dbeta * 180/pi; % sensitivity of thrust force to blade pitch angle (data are in deg, convert to rad)
    dTa_dU = 1/2*rho*pi*R^3*(2*Urel(ti-1) ...
        * interp2(data.betas, data.lambdas, data.Cq, beta(ti-1), TSR(ti-1)) ...
        - R*omega(ti-1)*dCq_dTSR); % sensitivity of aero torque to wind speed (kg m/s)
    dFa_dU = 1/2*rho*pi*R^2*(2*Urel(ti-1) ...
        * interp2(data.betas, data.lambdas, data.Ct, beta(ti-1), TSR(ti-1)) ...
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
controller_Tg = -(kp_Tg*(1+fz/s) + kr_Tg*s/(s^2+freq^2));
loop_bp = minreal(sys_beta_to_phi_OL * controller_beta);
sys_bp = minreal(loop_bp / (1 + loop_bp)); % y/r = PC / (1+PC), Rowley 2.16a
loop_to = minreal(sys_Tg_to_omega_OL * controller_Tg);
sys_to = minreal(loop_to / (1 + loop_to)); % y/r = PC / (1+PC), Rowley 2.16a
loop_bo = minreal(sys_beta_to_omega_OL * controller_beta);
sys_bo = minreal(loop_bo / (1 + loop_bo)); % y/r = PC / (1+PC), Rowley 2.16a
loop_tp = minreal(sys_Tg_to_phi_OL * controller_Tg);
sys_tp = minreal(loop_tp / (1 + loop_tp)); % y/r = PC / (1+PC), Rowley 2.16a