% PR gain sweep - phi controller stability frontier
% For each (kp, kr) combination, runs the v2 simulation and records:
%   - Did it run to completion (no NaN)?
%   - Did any actuator saturate?
%   - How much did rotor speed drop?
%   - How much did the platform drift?
% Outputs a printed table and a 2D heatmap of "instability index".
%
% Keeps the omega controller off (m_Tg=0) to isolate the phi path.
% Uses the same dynamics, base controllers, HP filter, and equilibrium offsets
% as PR_control_nonlinear_sim_v2.m.

clear; clc; close all;

%% Sweep grid
kp_values = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0];
kr_values = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2];
m_Tg      = 0.1;   % omega controller on at 10% authority (was 0 for phi-only)

n_kp = length(kp_values);
n_kr = length(kr_values);

%% Result storage
status        = strings(n_kp, n_kr);   % 'stable' | 'drift' | 'unstable'
min_TSR       = zeros(n_kp, n_kr);
spd_drop_pct  = zeros(n_kp, n_kr);     % percentage drop in rotor speed
ptfm_pp       = zeros(n_kp, n_kr);     % platform pitch peak-to-peak [deg]
ptfm_drift    = zeros(n_kp, n_kr);     % final mean platform tilt [deg]
pitch_max     = zeros(n_kp, n_kr);
pitch_saturated  = false(n_kp, n_kr);
torque_saturated = false(n_kp, n_kr);

fprintf('Sweeping %d x %d = %d cases...\n', n_kp, n_kr, n_kp*n_kr);
tic;
for i = 1:n_kp
    for j = 1:n_kr
        kp = kp_values(i);
        kr = kr_values(j);
        try
            r = run_PR_sim(kp, kr, m_Tg);
            min_TSR(i,j)         = r.min_TSR;
            spd_drop_pct(i,j)    = r.spd_drop_pct;
            ptfm_pp(i,j)         = r.ptfm_pp;
            ptfm_drift(i,j)      = r.ptfm_drift;
            pitch_max(i,j)       = r.pitch_max;
            pitch_saturated(i,j) = r.pitch_saturated;
            torque_saturated(i,j)= r.torque_saturated;

            % Classify
            if r.has_nan || r.min_TSR < 4 || r.pitch_saturated || r.torque_saturated
                status(i,j) = "unstable";
            elseif abs(r.spd_drop_pct) > 10 || abs(r.ptfm_drift) > 5
                status(i,j) = "drift";
            else
                status(i,j) = "stable";
            end
        catch ME
            warning('run failed kp=%g kr=%g: %s', kp, kr, ME.message);
            status(i,j) = "error";
        end
    end
end
fprintf('Sweep done in %.1f s\n\n', toc);

%% Print table
fprintf('STATUS (rows=kp, cols=kr):\n');
fprintf('%6s |', 'kp\kr');
for j = 1:n_kr, fprintf(' %9g', kr_values(j)); end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 7 + 10*n_kr));
for i = 1:n_kp
    fprintf('%6g |', kp_values(i));
    for j = 1:n_kr, fprintf(' %9s', status(i,j)); end
    fprintf('\n');
end

fprintf('\nROTOR SPEED DROP (%%):\n');
fprintf('%6s |', 'kp\kr');
for j = 1:n_kr, fprintf(' %9g', kr_values(j)); end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 7 + 10*n_kr));
for i = 1:n_kp
    fprintf('%6g |', kp_values(i));
    for j = 1:n_kr, fprintf(' %9.1f', spd_drop_pct(i,j)); end
    fprintf('\n');
end

fprintf('\nMIN TSR (failure ~3):\n');
fprintf('%6s |', 'kp\kr');
for j = 1:n_kr, fprintf(' %9g', kr_values(j)); end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 7 + 10*n_kr));
for i = 1:n_kp
    fprintf('%6g |', kp_values(i));
    for j = 1:n_kr, fprintf(' %9.2f', min_TSR(i,j)); end
    fprintf('\n');
end

%% Stability heatmap (3-color: stable=green, drift=yellow, unstable=red)
status_num = zeros(n_kp, n_kr);
status_num(status == "stable")   = 2;
status_num(status == "drift")    = 1;
status_num(status == "unstable") = 0;
status_num(status == "error")    = -1;

figure('Position', [100 100 1100 500]);

subplot(1,2,1);
imagesc(status_num);
colormap([0.5 0.5 0.5; 0.9 0.2 0.2; 0.95 0.85 0.2; 0.2 0.8 0.3]);  % error/unstable/drift/stable
clim([-1 2]);
set(gca, 'XTick', 1:n_kr, 'XTickLabel', kr_values);
set(gca, 'YTick', 1:n_kp, 'YTickLabel', kp_values);
xlabel('k_r'); ylabel('k_p');
title('Stability frontier (green=stable, yellow=drift, red=unstable)');
axis equal tight;
% Annotate cells with status text
for i = 1:n_kp
    for j = 1:n_kr
        text(j, i, status(i,j), 'HorizontalAlignment', 'center', ...
             'FontSize', 9, 'Color', 'k', 'FontWeight', 'bold');
    end
end

subplot(1,2,2);
imagesc(spd_drop_pct);
colormap(gca, flipud(parula));
colorbar;
set(gca, 'XTick', 1:n_kr, 'XTickLabel', kr_values);
set(gca, 'YTick', 1:n_kp, 'YTickLabel', kp_values);
xlabel('k_r'); ylabel('k_p');
title('Rotor speed drop (%)');
axis equal tight;
for i = 1:n_kp
    for j = 1:n_kr
        text(j, i, sprintf('%.1f', spd_drop_pct(i,j)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'k');
    end
end

sgtitle(sprintf('phi PR gain sweep (m_{Tg}=%g, UMaineSemi, U=10.74 m/s)', m_Tg));

%% Save results to .mat for later analysis
save('PR_gain_sweep_results.mat', 'kp_values', 'kr_values', 'm_Tg', ...
     'status', 'min_TSR', 'spd_drop_pct', 'ptfm_pp', 'ptfm_drift', ...
     'pitch_max', 'pitch_saturated', 'torque_saturated');
fprintf('\nResults saved to PR_gain_sweep_results.mat\n');


%% ===================================================================
%% Local function: single PR simulation, returns metrics
%% ===================================================================
function r = run_PR_sim(kp2, kr, m_Tg)
    % --- Setup (mirrors PR_control_nonlinear_sim_v2.m) ---
    TSR0  = 8.5;  beta0 = 0;  pds = 5;  dt = 1e-2;
    Uinf  = 10.74;
    Ng = 1;  Jr = 3.525e8;  ht = 150;  R = 120;  rho = 1.2;
    omega_rated = TSR0*Uinf/R;
    Tau_rated   = 1.979e7;

    % Platform (UMaineSemi)
    Jt = 1.25e10;  Dt = 2.7e8;  Kt = 5.674e8;

    % Base controllers
    PC_Kp = 1.143;  PC_Ki = 0.1196;
    PC_MinPit = 0;  PC_MaxPit = pi/2;
    VS_Rgn2K = 3.28052e7;  VS_MinTq = 0;  VS_MaxTq = Tau_rated;

    % PR settings
    freq = 0.213;
    t = (0 : dt : pds/freq*2*pi)';
    kTg_PR_ratio = 10;
    phi_amp = 5;  omega_amp = 0.01;  omega_phase = 90;

    % Aero lookup + omega-gain computation
    data = load('IEA15MW_Cp_Ct_Cq.mat');
    omega0 = TSR0*Uinf/R;
    [dCq_dbetas, dCq_dTSRs] = gradient(data.Cq, data.angles, data.TSRs);
    [dCt_dbetas, dCt_dTSRs] = gradient(data.Ct, data.angles, data.TSRs);
    dTa_dUs = 1/2*rho*pi*R^3*(2*Uinf*data.Cq - R*omega0*dCq_dTSRs);
    dTa_dU0 = interp2(data.betas, data.lambdas, dTa_dUs, beta0, TSR0, 'linear');
    k_Tg = m_Tg*ht/Ng*dTa_dU0;
    kr_Tg = k_Tg/kTg_PR_ratio;

    % Init arrays
    t = [0; 0; t];
    omega       = zeros(size(t));
    omega_ref   = omega_amp*sin(freq*t - omega_phase*pi/180);
    omega_error = zeros(size(t));
    phi_ref     = deg2rad(phi_amp)*sin(freq*t);
    phi         = zeros(size(t));
    phi_error   = zeros(size(t));
    phi_dot     = zeros(size(t));
    beta_PC     = ones(size(t))*beta0*pi/180;
    beta_PR     = zeros(size(t));
    beta        = ones(size(t))*beta0*pi/180;
    TSR         = ones(size(t))*TSR0;
    Tau_VS      = zeros(size(t));
    Tau_PR      = zeros(size(t));
    Tau_g       = zeros(size(t));
    Urel        = zeros(size(t));
    Uprime      = zeros(size(t));

    PC_Int = 0;
    phi_HP_corner = 0.02;  phi_LP = 0;  phi_HP = 0;
    alpha_HP = exp(-dt * phi_HP_corner);

    b0 = 4 + freq^2*dt^2;
    b1 = -8 + 2*freq^2*dt^2;
    b2 = 4 + freq^2*dt^2;
    a0 = b0*kp2 + 2*dt*kr;
    a1 = b1*kp2;
    a2 = b2*kp2 - 2*dt*kr;
    a0g = b0*k_Tg + 2*dt*kr_Tg;
    a1g = b1*k_Tg;
    a2g = b2*k_Tg - 2*dt*kr_Tg;

    beta_eq = beta0 * pi/180;
    Tau_g_eq = max(VS_MinTq, min(VS_MaxTq, VS_Rgn2K * omega0^2));

    for ti = 3 : length(t)
        x = [0; omega(ti-1); phi(ti-1); phi_dot(ti-1)];
        u = [beta(ti-1) - beta_eq; Tau_g(ti-1) - Tau_g_eq; Uprime(ti-1)];
        Urel(ti-1) = Uinf - ht*cos(phi(ti-1))*phi_dot(ti-1) + Uprime(ti-1);
        TSR(ti-1)  = (omega0 + omega(ti-1))*R/Urel(ti-1);

        dCq_dTSR  = interp2(data.betas, data.lambdas, dCq_dTSRs, beta(ti-1)*180/pi, TSR(ti-1));
        dCt_dTSR  = interp2(data.betas, data.lambdas, dCt_dTSRs, beta(ti-1)*180/pi, TSR(ti-1));
        dCq_dbeta = interp2(data.betas, data.lambdas, dCq_dbetas, beta(ti-1)*180/pi, TSR(ti-1));
        dCt_dbeta = interp2(data.betas, data.lambdas, dCt_dbetas, beta(ti-1)*180/pi, TSR(ti-1));
        dTa_domega = 1/2*rho*Urel(ti-1)*pi*R^4*dCq_dTSR;
        dFa_domega = 1/2*rho*Urel(ti-1)*pi*R^3*dCt_dTSR;
        dTa_dbeta  = 1/2*rho*Urel(ti-1)^2*pi*R^3*dCq_dbeta * 180/pi;
        dFa_dbeta  = 1/2*rho*Urel(ti-1)^2*pi*R^2*dCt_dbeta * 180/pi;
        dTa_dU = 1/2*rho*pi*R^3*(2*Urel(ti-1) ...
            * interp2(data.betas, data.lambdas, data.Cq, beta(ti-1)*180/pi, TSR(ti-1)) ...
            - R*omega(ti-1)*dCq_dTSR);
        dFa_dU = 1/2*rho*pi*R^2*(2*Urel(ti-1) ...
            * interp2(data.betas, data.lambdas, data.Ct, beta(ti-1)*180/pi, TSR(ti-1)) ...
            - R*omega(ti-1)*dCt_dTSR);

        A0 = [0 1 0 0;
            0 Ng/Jr*dTa_domega 0 -ht*Ng/Jr*dTa_dU;
            0 0 0 1;
            0 ht/Jt*dFa_domega -Kt/Jt -1/Jt*(Dt+ht^2*dFa_dU)];
        Bc = [0 0 0;
            Ng/Jr*dTa_dbeta -Ng^2/Jr Ng/Jr*dTa_dU;
            0 0 0;
            ht/Jt*dFa_dbeta 0 ht/Jt*dFa_dU];
        x_dot  = A0*x + Bc*u;
        x_next = x + x_dot*dt;
        omega(ti)   = x_next(2);
        phi(ti)     = x_next(3);
        phi_dot(ti) = x_next(4);

        phi_LP = alpha_HP * phi_LP + (1 - alpha_HP) * phi(ti);
        phi_HP = phi(ti) - phi_LP;
        phi_error(ti)   = phi_HP - phi_ref(ti);
        omega_error(ti) = omega(ti) - omega_ref(ti);

        omega_total = omega0 + omega(ti);
        Tau_VS(ti) = max(VS_MinTq, min(VS_MaxTq, VS_Rgn2K * omega_total^2));

        PC_SpdErr = omega_total - omega_rated;
        PC_Int = max(PC_MinPit, min(PC_MaxPit, PC_Int + PC_Ki * PC_SpdErr * dt));
        beta_PC(ti) = max(PC_MinPit, min(PC_MaxPit, PC_Kp * PC_SpdErr + PC_Int));

        beta_PR(ti) = 1/b0 * ( a0*phi_error(ti) + a1*phi_error(ti-1) + a2*phi_error(ti-2) ...
            - b1*beta_PR(ti-1) - b2*beta_PR(ti-2) );
        Tau_PR(ti) = 1/b0 * ( a0g*omega_error(ti) + a1g*omega_error(ti-1) + a2g*omega_error(ti-2) ...
            - b1*Tau_PR(ti-1) - b2*Tau_PR(ti-2) );

        beta(ti)  = max(PC_MinPit, min(PC_MaxPit, beta_PC(ti) + beta_PR(ti)));
        Tau_g(ti) = max(VS_MinTq, min(VS_MaxTq, Tau_VS(ti) + Tau_PR(ti)));

        % Abort if state diverges (NaN check)
        if ~isfinite(omega(ti)) || ~isfinite(phi(ti))
            break;
        end
    end

    % Remove padding
    omega(1:2)=[]; phi(1:2)=[]; beta(1:2)=[]; Tau_g(1:2)=[]; TSR(1:2)=[];

    % Compute metrics
    r.has_nan          = any(~isfinite(omega)) || any(~isfinite(phi));
    r.min_TSR          = min(TSR(isfinite(TSR)));
    if isempty(r.min_TSR), r.min_TSR = NaN; end
    spd_initial        = omega0*60/(2*pi);
    spd_final          = (omega0 + omega(end))*60/(2*pi);
    r.spd_drop_pct     = 100 * (spd_initial - spd_final) / spd_initial;
    r.ptfm_pp          = rad2deg(max(phi) - min(phi));
    r.ptfm_drift       = rad2deg(phi(end));    % final platform tilt = drift
    r.pitch_max        = rad2deg(max(beta));
    r.pitch_saturated  = any(beta >= PC_MaxPit - 1e-6);
    r.torque_saturated = any(Tau_g >= VS_MaxTq - 1e-3);
end
