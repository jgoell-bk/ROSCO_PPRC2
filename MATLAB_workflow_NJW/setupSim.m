% Nathan Wei / Claude Code
% Proportional-resonant control for FOWTs
%
% Set up one OpenFAST PPPR simulation case end-to-end:
%   1) Compute the beta->phi and Tau_g->omega gains via calculateGains.m
%      (zeta=1, stability_margin=0.02 by default; kp_kr_ratio used as
%      kp_kr_ratio_max).
%   2) Write the reference-waveform parameters and gains to the UMaineSemi
%      DISCON.IN file.
%   3) Set the ElastoDyn PtfmPitch/RotSpeed initial conditions to the
%      reference waveforms' values at t=0, and ensure blade root bending
%      moment (RootMyb1) and tip deflection (OoPDefl1) are in the OutList.
%   4) Set the .fst TMax to Npds periods of the forcing frequency.
%   5) Set the InflowFile's steady wind speed (HWindSpeed) to Uinf.
%
% Inputs:
%   phi_amp      = amplitude of the platform-pitch sinusoidal reference (deg)
%   phi_phase    = phase offset of the platform-pitch reference (deg)
%   phi_offset   = mean offset of the platform-pitch reference (deg)
%   omega_amp    = amplitude of the rotor-speed sinusoidal reference (rad/s)
%   omega_phase  = phase offset of the rotor-speed reference (deg)
%   TSR0         = design tip-speed ratio (-); sets omega_offset = TSR0*Uinf/R
%   beta0        = design/operating-point blade pitch angle (deg)
%   Uinf         = mean wind speed (m/s); deployed as the steady InflowFile speed
%   freq         = forcing frequency of both the phi and omega references (rad/s);
%                  same value used for PPPR_freq_phi and PPPR_freq_omega
%   Npds         = number of periods of freq to simulate (-); sets
%                  the .fst TMax = Npds*2*pi/freq (s)
%   kp_kr_ratio  = ceiling on |kp/kr| for the beta->phi gain search
%                  (kp_kr_ratio_max); also the fixed kp/kr ratio for the
%                  Tau_gen-omega loop (-) -- if unsure, set this to 10
% Created: 3 September 2026

% Default values and simulation guidelines (NJW):
%   freq         < 0.2187 rad/s (this is the resonant frequency of the
%                  system so forcing above this may not work well)
%   Uinf         < 10 m/s (rated wind speed is 10.75 m/s, with larger
%                  phi_amp or omega_amp the generator torque may saturate.
%                  Values of 8 m/s or lower are safer.
%   Npds         > 5-10 AFTER steady-state is reached. Use
%                  analyze_OpenFAST_phaseavg.m to determine settling time.
%   phi_offset   = 3 deg. This seems to be a decent equilibrium, but if it
%                  looks like the controller is fighting to maintain it,
%                  you may want to change this. One sign of that is if the
%                  mean blade pitch angle is significantly off the value
%                  you specified for beta0 as an input to this function.
%   kp_kr_ratio  = 10. This is safe; 5 should work as well. If kr gets too
%                  high, the system can go unstable.

function setupSim(phi_amp, phi_phase, phi_offset, omega_amp, omega_phase, ...
    TSR0, beta0, Uinf, freq, Npds, kp_kr_ratio)

% Hard-coded parameters for gain setting
zeta = 1; % damping ratio for generator-torque controller
stability_margin = 0.02; % minimum real component of any pole in the beta-phi controller system (>0 = LHP)

%% 1) Gains
params = calculateGains(TSR0, beta0, Uinf, freq, zeta, stability_margin, kp_kr_ratio);
omega_offset = params.omega0; % TSR0*Uinf/R

%% Locate the simulation's input files
repoRoot = fileparts(fileparts(mfilename('fullpath'))); % up from MATLAB_workflow_NJW to the repo root
caseDir  = fullfile(repoRoot, 'Examples', 'Test_Cases', 'IEA-15-240-RWT', 'IEA-15-240-RWT-UMaineSemi');
disconFile = fullfile(caseDir, 'IEA-15-240-RWT-UMaineSemi_DISCON.IN');
elastoFile = fullfile(caseDir, 'IEA-15-240-RWT-UMaineSemi_ElastoDyn.dat');
fstFile    = fullfile(caseDir, 'IEA-15-240-RWT-UMaineSemi.fst');
inflowFile = fullfile(fileparts(caseDir), 'IEA-15-240-RWT', 'IEA-15-240-RWT_InflowFile.dat');

%% 2) Reference-waveform parameters and gains -> DISCON.IN
writeDisconFile(disconFile, params, phi_amp, phi_phase, phi_offset, ...
    omega_amp, omega_phase, omega_offset, freq);

%% 3) Initial conditions from the reference waveforms at t=0 -> ElastoDyn
phi0_deg = phi_offset - phi_amp*sind(phi_phase); % phi_ref(0), deg
omega0_rpm = (omega_offset - omega_amp*sind(omega_phase)) * 30/pi; % omega_ref(0), rad/s -> rpm
writeOpenFASTValue(elastoFile, 'PtfmPitch', sprintf('%.6g', phi0_deg));
writeOpenFASTValue(elastoFile, 'RotSpeed', sprintf('%.6g', omega0_rpm));

%% 3.5) Ensure blade root bending moment and tip deflection are reported
% (matches the channel names analyze_OpenFAST_phaseavg.m looks for:
% RootMxb1 -- the FLAPWISE (out-of-plane, thrust/pitch-driven) root
% bending moment, not RootMyb1's edgewise (gravity-driven) one -- is its
% top-priority RootMoment candidate, OoPDefl1 is its only TipDefl
% candidate). No-op if already present in the OutList.
ensureOutListEntry(elastoFile, 'RootMxb1');
ensureOutListEntry(elastoFile, 'OoPDefl1');

%% 4) Simulation duration -> .fst
TMax = Npds*(2*pi/freq);
writeOpenFASTValue(fstFile, 'TMax', sprintf('%.6g', TMax));

%% 5) Inflow wind speed -> InflowFile
writeOpenFASTValue(inflowFile, 'HWindSpeed', sprintf('%.6g', Uinf));

end

%% ------------------------------------------------------------------
function writeDisconFile(disconFile, params, phi_amp, phi_phase, phi_offset, ...
    omega_amp, omega_phase, omega_offset, freq)
% Rewrite the PPPR reference-waveform and gain lines in disconFile,
% leaving everything else -- including the trailing "! ParamName -
% description" comments -- untouched.
lines = strsplit(fileread(disconFile), '\n', 'CollapseDelimiters', false);

lines = setDisconValue(lines, 'PPPR_amp_phi', sprintf('%.6g', phi_amp));
lines = setDisconValue(lines, 'PPPR_freq_phi', sprintf('%.6g', freq));
lines = setDisconValue(lines, 'PPPR_offset_phi', sprintf('%.6g', phi_offset));
lines = setDisconValue(lines, 'Phi_phaseoffset', sprintf('%.6g', phi_phase));
lines = setDisconValue(lines, 'PPPR_amp_omega', sprintf('%.6g', omega_amp));
lines = setDisconValue(lines, 'PPPR_freq_omega', sprintf('%.6g', freq));
lines = setDisconValue(lines, 'PPPR_offset_omega', sprintf('%.6g', omega_offset));
lines = setDisconValue(lines, 'Omega_phaseoffset', sprintf('%.6g', omega_phase));
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
    error('setupSim:notFound', 'Could not find "%s" in %s', paramName, lines{1});
end
commentStart = strfind(lines{idx}, '!');
lines{idx} = [newValueStr, '    ', lines{idx}(commentStart(1):end)];
end

%% ------------------------------------------------------------------
function ensureOutListEntry(elastoFile, channelName)
% Add "channelName" (quoted) to ElastoDyn's OutList if not already
% present, inserting it just before the "END of input file" sentinel
% line. No-op (and no file write) if the channel is already listed.
lines = strsplit(fileread(elastoFile), '\n', 'CollapseDelimiters', false);
quoted = ['"', channelName, '"'];
if any(strcmp(strtrim(lines), quoted))
    return;
end
endIdx = find(contains(lines, 'END of input file'), 1);
if isempty(endIdx)
    error('setupSim:noOutListEnd', 'Could not find the OutList END sentinel in %s', elastoFile);
end
lines = [lines(1:endIdx-1), {quoted}, lines(endIdx:end)];
fid = fopen(elastoFile, 'w');
fwrite(fid, strjoin(lines, '\n'));
fclose(fid);
end

%% ------------------------------------------------------------------
function writeOpenFASTValue(file, paramName, newValueStr)
% Replace the value portion of the line whose SECOND whitespace-separated
% token is paramName (the "Value ParamName - description" convention used
% by ElastoDyn/.fst/InflowFile primary input files, as opposed to
% DISCON.IN's "Value ! ParamName" convention), preserving everything from
% paramName onward -- spacing, description, etc. -- exactly.
lines = strsplit(fileread(file), '\n', 'CollapseDelimiters', false);
pattern = ['^\s*\S+(\s+', regexptranslate('escape', paramName), '\s.*)$'];
found = false;
for idx = 1:numel(lines)
    tok = regexp(lines{idx}, pattern, 'tokens', 'once');
    if ~isempty(tok)
        lines{idx} = [newValueStr, tok{1}];
        found = true;
        break;
    end
end
if ~found
    error('setupSim:notFound', 'Could not find parameter "%s" in %s', paramName, file);
end
fid = fopen(file, 'w');
fwrite(fid, strjoin(lines, '\n'));
fclose(fid);
end
