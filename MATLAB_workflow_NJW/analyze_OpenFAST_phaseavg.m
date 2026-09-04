%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nathan Wei
% Phase-average an OpenFAST PPPR (platform-pitch/rotor-speed resonant
% control) simulation once it has reached periodic steady state.
%
% Pipeline:
%   1) Import the .out time series (channel names read from the file's
%      own header, as in quick_analyze_OpenFAST.m, so column mapping
%      can't silently drift).
%   2) Read the PPPR reference-waveform parameters (amplitude, frequency,
%      offset, phase for both platform pitch and rotor speed) AND the
%      control gains (PPPR_CntrGains_phi/omega, PPPR_Mode) from the
%      companion *_DISCON.IN file in the same folder.
%   3) Compute phi/omega reference-tracking error vs time and find the
%      time at which each has settled to periodic steady state (per-
%      period RMS tracking error stays within settle_tol of its long-run
%      value -- a settling-time criterion, analogous to a control-systems
%      2%/5% settling band).
%   4)-5) Phase-average Power, GenTq, RotSpeed, BldPitch, blade root
%      bending moment, blade tip deflection, PtfmPitch, and PtfmSurge
%      over the steady-state span, using phase_average.m.
%   6) Cycle-to-cycle standard deviations from the same phase-averaging
%      call.
%   7) Amplitude, phase (via get_phase.m), and time-average, all computed
%      from the single-period phase-averaged waveform -- not the raw
%      steady-state segment, since t_steady:end isn't guaranteed to span
%      an integer number of periods and a partial trailing period would
%      skew a raw mean.
%   8) Save (results.mat alongside the .out file) and return as a struct.
%
% Inputs:
%   filename    = OpenFAST .out file name (with or without extension)
%   folder      = folder containing filename and its *_DISCON.IN file
%   settle_tol  = (optional) relative tolerance for the settling-time
%                 detector (default 0.05, i.e. 5%)
% Output:
%   results, struct containing:
%       filename, folder      = as given
%       ref                   = struct of PPPR reference-waveform params
%                                read from DISCON.IN
%       gains                 = struct of control gains read from
%                                DISCON.IN: PPPR_Mode, kp, kr, fz (phi
%                                loop), kp_Tg, kr_Tg, fz_Tg (omega loop)
%       Uinf_mean             = mean Wind1VelX (m/s) over the
%                                steady-state span (raw mean, not
%                                phase-averaged -- wind speed isn't part
%                                of the periodic response set)
%       dt                    = sample time (s)
%       T_period              = forcing period used for phase-averaging,
%                                2*pi/ref.freq_phi (s)
%       N_phaseavg            = samples per period
%       t_phaseavg            = (N_phaseavg x 1) time within one period (s)
%       t_steady              = simulation time at which steady state is
%                                judged reached (s)
%       n_periods             = number of complete periods spanned by
%                                t_steady:end (i.e. actually averaged over)
%       is_stable             = 1 if the phi/omega tracking error settled
%                                (or is at least not still growing) by the
%                                end of the run; 0 if either shows a
%                                non-finite value or a per-period RMS
%                                that's still increasing in the final
%                                quarter of the run relative to the
%                                quarter before it (i.e. appears to be
%                                diverging/unstable)
%       settle_tol            = tolerance used above
%       units                 = struct of unit strings per quantity
%       phaseavg              = struct of phase-averaged vectors (step 5)
%       phaseavg_std          = struct of cycle-to-cycle std vectors (step 6)
%       amp, phase, mean      = structs of scalar amplitude (units),
%                                phase (deg), and time-average (units)
%                                per quantity (step 7)
%   Quantities (struct fields of phaseavg/phaseavg_std/amp/phase/mean):
%       Power (GenPwr, kW), GenTq (kN-m), RotSpeed (rad/s),
%       BldPitch (BldPitch1, deg), RootMoment (blade root bending moment,
%       kN-m), TipDefl (blade tip deflection, OoPDefl1, m),
%       PtfmPitch (deg), PtfmSurge (m)
%   A quantity whose channel isn't present in this .out file's OutList
%   (e.g. blade root bending moment, if ElastoDyn wasn't configured to
%   output it) is filled with NaN and a warning is issued -- it is not a
%   fatal error, since which channels exist depends on the run's OutList.
%
% Dependencies: get_phase.m, phase_average.m
% Created: 3 September 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = analyze_OpenFAST_phaseavg(filename, folder, settle_tol)

% Nathan: hack for local testing
if nargin == 0
    filename = 'IEA-15-240-RWT-UMaineSemi.out';
    folder = '/Users/njwei/Documents/GitHub/ROSCO_PPRC2/Examples/Test_Cases/IEA-15-240-RWT/IEA-15-240-RWT-UMaineSemi/';
end

if nargin < 3 || isempty(settle_tol)
    settle_tol = 0.05; % 5% settling-time criterion
end
if ~endsWith(filename, '.out', 'IgnoreCase', true)
    filename = [filename, '.out'];
end

if isempty(which('phase_average')) || isempty(which('get_phase'))
    addpath(fileparts(mfilename('fullpath'))); % local copies live alongside this script
end

%% 1) Import the .out data
T = read_OpenFAST_out(fullfile(folder, filename));

%% 2) Reference-waveform parameters and control gains from the companion DISCON.IN
disconFile = find_DISCON_file(folder);
[ref, gains] = read_DISCON_ref_params(disconFile);
if ~(isfinite(ref.freq_phi) && ref.freq_phi > 0 && isfinite(ref.freq_omega) && ref.freq_omega > 0)
    error('analyze_OpenFAST_phaseavg:badFreq', ...
        'PPPR_freq_phi/PPPR_freq_omega not found or non-positive in %s.', disconFile);
end

%% 3) Reference-tracking error and steady-state detection
dt = median(diff(T.Time));
phi = T.PtfmPitch; % deg
omega = T.RotSpeed*pi/30; % rad/s
phi_ref = ref.phi_offset + ref.phi_amp*sin(ref.freq_phi*T.Time - deg2rad(ref.phi_phase));
omega_ref = ref.omega_offset + ref.omega_amp*sin(ref.freq_omega*T.Time - deg2rad(ref.omega_phase));
phi_err_norm = (phi - phi_ref) / max(ref.phi_amp, eps);
omega_err_norm = (omega - omega_ref) / max(ref.omega_amp, eps);

Nper_phi = round((2*pi/ref.freq_phi)/dt);
Nper_omega = round((2*pi/ref.freq_omega)/dt);
if Nper_phi ~= Nper_omega
    warning('analyze_OpenFAST_phaseavg:mismatchedPeriods', ...
        ['PPPR_freq_phi and PPPR_freq_omega give different periods (%.4g vs %.4g s); ', ...
        'phase-averaging uses the phi period. Each channel''s own settling time is still ', ...
        'computed against its own period.'], 2*pi/ref.freq_phi, 2*pi/ref.freq_omega);
end
Nper = Nper_phi;

[t_settle_phi, diverging_phi] = find_settle_time(T.Time, phi_err_norm, Nper_phi, settle_tol);
[t_settle_omega, diverging_omega] = find_settle_time(T.Time, omega_err_norm, Nper_omega, settle_tol);
t_steady = max(t_settle_phi, t_settle_omega);
is_stable = double(~(diverging_phi || diverging_omega));

%% 4)-5) Phase-average from t_steady to the end of the simulation
ssMask = T.Time >= t_steady;
t_phaseavg = (0:Nper-1)'*dt;
n_periods = floor(sum(ssMask)/Nper); % complete periods spanned by the steady-state window
                                      % (phase_average also folds in any leftover partial
                                      % period, but only into a subset of its phase bins)

channels = struct();
channels.Power      = getChannel(T, {'GenPwr'}, 'Power');
channels.GenTq      = getChannel(T, {'GenTq'}, 'GenTq');
channels.RotSpeed   = omega; % rad/s, already computed above
channels.BldPitch   = getChannel(T, {'BldPitch1'}, 'BldPitch');
% Flapwise (out-of-plane, thrust/pitch-driven) candidates first --
% RootMxb1 is what setupSim.m adds to the OutList. Edgewise (gravity-
% driven) names are kept only as a last-resort fallback.
channels.RootMoment = getChannel(T, {'RootMxb1', 'RootMFlp1', 'RootMOoP1', ...
                                      'RootMyb1', 'RootMEdg1', 'RootMIP1'}, 'RootMoment');
channels.TipDefl    = getChannel(T, {'OoPDefl1'}, 'TipDefl');
channels.PtfmPitch  = phi; % deg, already computed above
channels.PtfmSurge  = getChannel(T, {'PtfmSurge'}, 'PtfmSurge');

% Actual wind speed for this run (an operating condition, not part of the
% periodic response set, so reported as a plain steady-state mean rather
% than phase-averaged)
windSpeed = getChannel(T, {'Wind1VelX'}, 'Wind1VelX');
Uinf_mean = mean(windSpeed(ssMask), 'omitnan');

units = struct('Power', 'kW', 'GenTq', 'kN-m', 'RotSpeed', 'rad/s', ...
    'BldPitch', 'deg', 'RootMoment', 'kN-m', 'TipDefl', 'm', ...
    'PtfmPitch', 'deg', 'PtfmSurge', 'm');

fn = fieldnames(channels);
phaseavg = struct(); phaseavg_std = struct();
amp = struct(); phase = struct(); avg = struct();
for i = 1:numel(fn)
    y = channels.(fn{i})(ssMask);
    %% 5)-6) Phase-averaged vector and cycle-to-cycle standard deviation
    [pa, pa_std] = phase_average(y, Nper, 0);
    phaseavg.(fn{i}) = pa;
    phaseavg_std.(fn{i}) = pa_std;
    %% 7) Amplitude, phase, and time-average, all from the single-period
    %% phase-averaged waveform -- NOT the raw steady-state segment, since
    %% t_steady:end isn't guaranteed to span an integer number of periods
    %% and a partial trailing period would skew a raw mean
    [phase.(fn{i}), amp.(fn{i})] = get_phase(t_phaseavg, pa);
    avg.(fn{i}) = mean(pa, 'omitnan');
end

%% 8) Assemble and save
results = struct();
results.filename = filename;
results.folder = folder;
results.ref = ref;
results.gains = gains;
results.Uinf_mean = Uinf_mean;
results.dt = dt;
results.T_period = 2*pi/ref.freq_phi;
results.N_phaseavg = Nper;
results.t_phaseavg = t_phaseavg;
results.t_steady = t_steady;
results.n_periods = n_periods;
results.is_stable = is_stable;
results.settle_tol = settle_tol;
results.units = units;
results.phaseavg = phaseavg;
results.phaseavg_std = phaseavg_std;
results.amp = amp;
results.phase = phase;
results.mean = avg;

[~, baseName] = fileparts(filename);
save(fullfile(folder, [baseName, '_phaseavg.mat']), 'results', '-v7.3');

end

%% ------------------------------------------------------------------
function T = read_OpenFAST_out(outFile)
% Import an OpenFAST tabular .out file, reading channel names from the
% file's own header (row 7) rather than a hardcoded list.
fid = fopen(outFile, 'r');
if fid < 0
    error('analyze_OpenFAST_phaseavg:fileNotFound', 'Could not open %s', outFile);
end
for k = 1:6
    fgetl(fid); % title/description lines
end
varNames = matlab.lang.makeUniqueStrings(strsplit(strtrim(fgetl(fid)), '\t'));
fclose(fid);

opts = delimitedTextImportOptions("NumVariables", numel(varNames));
opts.DataLines = [9, Inf];
opts.Delimiter = "\t";
opts.VariableNames = varNames;
opts.VariableTypes = repmat({'double'}, 1, numel(varNames));
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, varNames, "ThousandsSeparator", ",");

T = readtable(outFile, opts);
end

%% ------------------------------------------------------------------
function disconFile = find_DISCON_file(folder)
% Locate the *DISCON*.IN file in folder (case-insensitive), erroring if
% none or more than one is found.
d = dir(folder);
names = {d.name};
isMatch = ~cellfun(@isempty, regexpi(names, 'DISCON.*\.IN$', 'once'));
matches = names(isMatch);
if isempty(matches)
    error('analyze_OpenFAST_phaseavg:noDISCON', 'No *DISCON*.IN file found in %s', folder);
elseif numel(matches) > 1
    error('analyze_OpenFAST_phaseavg:ambiguousDISCON', ...
        'Multiple DISCON.IN candidates found in %s: %s', folder, strjoin(matches, ', '));
end
disconFile = fullfile(folder, matches{1});
end

%% ------------------------------------------------------------------
function [ref, gains] = read_DISCON_ref_params(disconFile)
% Parse the PPPR reference-waveform parameters AND control gains out of a
% DISCON.IN file. Matches quick_analyze_OpenFAST.m's parsing convention:
% for the (single-value) reference-waveform parameters, the value is the
% first whitespace-separated token on the line whose trailing
% "! ParamName - description" comment names it. PPPR_CntrGains_phi/omega
% carry 3 values each (kp, kr, fz), so all tokens before "!" are read.
paramMap = struct('PPPR_amp_phi', 'phi_amp', 'PPPR_freq_phi', 'freq_phi', ...
    'PPPR_offset_phi', 'phi_offset', 'PPPR_amp_omega', 'omega_amp', ...
    'PPPR_freq_omega', 'freq_omega', 'PPPR_offset_omega', 'omega_offset', ...
    'Phi_phaseoffset', 'phi_phase', 'Omega_phaseoffset', 'omega_phase');
outFields = struct2cell(paramMap);
ref = cell2struct(num2cell(NaN(numel(outFields), 1)), outFields);
gains = struct('PPPR_Mode', NaN, 'kp', NaN, 'kr', NaN, 'fz', NaN, ...
    'kp_Tg', NaN, 'kr_Tg', NaN, 'fz_Tg', NaN);

fid = fopen(disconFile, 'r');
if fid < 0
    error('analyze_OpenFAST_phaseavg:fileNotFound', 'Could not open %s', disconFile);
end
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || line(1) == '!'
        continue
    end
    tokens = strsplit(line);
    bang = find(strcmp(tokens, '!'));
    if isempty(bang) || bang >= numel(tokens)
        continue
    end
    paramName = tokens{bang + 1};
    if isfield(paramMap, paramName)
        ref.(paramMap.(paramName)) = str2double(tokens{1});
    elseif strcmp(paramName, 'PPPR_Mode')
        gains.PPPR_Mode = str2double(tokens{1});
    elseif strcmp(paramName, 'PPPR_CntrGains_phi')
        vals = str2double(tokens(1:bang-1));
        gains.kp = vals(1); gains.kr = vals(2); gains.fz = vals(3);
    elseif strcmp(paramName, 'PPPR_CntrGains_omega')
        vals = str2double(tokens(1:bang-1));
        gains.kp_Tg = vals(1); gains.kr_Tg = vals(2); gains.fz_Tg = vals(3);
    end
end
fclose(fid);

missing = outFields(structfun(@(f) isnan(ref.(f)), paramMap));
if ~isempty(missing)
    warning('analyze_OpenFAST_phaseavg:missingDISCONParams', ...
        'Could not find the following PPPR parameters in %s: %s', disconFile, strjoin(missing, ', '));
end
missingGains = fieldnames(gains);
missingGains = missingGains(structfun(@isnan, gains));
if ~isempty(missingGains)
    warning('analyze_OpenFAST_phaseavg:missingDISCONGains', ...
        'Could not find the following control gain parameters in %s: %s', disconFile, strjoin(missingGains, ', '));
end
end

%% ------------------------------------------------------------------
function col = getChannel(T, candidates, label)
% Return the first column of T matching one of candidates (in priority
% order); NaN-fill with a warning if none of them are present, since
% which channels exist depends on the run's own OutList.
for i = 1:numel(candidates)
    if ismember(candidates{i}, T.Properties.VariableNames)
        col = T.(candidates{i});
        return;
    end
end
warning('analyze_OpenFAST_phaseavg:channelNotFound', ...
    '%s: none of {%s} found in this .out file''s OutList -- filling with NaN.', ...
    label, strjoin(candidates, ', '));
col = NaN(height(T), 1);
end

%% ------------------------------------------------------------------
function [t_settle, is_diverging] = find_settle_time(t, err_norm, Nper, tol)
% Settling-time detector: compute the per-period RMS of a (normalized)
% tracking-error signal, take the long-run ("final") RMS as the median
% over the last quarter of the periods, and find the last period whose
% RMS deviates from that final value by more than tol (relative). The
% period after that is judged the start of steady state.
%
% is_diverging: true if the response contains a non-finite value (the
% integration blew up), or if the per-period RMS is still growing in the
% final quarter of the run relative to the quarter before it -- i.e. the
% error hasn't leveled off by the end of the simulation, as opposed to a
% settled or bounded-but-still-oscillating (non-divergent) response.
if any(~isfinite(err_norm))
    t_settle = t(1);
    is_diverging = true;
    return;
end
nPeriods = floor(length(err_norm)/Nper);
if nPeriods < 4
    t_settle = t(1); % too short a run to judge -- treat the whole series as steady
    is_diverging = false;
    return;
end
rmsPerPeriod = zeros(nPeriods, 1);
for k = 1:nPeriods
    idx = (k-1)*Nper + (1:Nper);
    rmsPerPeriod(k) = sqrt(mean(err_norm(idx).^2));
end
nTail = max(1, round(0.25*nPeriods));
rmsFinal = median(rmsPerPeriod(end-nTail+1:end));
tolBand = max(tol*rmsFinal, 1e-6);
lastBad = find(abs(rmsPerPeriod - rmsFinal) > tolBand, 1, 'last');
if isempty(lastBad)
    ssPeriod = 1;
else
    ssPeriod = min(lastBad + 1, nPeriods);
end
startIdx = (ssPeriod - 1)*Nper + 1;
t_settle = t(startIdx);

nQ = max(1, round(0.25*nPeriods));
if nPeriods >= 2*nQ
    rmsLastQ = mean(rmsPerPeriod(end-nQ+1:end));
    rmsPrevQ = mean(rmsPerPeriod(end-2*nQ+1:end-nQ));
    is_diverging = rmsLastQ > rmsPrevQ*(1+tol) && rmsLastQ > 1e-6;
else
    is_diverging = false; % too short to judge a trend
end
end
