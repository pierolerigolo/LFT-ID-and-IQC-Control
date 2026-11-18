%% d_sswepcheck  -  gain between the designed controller and the one that ran
%
%   For every sine-sweep run (open-loop sweep + closed-loop sweep):
%     0. g_set: gain set in the controller build (1, or e.g. 5.6)
%     1. from the recorded command: u_rec = g_cmd * K(-y)   (utils.implgain, lsim)
%        -> g_cmd/g_set = scaling inside the controller block (logged y -> logged u)
%     1b. plant from the logged command to the logged acceleration in closed loop,
%        G_cl = (Y_cl - Y_ol)/U_cl, against G = GFRF(iOP): g_plant = |G_cl/G|
%        -> scaling between the logged signals and the identification signals
%     2. from the sensitivity only (does not use the recorded command):
%          S_meas = y_cl/y_ol = 1/(1 + g_S K G)  ->  g_S = (1/S_meas - 1)/(K G)
%        with G = GFRF(iOP) of the identification (median over a band)
%     3. plot S_meas against 1/(1 + K G) (as designed) and 1/(1 + g K G)
%   Loop as designed <=> g_S = 1, and g_S ~ g_cmd * g_plant.
%   g = 1 is expected. g_cmd and g_S agree when the plant seen by the controller
%   is the identified one; g_S alone also catches a scaling AFTER the logged
%   command (actuator side) or BEFORE the logged acceleration (sensor side).
%
clear; close all; clc;

%% Config
expdir  = fullfile('data', 'Experiments');
contdir = 'data';                                  % controller files relative to this folder
soa     = fullfile('finalrepportcont', 'controller_ellips_sus_7nm_order20_82hz_-30db_nomperf.mat');
% {experiment file, controller file (in contdir), controller variable, index (Inf = end; ignored if not a cell), gain set in the build}
runs = {
    fullfile('Old Experiments', 'gfrf0_ellips_80hz_sswep_order20.mat'),                           soa, 'K_itr', 5,   1
    'gfrf0_K_final_''controller_ellips_sus_7nm_order20_82hz_-30db_nomperf.mat',                    soa, 'K_itr', Inf, 1
    'gfrf0_K_final_''controller_ellips_sus_7nm_order20_82hz_-30db_nomperf_stbgain5.6.mat',         soa, 'K_itr', Inf, 5.6
    'gfrf0_K_final_DISK_1028_etfe_controller_order18_80hz_-34db_nperfhinf_stbg_4.5.mat', fullfile('controller', 'DISK_1028_etfe_controller_order18_80hz_-34db_nperfhinf.mat'), 'K_final', Inf, 4.5
    };
idfile  = fullfile('data', 'GFRF_iddata_1028_etfe.mat');   % identification FRFs
iOP     = 1;                                               % operating point of the runs (gfrf0)
f_band  = [20 200];        % band for g_S (Hz)
f_notch = 80;              % excluded +/- 3 Hz (controller peak, S ~ 0)

%% Identification
idd = load(idfile, 'GFRF');
G   = idd.GFRF(iOP);
w_id = idd.GFRF.Frequency;

%% Runs
nr  = size(runs, 1);
res = table('Size', [nr 9], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'run', 'g_set', 'g_cmd', 'lag', 'residual_pct', 'ctrl_block', 'g_plant', 'g_S', 'g_S_phase_deg'});
col = lines(nr);
figS = figure('Name', 'sensitivity: measured vs designed vs with gain'); hold on
for i = 1:nr
    expfile = fullfile(expdir, runs{i, 1});
    cs = load(fullfile(contdir, runs{i, 2}), runs{i, 3});
    Kc = cs.(runs{i, 3});
    if iscell(Kc)
        if isinf(runs{i, 4}), K = Kc{end}; else, K = Kc{runs{i, 4}}; end
    else
        K = Kc;
    end
    [~, name] = fileparts(runs{i, 1});

    % 1. gain from the recorded command
    r = utils.implgain(expfile, K);

    % 2. gain from the sensitivity
    [S, Ufrd, Ycfrd, Yofrd] = utils.sweepS(expfile);
    w  = S.Frequency;   f = w/2/pi;
    sel = f >= f_band(1) & f <= f_band(2) & abs(f - f_notch) > 3 & ...
          w >= min(w_id) & w <= max(w_id);
    Sm = squeeze(S.ResponseData);
    KG = squeeze(freqresp(K, w(sel))).*squeeze(interp(G, w(sel)).ResponseData);
    gS = (1./Sm(sel) - 1)./KG;
    gS_med = median(abs(gS));   gS_ph = rad2deg(median(angle(gS)));

    % 1b. plant seen from the logged command to the logged acceleration
    nb  = 61;                                       % bins averaged (~1 Hz)
    Uc  = squeeze(Ufrd.ResponseData);  Yc = squeeze(Ycfrd.ResponseData);  Yo = squeeze(Yofrd.ResponseData);
    Gcl = movsum((Yc - Yo).*conj(Uc), nb)./movsum(abs(Uc).^2, nb);
    g_plant = median(abs(Gcl(sel)./squeeze(interp(G, w(sel)).ResponseData)));

    g_set = runs{i, 5};
    res(i, :) = {string(name), g_set, r.gain, r.lag, 100*r.residual, r.gain/g_set, g_plant, gS_med, gS_ph};

    % 3. sensitivity plot
    fw = f(f >= 10 & w <= max(w_id));   ww = 2*pi*fw;
    KGw = squeeze(freqresp(K, ww)).*squeeze(interp(G, ww).ResponseData);
    figure(figS)
    semilogx(fw, mag2db(abs(Sm(f >= 10 & w <= max(w_id)))), 'Color', col(i, :), ...
        'DisplayName', sprintf('%s measured', name));
    semilogx(fw, mag2db(abs(1./(1 + r.gain*KGw))), '--', 'Color', col(i, :), ...
        'DisplayName', sprintf('1/(1 + %.3g K G)', r.gain));
    fprintf('%s:\n', name); disp(r.table)
end
semilogx(fw, mag2db(abs(1./(1 + KGw))), 'k:', 'LineWidth', 1.5, 'DisplayName', 'as designed 1/(1 + K G) (last K)');
set(gca, 'XScale', 'log'); grid on; xlim([10 500]);
xlabel('Frequency (Hz)'); ylabel('|S| (dB)'); legend('Location', 'best', 'Interpreter', 'none')
title(sprintf('Sweeps: measured S vs 1/(1 + g K G), G = GFRF(%d)', iOP))

disp(res)
fprintf(['g_cmd     : recorded command / designed K (lsim; lag 0 and residual ~0 = pure gain)\n' ...
         'ctrl_block: g_cmd/g_set, scaling inside the controller block (should be 1)\n' ...
         'g_plant   : |G_cl/GFRF(%d)|, logged command -> logged acceleration vs identification (should be 1)\n' ...
         'g_S       : loop gain / designed loop from the measured S (should be 1; ~ g_cmd*g_plant)\n' ...
         'gain to set in the build for the designed loop: g_set/g_S\n'], iOP);
res.g_to_set = res.g_set./res.g_S;
disp(res(:, {'run', 'g_set', 'g_S', 'g_to_set'}))

%% Recorded vs simulated command (last run, zoom)
figure('Name', 'command: recorded vs lsim(K)')
[~, kp] = max(abs(r.u_th));   zm = r.t(kp) + [-0.05 0.05];
plot(r.t, r.u_rec, 'b', r.t, r.gain*r.u_th, 'r--'); grid on; xlim(zm)
legend('recorded command', sprintf('%.3f x lsim(K, -y)', r.gain))
xlabel('time (s)'); ylabel('controller output'); title(sprintf('%s', name), 'Interpreter', 'none')
