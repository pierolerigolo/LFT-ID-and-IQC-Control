%% d_sswepcheck  -  gain between the designed controller and the one that ran
%
%   For every sine-sweep run (open-loop sweep + closed-loop sweep):
%     1. from the recorded command: u_rec = g_cmd * K(-y)   (utils.implgain, lsim)
%     2. from the sensitivity only (does not use the recorded command):
%          S_meas = y_cl/y_ol = 1/(1 + g_S K G)  ->  g_S = (1/S_meas - 1)/(K G)
%        with G = GFRF(iOP) of the identification (median over a band)
%     3. plot S_meas against 1/(1 + K G) (as designed) and 1/(1 + g K G)
%   g = 1 is expected. g_cmd and g_S agree when the plant seen by the controller
%   is the identified one; g_S alone also catches a scaling AFTER the logged
%   command (actuator side) or BEFORE the logged acceleration (sensor side).
%
clear; close all; clc;

%% Config
expdir  = fullfile('data', 'Experiments');
contdir = fullfile('data', 'finalrepportcont');
% {experiment file, controller file, controller variable, index (Inf = end)}
runs = {
    fullfile('Old Experiments', 'gfrf0_ellips_80hz_sswep_order20.mat'),                           'controller_ellips_sus_7nm_order20_82hz_-30db_nomperf.mat', 'K_itr', 5
    'gfrf0_K_final_''controller_ellips_sus_7nm_order20_82hz_-30db_nomperf.mat',                    'controller_ellips_sus_7nm_order20_82hz_-30db_nomperf.mat', 'K_itr', Inf
    'gfrf0_K_final_''controller_ellips_sus_7nm_order20_82hz_-30db_nomperf_stbgain5.6.mat',         'controller_ellips_sus_7nm_order20_82hz_-30db_nomperf.mat', 'K_itr', Inf
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
res = table('Size', [nr 6], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'run', 'g_cmd', 'lag', 'residual_pct', 'g_S', 'g_S_phase_deg'});
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
    S  = utils.sweepS(expfile);
    w  = S.Frequency;   f = w/2/pi;
    sel = f >= f_band(1) & f <= f_band(2) & abs(f - f_notch) > 3 & ...
          w >= min(w_id) & w <= max(w_id);
    Sm = squeeze(S.ResponseData);
    KG = squeeze(freqresp(K, w(sel))).*squeeze(interp(G, w(sel)).ResponseData);
    gS = (1./Sm(sel) - 1)./KG;
    gS_med = median(abs(gS));   gS_ph = rad2deg(median(angle(gS)));

    res(i, :) = {string(name), r.gain, r.lag, 100*r.residual, gS_med, gS_ph};

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
fprintf(['g_cmd: recorded command / designed K (lsim, lag 0 and residual ~0 = pure gain).\n' ...
         'g_S  : gain explaining the measured sensitivity with G = GFRF(%d); ~g_cmd if the\n' ...
         '       plant seen in closed loop is the identified one, 1 when the loop is as designed.\n'], iOP);

%% Recorded vs simulated command (last run, zoom)
figure('Name', 'command: recorded vs lsim(K)')
[~, kp] = max(abs(r.u_th));   zm = r.t(kp) + [-0.05 0.05];
plot(r.t, r.u_rec, 'b', r.t, r.gain*r.u_th, 'r--'); grid on; xlim(zm)
legend('recorded command', sprintf('%.3f x lsim(K, -y)', r.gain))
xlabel('time (s)'); ylabel('controller output'); title(sprintf('%s', name), 'Interpreter', 'none')
