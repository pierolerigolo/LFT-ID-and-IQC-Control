%% d_gaincheck_W2  -  W2 controller (1110, order 18): gain on the setup, sensitivity, Nyquist
%
%   Runs at 0 kg with the W2 controller (etaf 0.36): sweeps at build gain 1 and 4.8,
%   80 Hz sine at build gain 1.
%     g_cmd : recorded command / lsim(K, -y)  -> factor inside the controller block (0.1)
%     g_S   : loop / designed loop from the measured S, (1/S - 1)/(K G_old)  (sweeps)
%     att   : attenuation of the 80 Hz sine (utils.sineatt)
%   Old id: GFRF_cell{iOP} of idfile (designed loop K G_old, built with build gain g_ref).
%   New id: sweep sent through x_command with the controller model; it contains the
%   actuator-side factor (~2), so the loop at build gain g is g_block*g*K*G_new.
%
clear; close all; clc;

%% Config
contfn  = 'LFT_GFRF_GdmGLTIW2opt_1110_etfe_80hz_fmin10_hyb_controller_order18_80hz_-34db_nperfoptiIQC_notconv_etai0.35_etaf0.36_multi_freq';
expdir  = fullfile('data', 'Experiments', [contfn, '_K_final']);
runs    = {'gftf0_ssweep_stbg1.mat', 1;  'gftf0_ssweep_stbg4.8.mat', 4.8};   % sweeps {file, build gain}
sinefn  = 'gfrf0_80hz_stbg1.mat';                                        % 80 Hz sine, build gain 1
idfile  = fullfile('data', 'GFRF_iddata_1110_etfe.mat');                 % old id
newidfile = fullfile('data', 'MIVIDA 261118 sswp new simulink', 'gfrf0_ssweep_id.mat');
iOP     = 1;          % 0 kg
g_ref   = 4.8;        % build gain giving the designed loop (old id)
g_block = 0.1;        % factor inside the controller block
f_newid = [20 300];   % band where the new id is reliable (Hz)
f_band  = [20 200];   % band for g_S (Hz)

bo = bodeoptions;     bo.FreqUnits = 'Hz';  bo.XLim = [10 500];  bo.Grid = 'on';
bo.Title.String = 'W2 controller, 0 kg: sensitivity';
no = nyquistoptions;  no.FreqUnits = 'Hz';  no.ShowFullContour = 'off';  no.Grid = 'on';
no.XLim = [-3 1];     no.YLim = [-2 2];     % zoom around -1

%% Controller and identifications
K   = getfield(load(fullfile('data', 'controller', [contfn, '.mat']), 'K_final'), 'K_final');
idd = load(idfile, 'GFRF_cell', 'realmass');
wf  = 2*pi*linspace(10, 500, 5000).';                       % fine grid for the plots
G_old = interp(idd.GFRF_cell{iOP}, wf(wf <= max(idd.GFRF_cell{iOP}.Frequency)));

[y, ~, u] = utils.loadexp(newidfile);                       % new id: x_command -> Payload_acc1
seg = find(sqrt(movmean(u.^2, 500)) > 0.05*max(sqrt(movmean(u.^2, 500))));
d_new = detrend(iddata(y(seg(1):seg(end)), u(seg(1):seg(end)), K.Ts));   % as in a_MultimodelFRFID
G_new = etfe(d_new, 1000, 2000);                           % lag 1000, 2000 points (0.25 Hz), as in a_MultimodelFRFID
G_new = frd(G_new.ResponseData, G_new.Frequency, K.Ts);
G_new = fselect(G_new, 2*pi*f_newid(1), 2*pi*f_newid(2));   % reliable band

%% Gain of each run
S = cell(size(runs, 1), 1);   g_cmd = zeros(size(runs, 1), 1);   g_S = g_cmd;
for i = 1:size(runs, 1)
    fn = fullfile(expdir, runs{i, 1});
    g_cmd(i) = utils.implgain(fn, K).gain;
    S{i} = utils.sweepS(fn);
    w  = S{i}.Frequency(S{i}.Frequency >= 2*pi*f_band(1) & S{i}.Frequency <= 2*pi*f_band(2));
    gs = (1./squeeze(interp(S{i}, w).ResponseData) - 1)./squeeze(freqresp(K*interp(idd.GFRF_cell{iOP}, w), w));
    g_S(i) = median(abs(gs));
end
disp(table(string(runs(:, 1)), [runs{:, 2}].', g_cmd, g_S, [runs{:, 2}].'./g_S, ...
    'VariableNames', {'run', 'build_gain', 'g_cmd', 'g_S', 'build_gain_for_design'}))
fprintf('80 Hz sine (build gain 1): attenuation %.1f dB, predicted %.1f dB (old id, K/%.3g)\n', ...
    utils.sineatt(fullfile(expdir, sinefn)), ...
    mag2db(abs(1/(1 + squeeze(freqresp(K*interp(idd.GFRF_cell{iOP}, 2*pi*80)/g_ref, 2*pi*80))))), g_ref);

%% Sensitivity: measured vs old id and new id
figure
bodemag(S{1}, 'b', S{2}, 'r', ...
    feedback(1, K*G_old/g_ref), 'k--', feedback(1, K*G_old), 'k', ...
    feedback(1, g_block*K*G_new), 'm--', feedback(1, g_block*g_ref*K*G_new), 'm', bo)
legend(sprintf('measured, build gain %g', runs{1, 2}), sprintf('measured, build gain %g', runs{2, 2}), ...
    sprintf('old id, 1/(1+K G/%g)  (build gain 1)', g_ref), 'old id, 1/(1+K G)  (build gain 4.8)', ...
    sprintf('new id, 1/(1+%g K G)  (build gain 1)', g_block), sprintf('new id, 1/(1+%g K G)  (build gain %g)', g_block*g_ref, g_ref), ...
    'Location', 'southwest')
exportgraphics(gcf, fullfile('plots', 'd_gaincheck_W2_bode.pdf'), 'ContentType', 'vector')

%% Nyquist: old id (gain 1, every operating point), new id, new id with the gain correction
L_old = cellfun(@(Gi) K*interp(Gi, wf(wf <= max(Gi.Frequency))), idd.GFRF_cell, 'UniformOutput', false);
L_new = {K*G_new, g_block*g_ref*K*G_new};
figure; nyquistplot(L_old{:}, no); title('L = K G, old id (gain 1)')
legend(compose('%.1f kg', idd.realmass), 'Location', 'best')
exportgraphics(gcf, fullfile('plots', 'd_gaincheck_W2_nyq_oldid.pdf'), 'ContentType', 'vector')
figure; nyquistplot(L_new{:}, no); title(sprintf('L = K G, new id (%g-%g Hz)', f_newid))
legend('new id, gain 1', sprintf('new id, %g K G (build gain %g)', g_block*g_ref, g_ref), 'Location', 'best')
exportgraphics(gcf, fullfile('plots', 'd_gaincheck_W2_nyq_newid.pdf'), 'ContentType', 'vector')

% closest point to -1
L_all = [L_old(:); L_new(:)];
names = [compose('old id, %.1f kg', idd.realmass(:)); "new id, gain 1"; sprintf("new id, build gain %g", g_ref)];
[m1pL, kmin] = cellfun(@(L) min(abs(1 + squeeze(L.ResponseData))), L_all);
f_min = cellfun(@(L, k) L.Frequency(k)/2/pi, L_all, num2cell(kmin));
disp(table(names, f_min, m1pL, -mag2db(m1pL), 'VariableNames', {'loop', 'f_Hz', 'min_1pL', 'peak_S_dB'}))
