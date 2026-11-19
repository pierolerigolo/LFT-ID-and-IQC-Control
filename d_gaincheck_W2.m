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

% FRF frequency vector (rad/s), as in a_MultimodelFRFID: log f_lo..f_split, linear above
Ts      = 0.001;
f_lo    = 10;          % lowest frequency (Hz)
f_split = 20;          % log below, linear above (Hz)
n_log   = 50;          % points f_lo..f_split
n_lin   = 1000;        % points f_split..Nyquist (~0.5 Hz spacing)
n_etfe  = 2000;        % etfe frequency points (0..Nyquist, ~0.25 Hz)
w_lin = linspace(2*pi*f_split, pi/Ts, n_lin + 1);
w   = [logspace(log10(2*pi*f_lo), log10(2*pi*f_split), n_log), w_lin(2:end)];

%% Controller and identifications
K   = getfield(load(fullfile('data', 'controller', [contfn, '.mat']), 'K_final'), 'K_final');
idd = load(idfile, 'GFRF_cell', 'realmass');
wf  = 2*pi*linspace(10, 500, 5000).';                       % fine grid for the plots
G_old = interp(idd.GFRF_cell{iOP}, wf(wf <= max(idd.GFRF_cell{iOP}.Frequency)));

[y, ~, u] = utils.loadexp(newidfile);                       % new id: x_command -> Payload_acc1
seg = find(sqrt(movmean(u.^2, 500)) > 0.05*max(sqrt(movmean(u.^2, 500))));
d_new = detrend(iddata(y(seg(1):seg(end)), u(seg(1):seg(end)), K.Ts));   % as in a_MultimodelFRFID
G_new = etfe(d_new, 1000, n_etfe);                         % lag 1000, as in a_MultimodelFRFID
G_new = interp(frd(G_new.ResponseData, G_new.Frequency, Ts), ...
    w(w >= 2*pi*f_newid(1) & w <= 2*pi*f_newid(2)));        % grid w, reliable band only

%% Plant: old id vs new id (the new id contains the actuator-side factor ~2)
bo = bodeoptions;  bo.FreqUnits = 'Hz';  bo.XLim = [10 500];  bo.Grid = 'on';
bo.Title.String = sprintf('Plant at 0 kg: old id vs new id (new id on %g-%g Hz)', f_newid);
figure
bodemag(G_old, 'k', G_new, 'm', 2*G_old, 'k:', bo)
legend('old id G_{old} (Out1 command \rightarrow Wilcoxon1)', 'new id G_{new} (x\_command \rightarrow Payload\_acc1)', ...
    '2 G_{old}', 'Location', 'southeast')
exportgraphics(gcf, fullfile('plots', 'd_gaincheck_W2_plant_oldnew.pdf'), 'ContentType', 'vector')

%% Gain of each run
S = cell(size(runs, 1), 1);   g_cmd = zeros(size(runs, 1), 1);   g_S = g_cmd;
for i = 1:size(runs, 1)
    fn = fullfile(expdir, runs{i, 1});
    g_cmd(i) = utils.implgain(fn, K).gain;
    S{i} = utils.sweepS(fn);
    ws  = S{i}.Frequency(S{i}.Frequency >= 2*pi*f_band(1) & S{i}.Frequency <= 2*pi*f_band(2));
    gs = (1./squeeze(interp(S{i}, ws).ResponseData) - 1)./squeeze(freqresp(K*interp(idd.GFRF_cell{iOP}, ws), ws));
    g_S(i) = median(abs(gs));
end
disp(table(string(runs(:, 1)), [runs{:, 2}].', g_cmd, g_S, [runs{:, 2}].'./g_S, ...
    'VariableNames', {'run', 'build_gain', 'g_cmd', 'g_S', 'build_gain_for_design'}))
fprintf('80 Hz sine (build gain 1): attenuation %.1f dB, predicted %.1f dB (old id, K/%.3g)\n', ...
    utils.sineatt(fullfile(expdir, sinefn)), ...
    mag2db(abs(1/(1 + squeeze(freqresp(K*interp(idd.GFRF_cell{iOP}, 2*pi*80)/g_ref, 2*pi*80))))), g_ref);

%% Sensitivity: measured vs old id and new id (semilogx: line widths and ylim are kept)
th = {feedback(1, K*G_old/g_ref),  'r--', sprintf('old id, 1/(1+K G/%g)  (build gain 1)', g_ref)
      feedback(1, K*G_old),        'r-',  sprintf('old id, 1/(1+K G)  (build gain %g)', g_ref)
      feedback(1, g_block*K*G_new),       'b--', sprintf('new id, 1/(1+%g K G)  (build gain 1)', g_block)
      feedback(1, g_block*g_ref*K*G_new), 'b-',  sprintf('new id, 1/(1+%g K G)  (build gain %g)', g_block*g_ref, g_ref)};
fHz = @(sys) sys.Frequency/2/pi;   dB = @(sys) mag2db(abs(squeeze(sys.ResponseData)));
figure; hold on
semilogx(fHz(S{2}), dB(S{2}), 'Color', [0.85 0.33 0.10 0.6], 'DisplayName', sprintf('measured, build gain %g', runs{2, 2}))
semilogx(fHz(S{1}), dB(S{1}), 'Color', [0 0.45 0.74 0.6],    'DisplayName', sprintf('measured, build gain %g', runs{1, 2}))
for i = 1:size(th, 1)
    semilogx(fHz(th{i, 1}), dB(th{i, 1}), th{i, 2}, 'LineWidth', 1.5, 'DisplayName', th{i, 3})
end
set(gca, 'XScale', 'log'); grid on; xlim([10 500]); ylim([-60 30])
xlabel('Frequency (Hz)'); ylabel('|S| (dB)'); title('W2 controller, 0 kg: sensitivity')
legend('Location', 'southwest')
exportgraphics(gcf, fullfile('plots', 'd_gaincheck_W2_bode.pdf'), 'ContentType', 'vector')
xlim([40 120])
exportgraphics(gcf, fullfile('plots', 'd_gaincheck_W2_bode_zoomed.pdf'), 'ContentType', 'vector')

%% Nyquist (utils.nyquistcheck on the grid w, refined at K's narrow peaks)
% old id (gain 1, every operating point); new id with and without the gain of g_ref
nyq_old = utils.nyquistcheck(K, [idd.GFRF_cell{:}], 'realmass', idd.realmass, 'w', w, ...
    'name', 'L = K G, old id (gain 1)', 'file', fullfile('plots', 'd_gaincheck_W2_nyq_oldid'));
nyq_new = utils.nyquistcheck(K, [G_new, g_block*g_ref*G_new], 'w', w, 'f_range', f_newid, ...
    'labels', ["new id, K G (gain 1)"; sprintf("new id, %g K G (build gain %g)", g_block*g_ref, g_ref)], ...
    'name', sprintf('L = K G, new id (%g-%g Hz)', f_newid), 'file', fullfile('plots', 'd_gaincheck_W2_nyq_newid'));
