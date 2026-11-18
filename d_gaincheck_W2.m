%% d_gaincheck_W2  -  W2 controller (1110, order 18): gain on the setup and attenuation
%
%   Runs with the W2 controller after the sign fix (etaf 0.36), 0 kg:
%     sweep, build gain 1      sweep, build gain 4.8      80 Hz sine, build gain 1
%   For every run:
%     g_cmd   recorded command / lsim(K, -y)         -> scaling inside the controller block
%   For the sweeps (open-loop + closed-loop sweep):
%     g_plant |G_cl / GFRF(iOP)|, G_cl = (Y_cl - Y_ol)/U_cl  -> logged command -> logged acc.
%             against the identification (should be 1)
%     g_S     loop / designed loop from the measured S: (1/S - 1)/(K G)   (should be 1)
%   For the sine: attenuation (utils.sineatt) against 1/|1 + g K G| at the tone, with
%     g = 1 (design) and g = r*g_set, r = mean(g_S/g_set) of the sweeps.
%   Then the Nyquist check of the loop as it runs (r*g_set*K) and as designed (K).
%
clear; close all; clc;

%% Config
contfn  = 'LFT_GFRF_GdmGLTIW2opt_1110_etfe_80hz_fmin10_hyb_controller_order18_80hz_-34db_nperfoptiIQC_notconv_etai0.35_etaf0.36_multi_freq';
expdir  = fullfile('data', 'Experiments', [contfn, '_K_final']);
% {file, 'sweep' | 'sine', build gain}
runs = {
    'gftf0_ssweep_stbg1.mat',   'sweep', 1
    'gftf0_ssweep_stbg4.8.mat', 'sweep', 4.8
    'gfrf0_80hz_stbg1.mat',     'sine',  1
    };
idfile  = fullfile('data', 'GFRF_iddata_1110_etfe.mat');   % identification of the controller (1110)
iOP     = 1;                       % 0 kg
f_band  = [20 200];                % band for g_plant and g_S (Hz)
f_notch = 80;                      % excluded +/- 3 Hz (S ~ 0)
trans   = 3000;                    % sineatt transients

%% Controller and identification
cs  = load(fullfile('data', 'controller', [contfn, '.mat']), 'K_final');
K   = cs.K_final;
idd = load(idfile, 'GFRF', 'realmass');
G   = idd.GFRF(iOP);
w_id = idd.GFRF.Frequency;

%% Runs
nr  = size(runs, 1);
res = table('Size', [nr 9], 'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'run', 'type', 'g_set', 'g_cmd', 'cmd_resid_pct', 'g_plant', 'g_S', 'S80_dB', 'att_dB'});
res{:, 4:end} = NaN;
Sm = cell(nr, 1);
for i = 1:nr
    f_exp = fullfile(expdir, runs{i, 1});
    r = utils.implgain(f_exp, K);
    res.run(i) = string(runs{i, 1});   res.type(i) = runs{i, 2};   res.g_set(i) = runs{i, 3};
    res.g_cmd(i) = r.gain;   res.cmd_resid_pct(i) = 100*r.residual;

    switch runs{i, 2}
        case 'sweep'
            [S, Ufrd, Ycfrd, Yofrd] = utils.sweepS(f_exp);
            w = S.Frequency;   f = w/2/pi;
            sel = f >= f_band(1) & f <= f_band(2) & abs(f - f_notch) > 3 & w <= max(w_id);
            Ss  = squeeze(S.ResponseData);
            Gw  = squeeze(interp(G, w(sel)).ResponseData);
            KG  = squeeze(freqresp(K, w(sel))).*Gw;
            res.g_S(i) = median(abs((1./Ss(sel) - 1)./KG));
            Uc  = squeeze(Ufrd.ResponseData);  Yc = squeeze(Ycfrd.ResponseData);  Yo = squeeze(Yofrd.ResponseData);
            Gcl = movsum((Yc - Yo).*conj(Uc), 61)./movsum(abs(Uc).^2, 61);
            res.g_plant(i) = median(abs(Gcl(sel)./Gw));
            [~, k80] = min(abs(f - f_notch));
            res.S80_dB(i) = mag2db(abs(Ss(k80)));
            Sm{i} = S;
        case 'sine'
            res.att_dB(i) = utils.sineatt(f_exp, trans);
    end
end

% loop gain per unit of build gain, from the sweeps
isw = res.type == "sweep";
rr  = mean(res.g_S(isw)./res.g_set(isw));
fprintf('loop / designed loop = %.3f x build gain  ->  build gain for the designed loop: %.2f\n', rr, 1/rr);

% predicted attenuation of the sine runs
w80 = 2*pi*f_notch;
KG80 = squeeze(freqresp(K, w80))*squeeze(interp(G, w80).ResponseData);
res.att_pred_design_dB = NaN(nr, 1);   res.att_pred_asrun_dB = NaN(nr, 1);
for i = find(res.type == "sine").'
    res.att_pred_design_dB(i) = mag2db(abs(1/(1 + KG80)));
    res.att_pred_asrun_dB(i)  = mag2db(abs(1/(1 + rr*res.g_set(i)*KG80)));
end
disp(res)

%% Sensitivities: measured vs 1/(1 + g_S K G) and the design
figure('Name', 'W2 controller: measured S'); hold on
col = lines(nr);
for i = find(isw).'
    S = Sm{i};   w = S.Frequency;   f = w/2/pi;   m = f >= 10 & w <= max(w_id);
    KGw = squeeze(freqresp(K, w(m))).*squeeze(interp(G, w(m)).ResponseData);
    semilogx(f(m), mag2db(abs(squeeze(S.ResponseData(1, 1, m)))), 'Color', col(i, :), ...
        'DisplayName', sprintf('%s measured', res.run(i)));
    semilogx(f(m), mag2db(abs(1./(1 + res.g_S(i)*KGw))), '--', 'Color', col(i, :), ...
        'DisplayName', sprintf('1/(1 + %.2f K G)', res.g_S(i)));
end
semilogx(f(m), mag2db(abs(1./(1 + KGw))), 'k:', 'LineWidth', 1.5, 'DisplayName', 'design 1/(1 + K G)');
set(gca, 'XScale', 'log'); grid on; xlim([10 500]); legend('Location', 'best', 'Interpreter', 'none')
xlabel('Frequency (Hz)'); ylabel('|S| (dB)'); title(sprintf('W2 controller, GFRF(%d) of %s', iOP, idfile), 'Interpreter', 'none')
exportgraphics(gcf, fullfile('plots', 'd_gaincheck_W2_sens.pdf'), 'ContentType', 'vector')

%% Nyquist of the loop as designed and as it ran with each build gain
nyq_design = utils.nyquistcheck(K, idd.GFRF, 'realmass', idd.realmass, ...
    'name', 'W2: L = G K (design)', 'file', fullfile('plots', 'd_gaincheck_W2_nyq_design'));
for gset = unique(res.g_set).'
    utils.nyquistcheck(rr*gset*K, idd.GFRF, 'realmass', idd.realmass, ...
        'name', sprintf('W2: L = %.2f G K (build gain %.3g)', rr*gset, gset), ...
        'file', fullfile('plots', sprintf('d_gaincheck_W2_nyq_stbg%.3g', gset)));
end
writetable(res, fullfile('plots', 'd_gaincheck_W2.csv'));
