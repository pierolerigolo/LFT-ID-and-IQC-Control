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
%   New identification (sweep sent through x_command with the controller Simulink model):
%     G_new, compared with the old id; it contains the actuator-side factor, so the loop
%     at build gain g is g_block*g*K*G_new (g_block = 0.1, inside the controller block).
%   Plots (magnitude): measured S vs theory with the old id; measured S vs theory with
%   the old and the new id. Nyquist: K G_old (gain 1, old id), K G_new (gain 1, new id),
%   g_block*g_ref*K G_new (new id with the gain correction).
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
g_ref   = 4.8;                     % build gain giving the designed loop with the old id
g_block = 0.1;                     % factor inside the controller block (logged acc. -> logged command)
newidfile = fullfile('data', 'MIVIDA 261118 sswp new simulink', 'gfrf0_ssweep_id.mat');  % sweep through x_command
Ts_exp  = 0.001;

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

%% New identification: plant seen from x_command (controller Simulink model)
% open-loop sweep sent through x_command -> includes the actuator-side factor (~2)
[y_id, ~, u_id] = utils.loadexp(newidfile);
m_id  = sqrt(movmean(u_id.^2, 500));   act = find(m_id > 0.05*max(m_id));
seg   = act(1):act(end);                                % sweep part
nwin  = 4096;
[H_new, f_new] = tfestimate(u_id(seg), y_id(seg), hann(nwin), nwin/2, nwin, 1/Ts_exp);
keep  = f_new >= 5 & f_new <= 1/(2*Ts_exp);
G_new = frd(reshape(H_new(keep), 1, 1, []), 2*pi*f_new(keep), Ts_exp);
wcmp  = w_id(w_id >= min(G_new.Frequency) & w_id <= max(G_new.Frequency));
r_new = squeeze(interp(G_new, wcmp).ResponseData)./squeeze(interp(G, wcmp).ResponseData);
bd    = wcmp/2/pi >= f_band(1) & wcmp/2/pi <= f_band(2);
fprintf('new id / old id GFRF(%d): |ratio| median %.3f, phase %+.1f deg (20-200 Hz)\n', ...
    iOP, median(abs(r_new(bd))), rad2deg(median(angle(r_new(bd)))));

%% Bode (magnitude): measured sensitivities against the theory (old id)
% g_ref = build gain that gives the designed loop (4.8 here): the run at build gain 1
% then has the loop K/g_ref, the run at g_ref the designed loop K.
wb  = linspace(2*pi*10, max(w_id), 4000).';
KGb = squeeze(freqresp(K, wb)).*squeeze(interp(G, wb).ResponseData);
fig = figure('Name', 'Bode: measured vs theoretical S (old id)'); hold on
cs = lines(2);   k = 0;
for i = find(isw).'
    k = k + 1;   S = Sm{i};   w = S.Frequency;   f = w/2/pi;   m = f >= 10 & w <= max(w_id);
    semilogx(f(m), mag2db(abs(squeeze(S.ResponseData(1, 1, m)))), 'Color', [cs(k, :) .6], ...
        'DisplayName', sprintf('measured, build gain %.3g', res.g_set(i)));
end
semilogx(wb/2/pi, mag2db(abs(1./(1 + KGb))),       'k',   'LineWidth', 1.5, 'DisplayName', '1/(1 + K G)  (design gain)');
semilogx(wb/2/pi, mag2db(abs(1./(1 + KGb/g_ref))), 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('1/(1 + K G/%.3g)  (build gain 1)', g_ref));
set(gca, 'XScale', 'log'); xlim([10 500]); grid on; legend('Location', 'southwest')
xlabel('Frequency (Hz)'); ylabel('|S| (dB)'); title(sprintf('W2 controller, 0 kg: measured S and theory, old id GFRF(%d)', iOP))
exportgraphics(fig, fullfile('plots', 'd_gaincheck_W2_bode.pdf'), 'ContentType', 'vector')

%% Bode (magnitude): old id vs new id
% old id: loop = (build gain/g_ref) K G_old;  new id: loop = g_block * build gain * K G_new
wn   = linspace(max(2*pi*10, min(G_new.Frequency)), min(max(w_id), max(G_new.Frequency)), 4000).';
KGn  = squeeze(freqresp(K, wn)).*squeeze(interp(G_new, wn).ResponseData);
KGo  = squeeze(freqresp(K, wn)).*squeeze(interp(G, wn).ResponseData);
fig = figure('Name', 'Bode: old id vs new id'); hold on
k = 0;
for i = find(isw).'
    k = k + 1;   S = Sm{i};   w = S.Frequency;   f = w/2/pi;   m = f >= 10 & w <= max(w_id);
    semilogx(f(m), mag2db(abs(squeeze(S.ResponseData(1, 1, m)))), 'Color', [cs(k, :) .5], ...
        'DisplayName', sprintf('measured, build gain %.3g', res.g_set(i)));
end
for gset = [1 g_ref]
    if gset == 1, ls = '--'; else, ls = '-'; end
    semilogx(wn/2/pi, mag2db(abs(1./(1 + gset/g_ref*KGo))), ['k', ls], 'LineWidth', 1.3, ...
        'DisplayName', sprintf('old id: 1/(1 + %.3g K G_{old})  (build gain %.3g)', gset/g_ref, gset));
    semilogx(wn/2/pi, mag2db(abs(1./(1 + g_block*gset*KGn))), ['r', ls], 'LineWidth', 1.3, ...
        'DisplayName', sprintf('new id: 1/(1 + %.3g K G_{new})  (build gain %.3g)', g_block*gset, gset));
end
set(gca, 'XScale', 'log'); xlim([10 500]); grid on; legend('Location', 'southwest')
xlabel('Frequency (Hz)'); ylabel('|S| (dB)'); title('W2 controller, 0 kg: measured S, theory with the old and the new id')
exportgraphics(fig, fullfile('plots', 'd_gaincheck_W2_bode_oldnew.pdf'), 'ContentType', 'vector')

%% Nyquist: design with the old id, new id, new id with the gain correction
nyq_old = utils.nyquistcheck(K, idd.GFRF, 'realmass', idd.realmass, ...
    'name', 'W2: L = K G_{old} (gain 1, old id)', 'file', fullfile('plots', 'd_gaincheck_W2_nyq_oldid'));
nyq_new = utils.nyquistcheck(K, G_new, 'realmass', 0, ...
    'name', 'W2: L = K G_{new} (gain 1, new id)', 'file', fullfile('plots', 'd_gaincheck_W2_nyq_newid'));
nyq_cor = utils.nyquistcheck(g_block*g_ref*K, G_new, 'realmass', 0, ...
    'name', sprintf('W2: L = %.3g K G_{new} (new id, build gain %.3g)', g_block*g_ref, g_ref), ...
    'file', fullfile('plots', 'd_gaincheck_W2_nyq_newid_corr'));
writetable(res, fullfile('plots', 'd_gaincheck_W2.csv'));
