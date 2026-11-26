%% a_GainCorrectID  -  Step 1b: gain-corrected 1110 identification (until the full re-id)
%
%   The 1110 identification (Out1 command -> Wilcoxon1, identification Simulink model)
%   differs from the plant the controller sees (x_command -> Payload_acc1, controller
%   Simulink model, times the 0.1 inside the controller block) by an (almost) constant
%   real gain. It is estimated from the two 0 kg sweeps, computed exactly as in
%   a_MultimodelFRFID:
%       G_old = 261110 sweep (identification model)
%       G_new = 0.1 * 261118 sweep through x_command (controller model, 0.1 folded in)
%       c     = median |G_new/G_old| over f_band where both coherences > coh_min
%   Every FRF of GFRF_iddata_1110_etfe is multiplied by c (and the inputs of d by 1/c)
%   and saved as GFRF_iddata_1110gc_etfe.mat: run b_LFT2ndorderID with
%   expname = '1110gc_etfe'. A controller designed on it runs with build gain 1.
%
clear; close all; clc

%% Config
idfile_in  = fullfile('Data', 'GFRF_iddata_1110_etfe.mat');
idfile_out = fullfile('Data', 'GFRF_iddata_1110gc_etfe.mat');
oldidfile  = fullfile('Data', 'MIVIDA 261110', '261110_gfrf0kg_xsswp.mat');                 % Out1 -> Wilcoxon1
newidfile  = fullfile('Data', 'MIVIDA 261118 sswp new simulink', 'gfrf0_ssweep_id.mat');   % x_command -> Payload_acc1
g_block = 0.1;         % factor inside the controller block (folded into G_new)
f_band  = [20 200];    % band for the estimate (Hz)
coh_min = 0.9;         % minimum coherence of both sweeps

% FRF exactly as in a_MultimodelFRFID
Ts      = 0.001;
f_lo    = 10;  f_split = 20;  n_log = 50;  n_lin = 1000;  n_etfe = 2000;
w_lin = linspace(2*pi*f_split, pi/Ts, n_lin + 1);
w     = [logspace(log10(2*pi*f_lo), log10(2*pi*f_split), n_log), w_lin(2:end)];
frf   = @(u, y) interp(etfe(detrend(iddata(y(:), u(:), Ts)), 1000, n_etfe), w);
coh   = @(u, y) interp1(mscohere(detrend(u(:)), detrend(y(:)), hann(4096), 2048, 4096, 1/Ts), ...
                        (0:2048).'/4096/Ts, w/2/pi);

%% Old and new 0 kg identification
[u, y] = utils.openIOtemp(oldidfile);      G_old = frf(u, y);           C_old = coh(u, y);
[y, ~, u] = utils.loadexp(newidfile);      G_new = g_block*frf(u, y);   C_new = coh(u, y);

%% Gain estimate
r   = squeeze(G_new.ResponseData)./squeeze(G_old.ResponseData);
sel = w(:)/2/pi >= f_band(1) & w(:)/2/pi <= f_band(2) & C_old(:) > coh_min & C_new(:) > coh_min;
c   = median(abs(r(sel)));
ph  = rad2deg(median(angle(r(sel))));
Go  = squeeze(G_old.ResponseData);  Gn = squeeze(G_new.ResponseData);
c_ls = (Go(sel)'*Gn(sel))/(Go(sel)'*Go(sel));                    % complex least squares
fprintf(['gain G_new/G_old (%g-%g Hz, %d points): c = %.4f (= 1/%.3f), IQR %.4f-%.4f, ' ...
    'median phase %+.1f deg; least squares %.4f at %+.1f deg\n'], f_band, nnz(sel), c, 1/c, ...
    prctile(abs(r(sel)), 25), prctile(abs(r(sel)), 75), ph, abs(c_ls), rad2deg(angle(c_ls)));
if abs(ph) > 10
    warning('a:phase', 'The ratio has a phase of %.1f deg: not a pure gain.', ph);
end

figure
tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile
bo = bodeoptions;  bo.FreqUnits = 'Hz';  bo.XLim = [10 500];  bo.Grid = 'on';  bo.Title.String = '';
bodemag(G_old, 'k', c*G_old, 'k--', G_new, 'm', bo)
legend('G_{old} (1110)', sprintf('c G_{old}, c = %.3f', c), 'G_{new} (0.1 x\\_command \\rightarrow Payload\\_acc1)', ...
    'Location', 'southeast')
title('0 kg: old id, gain-corrected old id, new id')
nexttile
semilogx(w/2/pi, abs(r), 'Color', [0.6 0.6 0.6]); hold on
semilogx(w(sel)/2/pi, abs(r(sel)), 'b.', 'MarkerSize', 4)
yline(c, 'r', sprintf('c = %.3f', c), 'LineWidth', 1.5)
grid on; xlim([10 500]); ylim([0 2*c]); xlabel('Frequency (Hz)'); ylabel('|G_{new}/G_{old}|')
legend('all', sprintf('used (%g-%g Hz, coherence > %.2f)', f_band, coh_min), 'Location', 'best')
exportgraphics(gcf, fullfile('plots', 'gaincorrect_1110.pdf'), 'ContentType', 'vector')

%% Gain-corrected 1110 identification
load(idfile_in, 'GFRF_cell', 'GFRF', 'd', 'realmass')
GFRF = c*GFRF;
for i = 1:numel(GFRF_cell)
    GFRF_cell{i}.ResponseData = c*GFRF_cell{i}.ResponseData;
    d{i}.InputData = d{i}.InputData/c;                        % y/u scaled by c as well
end
gain_corr = c;   gain_band = f_band;
gain_note = sprintf(['1110 id times c = %.4f: plant from the controller output (build gain 1, ' ...
    '0.1 of the controller block included) to Payload_acc1, from the 261118 sweep through x_command.'], c);
save(idfile_out, 'GFRF_cell', 'GFRF', 'd', 'realmass', 'gain_corr', 'gain_band', 'gain_note')
fprintf('saved %s  (b_LFT2ndorderID: expname = ''1110gc_etfe'')\n', idfile_out);
