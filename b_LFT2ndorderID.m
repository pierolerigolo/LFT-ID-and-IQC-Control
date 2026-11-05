%% b_LFT2ndorderID  -  Step 2: uncertain LFT model  G = G_dm(dm) * G_LTI * (1 + W2*Delta)
%
%   1. Grey-box fit of the mass-dependent part G_dm on 25-75 Hz
%      (one mass per operating point, common c and k).
%   2. Residual LTI part of each operating point  G_LTI,i = G_i / G_dm,i
%      and nominal G_LTI = mean of the residuals.
%   3. Multiplicative uncertainty weight W2 covering
%      |(G_LTI,i - G_LTI)/G_LTI| (ucover), optionally times a high-pass
%      filter (use_hpf): the residual uncertainty below f_hpf is then dropped,
%      so the robustness certificate no longer covers it there.
%   4. LFT:  dm in [-dm_bound, dm_bound] (ureal, centred, not normalised)
%            Delta: 1x1 LTI, ||Delta||_inf <= 1 (ultidyn)
%
%   Input : Data/GFRF_iddata_<expname>.mat              (from a_MultimodelFRFID)
%   Output: Data/LFT_GFRF_GdmGLTIW2n<nW2>[hpf<f_hpf>]_<expname>.mat (for c_IQCcontrolSynW2)
%
clear; close all; clc

%% Config
expname = '1028_etfe';
fmin_dm = 25;          % band (Hz) used to fit G_dm
fmax_dm = 75;
nW2     = 25;          % order of the W2 fit (ucover)
f_dist  = 80;        % disturbance frequency (Hz) targeted by W1: kept ON the grid
use_hpf = true;      % W2 <- Hpf*W2 (Butterworth high-pass)
n_hpf   = 6;         % order of the high-pass filter
f_hpf   = 100;       % cut-off frequency (Hz)
w2_floor_db = -40;   % lower bound on |W2| in dB (-Inf = none). Keeps the uncertainty
                     % channel from vanishing below f_hpf (6th order: -120 dB at 10 Hz),
                     % which makes the IQC LMIs badly scaled. Raising |W2| only adds
                     % uncertainty (conservative); the phase of W2 is kept.

if ~exist('plots', 'dir'), mkdir('plots'); end
load(fullfile('Data', sprintf('GFRF_iddata_%s.mat', expname)), 'GFRF', 'GFRF_cell', 'realmass');
nOP = numel(GFRF_cell);

figure; bode(GFRF_cell{:})
exportgraphics(gcf, fullfile('plots', 'GFRFbode.jpg'), 'Resolution', 300)

%% 1. Mass-dependent part G_dm
[Gdm, Gdmnom, sys_est_Gdm, dm_bound, m0, m] = utils.Gdm_greyestid(GFRF, fmin_dm, fmax_dm);
exportgraphics(gcf, fullfile('plots', 'gdmcompare.jpg'), 'Resolution', 300)

T_Gdm = param_table(sys_est_Gdm, [compose("m_%d", 1:nOP), "c", "k"]);
disp(T_Gdm)
writetable(T_Gdm, fullfile('plots', 'params_Gdm.xlsx'))

%% 2. Residual LTI part and its nominal (mean)
clear GLTI
figure; hold on
for i = 1:nOP
    GLTI(1, i) = GFRF(i)/Gdm(i); %#ok<SAGROW>
    bode(GLTI(1, i))
end
GLTImean = frd(mean(GLTI.ResponseData, 2), GLTI.Frequency, GLTI.Ts);
bode(GLTImean)
exportgraphics(gcf, fullfile('plots', 'GLTIbode.jpg'), 'Resolution', 300)

%% 3. Multiplicative uncertainty weight W2
% ucover is run on a frequency axis scaled by a (better conditioning) and the
% continuous-time fit is scaled back: W2(s) = Ws(s/a).
a  = 1e3;
ws = GLTImean.Frequency/a;
Pn_s     = frd(GLTImean.ResponseData, ws);
Parray_s = cell(1, nOP);
Parray   = cell(1, nOP);
for i = 1:nOP
    Parray_s{i} = frd(GLTI(1, i).ResponseData, ws);
    Parray{i}   = GLTI(1, i);
end
[~, info_s] = ucover(stack(1, Parray_s{:}), Pn_s, nW2);
Ws = ss(info_s.W1);
W2 = ss(a*Ws.A, a*Ws.B, Ws.C, Ws.D);

relerr = GLTImean\(GLTImean - stack(1, Parray{:}));
bo = bodeoptions; bo.FreqUnits = 'Hz';
figure; bodemag(relerr, 'b--', W2, 'r', bo);
exportgraphics(gcf, fullfile('plots', 'relerrW2bode.jpg'), 'Resolution', 300)

% Optional high-pass filter on W2 (discrete Butterworth)
if use_hpf
    [b_hpf, a_hpf] = butter(n_hpf, f_hpf/(1/GFRF.Ts/2), 'high');
    Hpf = tf(b_hpf, a_hpf, GFRF.Ts);
    W2plot    = frd(W2, GLTImean.Frequency);
    W2plot.Ts = GFRF.Ts;
    figure; bodemag(relerr, 'b--', W2plot, 'r', frd(Hpf, GLTImean.Frequency), 'k:', ...
        frd(Hpf, GLTImean.Frequency)*W2plot, 'm', bo);
    legend('relative error', 'W2', 'Hpf', 'Hpf*W2', 'Location', 'best')
    xlim([10 500])
    exportgraphics(gcf, fullfile('plots', 'w2relerrhpf.jpg'), 'Resolution', 300)
    W2tag = sprintf('n%dhpf%d', nW2, f_hpf);
else
    Hpf   = tf(1, 1, GFRF.Ts);
    W2tag = sprintf('n%d', nW2);
end

%% 4. LFT on the synthesis frequency grid
omega_global     = logspace(log10(2*pi*10), log10(pi/GFRF.Ts), 150);   % 10 Hz .. Nyquist
omega_peak_dense = logspace(log10(25*2*pi), log10(100*2*pi), 150);     % dense 25-100 Hz
w_lft = unique([omega_global(:); omega_peak_dense(:); 2*pi*f_dist]);

% G_LTI*(1 + W2*Delta) = Fu([0, W2; G_LTI, G_LTI], Delta)
GLTI_selected = interp(GLTImean, w_lft);
W2frd    = frd(W2, w_lft);      % continuous-time W2 sampled in rad/s
W2frd.Ts = GFRF.Ts;             % same Ts as GLTI_selected
W2frd    = frd(Hpf, w_lft)*W2frd;   % high-pass (identity if ~use_hpf)
if isfinite(w2_floor_db)
    R     = W2frd.ResponseData;
    magR  = abs(R);
    floorv = db2mag(w2_floor_db);
    low   = magR < floorv;
    R(low) = floorv*exp(1j*angle(R(low)));   % same phase, magnitude = floor
    W2frd.ResponseData = R;
    fprintf('W2 floor %.0f dB applied at %d of %d frequencies\n', w2_floor_db, nnz(low), numel(R));
    W2tag = sprintf('%sfl%.0f', W2tag, -w2_floor_db);
end
GLTI_H   = [0,             W2frd;
            GLTI_selected, GLTI_selected];
GLTILFT  = lft(ultidyn('d', 1), GLTI_H);

% G_dm(m0 - dm) = Gdmnom/(1 + Gdmnom*dm) = Fu([-Gdmnom, Gdmnom; -Gdmnom, Gdmnom], dm)
Gdmnom = frd(Gdmnom, w_lft);
Gdm_H  = [-Gdmnom, Gdmnom;
          -Gdmnom, Gdmnom];
GdmLFT = lft(ureal('dm', 0, 'PlusMinus', dm_bound), Gdm_H);

GLFT_delta    = GdmLFT*GLTILFT;
[GLFT, Delta] = lftdata(GLFT_delta)

save(fullfile('Data', sprintf('LFT_GFRF_GdmGLTIW2%s_%s.mat', W2tag, expname)), ...
    "GLFT_delta", "GLFT", "Delta", "dm_bound", "m0", "m", "GFRF", "GFRF_cell", ...
    "realmass", "W2", "Hpf", "W2frd", "GLTImean")

%% Local functions
function T = param_table(sys_est, names)
%PARAM_TABLE Parameter estimates and standard deviations of a greyest model.
%   names: one entry per estimated parameter (ParVector order)
    P   = sys_est.Report.Parameters;
    val = P.ParVector(:);
    C   = P.FreeParCovariance;
    if isempty(C) || all(isnan(C(:)))
        sd = nan(size(val));
        warning('param_table:noCov', 'No covariance available for this model.')
    else
        sd = sqrt(diag(C));
    end
    names = string(names(:));
    assert(numel(names) == numel(val), 'Expected %d names, got %d.', numel(val), numel(names));
    T = table(names, val, sd, 100*sd./abs(val), ...
        'VariableNames', {'Parameter', 'Value', 'StdDev', 'RelStd_pct'});
end
