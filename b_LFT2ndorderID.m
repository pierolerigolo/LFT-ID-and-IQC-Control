%% b_LFT2ndorderID  -  Step 2: uncertain LFT model  G = G_dm(dm) * G_LTI * (1 + W2*Delta)
%
%   1. Grey-box fit of the mass-dependent part G_dm on fmin_dm..fmax_dm Hz
%      (one mass per operating point, common c and k).
%   2. Synthesis frequency grid w_lft (everything below is computed on it).
%   3. Residual LTI part of each operating point  G_LTI,i = G_i / G_dm,i
%      and nominal G_LTI = mean of the residuals.
%   4. Multiplicative uncertainty weight W2 covering |(G_LTI,i - G_LTI)/G_LTI|
%      (ucover): pointwise bound W1opt (w2_type = "frequency") or rational fit
%      W1 of order nW2 ("fit"). Optional high-pass (use_hpf: the uncertainty
%      below f_hpf is then dropped, so the certificate no longer covers it) and
%      lower bound on |W2| (w2_floor_db).
%   5. LFT:  dm in [-dm_bound, dm_bound] (ureal, centred, not normalised)
%            Delta: 1x1 LTI, ||Delta||_inf <= 1 (ultidyn)
%
%   Input : Data/GFRF_iddata_<expname>.mat   (from a_MultimodelFRFID)
%   Output: Data/LFT_GFRF_GdmGLTIW2<opt|n<nW2>>[hpf<f_hpf>][fl<dB>]_<expname>_<f_dist>hz_fmin<f_min>.mat
%           (for c_IQCcontrolSynW2)
%
clear; close all; clc

%% Config
expname = '1028_etfe';
fmin_dm = 25;            % band (Hz) used to fit G_dm
fmax_dm = 75;
f_dist  = 50;            % disturbance frequency (Hz) targeted by W1: a grid point
f_min   = 1;             % lower end (Hz) of the synthesis grid (upper end: Nyquist)

w2_type = "frequency";   % "frequency": pointwise bound W1opt, "fit": rational fit W1
nW2     = 2;             % order of the ucover fit (used by "fit"; ucover always fits)

use_hpf     = true;     % W2 <- Hpf*W2 (Butterworth high-pass)
n_hpf       = 6;         % order of the high-pass filter
f_hpf       = 75;       % cut-off frequency (Hz)
w2_floor_db = -40;      % lower bound on |W2| in dB (-Inf = none, phase kept)

if ~exist('plots', 'dir'), mkdir('plots'); end
load(fullfile('Data', sprintf('GFRF_iddata_%s.mat', expname)), 'GFRF', 'GFRF_cell', 'realmass');
nOP = numel(GFRF_cell);
Ts  = GFRF.Ts;

% Keep only frequencies where every FRF is finite (interp in step 1 returns
% NaN outside the original range, e.g. at Nyquist); NaN breaks ucover's fit.
ok = squeeze(all(isfinite(GFRF.ResponseData), [1 2]));
if ~all(ok)
    fprintf('Dropping %d frequencies with non-finite FRF data: %s Hz\n', ...
        nnz(~ok), mat2str(GFRF.Frequency(~ok).'/2/pi, 5));
    GFRF = frd(GFRF.ResponseData(:, :, ok), GFRF.Frequency(ok), Ts);
end

figure; bode(GFRF_cell{:})
exportgraphics(gcf, fullfile('plots', 'GFRFbode.jpg'), 'Resolution', 300)

%% 1. Mass-dependent part G_dm
[Gdm, Gdmnom, sys_est_Gdm, dm_bound, m0, m] = utils.Gdm_greyestid(GFRF, fmin_dm, fmax_dm);
exportgraphics(gcf, fullfile('plots', 'gdmcompare.jpg'), 'Resolution', 300)

T_Gdm = param_table(sys_est_Gdm, [compose("m_%d", 1:nOP), "c", "k"]);
disp(T_Gdm)
writetable(T_Gdm, fullfile('plots', 'params_Gdm.xlsx'))

%% 2. Synthesis frequency grid
% Log grid f_min .. Nyquist merged with a denser log grid on 25-100 Hz
% (300 points, the practical maximum for YALMIP); the grid point closest to
% f_dist is moved onto it (no extra point).
w_glob  = logspace(log10(2*pi*f_min), log10(pi/Ts), 150);
w_dense = logspace(log10(2*pi*25), log10(2*pi*100), 150);
w_lft   = unique([w_glob(:); w_dense(:)]);
[~, k]  = min(abs(w_lft - 2*pi*f_dist));
w_lft(k) = 2*pi*f_dist;
w_lft   = unique(w_lft);
inrange = w_lft >= min(GFRF.Frequency) & w_lft <= max(GFRF.Frequency);   % no extrapolation
if ~all(inrange)
    warning('b:grid', ['%d grid points outside the measured FRF range (%.3g-%.3g Hz) dropped: ' ...
        'extend the frequency vector w in a_MultimodelFRFID.'], nnz(~inrange), ...
        min(GFRF.Frequency)/2/pi, max(GFRF.Frequency)/2/pi);
end
w_lft   = w_lft(inrange);

%% 3. Residual LTI part and its nominal (mean), on w_lft
GFRF_w = interp(GFRF, w_lft);
clear GLTI
for i = 1:nOP
    GLTI(1, i) = GFRF_w(i)/Gdm(i); %#ok<SAGROW>
end
GLTImean  = frd(mean(GLTI.ResponseData, 2), w_lft, Ts);
GLTIarray = stack(1, GLTI(1, 1), GLTI(1, 2), GLTI(1, 3), GLTI(1, 4));

figure; bode(GLTIarray, GLTImean)
exportgraphics(gcf, fullfile('plots', 'GLTIbode.jpg'), 'Resolution', 300)

%% 4. Multiplicative uncertainty weight W2
if any(~isfinite(GLTIarray.ResponseData(:))) || any(~isfinite(GLTImean.ResponseData(:)))
    error('G_LTI has non-finite values on w_lft: ucover cannot fit W2.');
end
relerr    = GLTImean\(GLTImean - GLTIarray);
try
    [~, info] = ucover(GLTIarray, GLTImean, nW2);
catch ME
    % ucover also fits a rational weight (fitmagfrd), which can fail (e.g. with
    % the Nyquist point: "Index in position 2 exceeds array bounds").
    if w2_type == "fit", rethrow(ME); end
    warning('b:ucover', 'ucover failed (%s): W1opt computed directly.', ME.message);
    info.W1opt = frd(max(abs(relerr.ResponseData), [], 4), w_lft, Ts);   % max over the models
end

switch w2_type
    case "frequency"
        W2    = info.W1opt;          % pointwise max_i |relative error_i|
        W2frd = W2;
        W2tag = 'opt';
    case "fit"
        W2    = info.W1;             % rational fit of order nW2
        W2frd = frd(W2, w_lft);
        W2tag = sprintf('n%d', nW2);
end
W2frd.Ts = Ts;
W2_plain = W2frd;

if use_hpf
    [b_hpf, a_hpf] = butter(n_hpf, f_hpf/(1/Ts/2), 'high');
    Hpf   = tf(b_hpf, a_hpf, Ts);
    W2frd = frd(Hpf, w_lft)*W2frd;
    W2tag = sprintf('%shpf%d', W2tag, f_hpf);
else
    Hpf = tf(1, 1, Ts);
end

if isfinite(w2_floor_db)
    R      = W2frd.ResponseData;
    low    = abs(R) < db2mag(w2_floor_db);
    R(low) = db2mag(w2_floor_db)*exp(1j*angle(R(low)));
    W2frd.ResponseData = R;
    W2tag = sprintf('%sfl%.0f', W2tag, -w2_floor_db);
end

bo = bodeoptions; bo.FreqUnits = 'Hz';
figure; bodemag(relerr, 'b--', W2_plain, 'r', W2frd, 'm:', bo);
legend('relative error', 'W2', 'W2 used (hpf/floor)', 'Location', 'best')
xlim([10 500])
exportgraphics(gcf, fullfile('plots', 'relerrW2bode.jpg'), 'Resolution', 300)

%% 5. LFT
% G_LTI*(1 + W2*Delta) = Fu([0, W2; G_LTI, G_LTI], Delta)
GLTI_H  = [0,        W2frd;
           GLTImean, GLTImean];
GLTILFT = lft(ultidyn('d', 1), GLTI_H);

% G_dm(m0 - dm) = Gdmnom/(1 + Gdmnom*dm) = Fu([-Gdmnom, Gdmnom; -Gdmnom, Gdmnom], dm)
Gdmnom_w = frd(Gdmnom, w_lft);
Gdm_H    = [-Gdmnom_w, Gdmnom_w;
            -Gdmnom_w, Gdmnom_w];
GdmLFT   = lft(ureal('dm', 0, 'PlusMinus', dm_bound), Gdm_H);

GLFT_delta    = GdmLFT*GLTILFT;
[GLFT, Delta] = lftdata(GLFT_delta)

save(fullfile('Data', sprintf('LFT_GFRF_GdmGLTIW2%s_%s_%dhz_fmin%g.mat', W2tag, expname, f_dist, f_min)), ...
    "GLFT_delta", "GLFT", "Delta", "dm_bound", "m0", "m", "GFRF", "GFRF_cell", ...
    "realmass", "W2", "W2frd", "Hpf", "GLTI", "GLTImean", "relerr")

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
