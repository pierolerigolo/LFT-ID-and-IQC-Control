%% b_LFTdmhatEstimate  -  Step 2b: LFT of the plant with the mass-estimate compensation
%
%   On the setup the controller output is multiplied by (1 + dmhat*G_nom) (Simulink:
%   "Nominal 2nd Order Model" -> dmhat -> sum) before the plant. With the LFT of
%   b_LFT2ndorderID,  G = G_LTI (1 + W2 Delta) * G_nom/(1 + G_nom dm),  the plant seen by K is
%
%       G_eff = G_LTI (1 + W2 Delta) * G_nom (1 + G_nom dmhat)/(1 + G_nom dm)          (exact)
%             = G_LTI (1 + W2 Delta) * G_nom / (1 + eps G'),   G' = G_nom/(1 + G_nom dmhat)
%
%   Estimation error  eps = dm - dmhat  (same units and sign as dm in b: the mass offset
%   of the grey-box, dm = 1/G_dm(m) - 1/G_nom, identified-mass units; dmhat must be given
%   to the compensator with that same convention).
%   eps is the new uncertain parameter. It is declared over the same range as dm in b,
%   eps in [-dm_bound, dm_bound], so that after the normalisation of lftdata
%       delta_eps = eps/dm_bound   (the mass channel of GLFT)
%   and eta in the c scripts IS the estimation error as a fraction of the half-range of
%   identified masses:
%       eta = |dm - dmhat|/dm_bound      (eta = 0.2: estimate within 20 % of the half-range;
%                                          eta = 1: error as large as the whole mass spread)
%   Set eta0_list in the c script to the accuracy you expect from the estimator (e.g. 0.2);
%   the certified eta is the largest estimation error (as that fraction) the controller
%   tolerates. In mass units: eps = eta*dm_bound (dm_bound saved unchanged).
%
%   dmhat_mode
%     "fixed" : dmhat = dmhat_design (a number), LFT in (Delta_W2, eps) -> same structure as
%               b, usable as is in c_IQCcontrolSynW2 / W2RP / NoW2. With dmhat_design = 0
%               the eps channel goes through G_nom: exact for an estimate near the nominal
%               mass, approximate otherwise (G' = G_dm(m0 + dmhat) varies with dmhat).
%     "range" : dmhat in [-dm_bound, dm_bound] as a second real uncertain parameter (K does
%               not know dmhat): exact for any estimate, LFT in (Delta_W2, dmhat x2, eps).
%               Needs a multiplier with separate blocks for dmhat and eps (not yet in c).
%
%   Measured FRFs for the checks of the c scripts (encirclements, Nyquist, tables): each
%   operating point is compensated with its own estimate dmhat_i = dm_i + eps_test*dm_bound
%   (eps_test = 0: perfect estimate), GFRF = GFRF_raw (1 + dmhat_i G_nom).
%
%   Input : Data/<lftname>.mat from b_LFT2ndorderID
%   Output: Data/<lftname>_est<mode>.mat (GLFT_delta, GLFT, Delta, dm_bound,
%           GFRF, ... ) -> idname of the c scripts.
%
clear; close all; clc

%% Config
lftname      = 'LFT_GFRF_GdmGLTIW2opt_1110gc_etfe_80hz_fmin10_hyb';
dmhat_mode   = "fixed";    % "fixed" | "range"
dmhat_design = 0;          % dmhat of the "fixed" LFT (identified-mass units)
eps_test     = 0;          % estimate error used to compensate the measured FRFs (fraction of dm_bound)

%% Load b
L = load(fullfile('Data', [lftname, '.mat']));
w_lft = L.GLFT.Frequency;   Ts = L.GLFT.Ts;   dm_bound = L.dm_bound;
GLTIwW2 = lft(ultidyn('d', 1), [frd(zeros(1, 1, numel(w_lft)), w_lft, Ts), L.W2frd;
                                 L.GLTImean,                               L.GLTImean]);
if isfield(L, 'Gdmnom')
    Gn = frd(L.Gdmnom, w_lft);   Gn.Ts = Ts;
else
    % older b files: G_nom = nominal plant / G_LTI (exact on w_lft)
    Gn = L.GLFT_delta.NominalValue / L.GLTImean;
end

%% Mass offset of each operating point in the convention of b: dm_i = 1/G_dm(m_i) - 1/G_nom
% (G_dm(m_i) = GFRF_i / G_LTI_i; median over frequency, sign-safe whatever the phase convention)
nOP = size(L.GFRF, 2);
dm_i = zeros(1, nOP);
for i = 1:nOP
    Gdm_i  = squeeze(interp(L.GFRF(i), w_lft).ResponseData)./squeeze(L.GLTI(i).ResponseData);
    dm_i(i) = median(real(1./Gdm_i - 1./squeeze(Gn.ResponseData)));
end
fprintf('dm_bound = %.4g, mass offsets of the operating points dm_i = %s\n', dm_bound, mat2str(dm_i, 3));

%% Compensated LFT in the estimation error eps (declared over +/- dm_bound: eta = |eps|/dm_bound)
epsu    = ureal('eps', 0, 'PlusMinus', dm_bound);
switch dmhat_mode
    case "fixed"
        dmh = dmhat_design;
        tag = sprintf('estfix%g', dmhat_design);
    case "range"
        dmh = ureal('dmhat', 0, 'PlusMinus', dm_bound);
        tag = 'estrange';
end
Gcomp      = Gn*(1 + Gn*dmh)/(1 + Gn*(dmh + epsu));       % exact compensated mass part
GLFT_delta = GLTIwW2*Gcomp;
[GLFT, Delta, BlkStruct] = lftdata(GLFT_delta);
disp(GLFT_delta)
fprintf('lftdata blocks: %s\n', strjoin(string({BlkStruct.Name}) + "(" + string([BlkStruct.Occurrences]) + ")", ', '));
fprintf('mass channel of GLFT = eps/dm_bound: eta in the c scripts = |dm - dmhat|/dm_bound (dm_bound = %.4g)\n', dm_bound);

%% Compensated measured FRFs (each with its own estimate)
GFRF_raw = L.GFRF;
GFRF     = GFRF_raw;
wF       = GFRF_raw.Frequency;
Gn_F     = squeeze(interp(Gn, wF).ResponseData);     % G_nom on the FRF grid (linear interp.)
for i = 1:nOP
    dmhat_i = dm_i(i) + eps_test*dm_bound;
    GFRF.ResponseData(1, i, :) = squeeze(GFRF_raw.ResponseData(1, i, :)).*(1 + dmhat_i*Gn_F);
end

figure; hold on
for i = 1:nOP
    utils.magplot(GFRF_raw(i), 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    utils.magplot(GFRF(i), 'DisplayName', sprintf('compensated OP%d', i));
end
utils.magplot(GLFT_delta.NominalValue, 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'nominal');
title(sprintf('Measured FRFs (grey) and with the estimate compensation (eps\\_test = %g)', eps_test))

%% Save (names expected by the c scripts)
realmass = L.realmass;   m0 = L.m0;   m = L.m;   W2frd = L.W2frd;   GLTI = L.GLTI;   GLTImean = L.GLTImean;
outname = sprintf('%s_%s', lftname, tag);
save(fullfile('Data', [outname, '.mat']), 'GLFT_delta', 'GLFT', 'Delta', 'dm_bound', ...
    'dmhat_mode', 'dmhat_design', 'dm_i', 'eps_test', 'GFRF', 'GFRF_raw', 'realmass', ...
    'm0', 'm', 'W2frd', 'GLTI', 'GLTImean', 'Gn')
fprintf('saved Data/%s.mat  (c scripts: idname = ''%s'')\n', outname, outname);
if dmhat_mode == "range"
    warning('b:range', ['"range" LFT has the blocks d, dmhat (x2) and eps: the c scripts treat all ' ...
        'real channels as ONE block with one eta, so they need a block-diagonal multiplier first.']);
end
