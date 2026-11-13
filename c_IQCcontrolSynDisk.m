%% c_IQCcontrolSynDisk  -  Step 3 (alternative): robust controller with a disk uncertainty set
%
%   Additive disk uncertainty on the measured FRFs (no mass block, no W2),
%   Hindi, Seong, Boyd (2002), "Computing optimal uncertainty models from
%   frequency domain data", CDC  ->  MATLAB ucover (Robust Control Toolbox):
%       G(jw) in { G_hat(jw) + r(jw) delta : |delta| <= 1, delta complex }
%   Nominal model G_hat (nominal_type):
%     "chebyshev"  centre of the smallest disk holding all the operating points at
%                  every frequency (optimal nominal, smallest radius)
%     "mean"       average of the measured FRFs
%     "elliptical" nominal of find_elliptical_uncertainty (same as c_IQCcontrolSynElli)
%   Radius r(jw): ucover(..., 'Additive') around G_hat; disk_radius = "pointwise"
%   uses its pointwise bound W1opt = max_i |G_i - G_hat|, "fitted" the rational
%   weight of order nWdisk (always >= W1opt).
%   The disk is the elliptical set with A(jw) = I/r(jw): same IQC, multiplier
%   PI = diag(1, -1/r^2). Same frequency grid, measured FRFs, W1 and stages as
%   c_IQCcontrolSynElli:
%     1  robust stability only (slack minimised)
%     2  min nominal H-inf  s.t. robust stability
%   with the same open-loop-stable fallbacks.
%
%   Output: Data/controller/DISK_<expname>_controller_... .mat and the run diary.
%
close all; clc; clear; clear("yalmip");

%% Config
idname        = 'LFT_GFRF_GdmGLTIW2opt_1028_etfe';  % only its frequency grid and measured FRFs are used
orderK_list   = 10:2:20;     % controller orders to synthesise
Ts            = 0.001;
seed          = 22;

% Disk uncertainty
nominal_type  = "chebyshev"; % "chebyshev" | "mean" | "elliptical"
disk_radius   = "pointwise"; % "pointwise" (ucover W1opt) | "fitted" (ucover rational fit)
nWdisk        = 4;           % order of the ucover fit

% Performance weight W1 = floor_gain + peak_extra * band-pass at omega_d
omega_d       = 80*2*pi;   % disturbance frequency (rad/s)
floor_gain    = 0.5;
peak_extra    = 500;
zeta_d        = 0.001;
zeta_n        = 0.0001;

% IQC
alpha_PI        = 0;         % split of the multiplier (Gupta et al. 2024, Cor. 1)
beta_PI         = 0.5;

% Optimisation
solver          = 'mosek';
mosek_near_rel  = 1e6;       % MOSEK default 1e3
max_iter        = 30;
max_iter_stage2 = 20;        % stage 2 ends after this many iterations even if the
                             % objective still oscillates (Inf = only on convergence)
obj_err         = [5e-3, 1e-3];   % convergence tolerance per stage

% Solver paths (machine specific)
addpath(genpath("C:\Users\admin\mosek"))
addpath(genpath("C:\Program Files\Mosek\10.1"))
addpath(genpath("C:\Users\admin\Documents\MATLAB\tbxmanager"))
% addpath(genpath("/Users/pierre/mosek"))
% addpath(genpath("/Users/pierre/Documents/MATLAB/tbxmanager"))

%% Load the measured FRFs and the frequency grid
load(fullfile('Data', sprintf('%s.mat', idname)), 'GLFT', 'GFRF')
W_lft   = GLFT.Frequency;
nOP     = size(GFRF, 2);
GFRF_w  = interp(GFRF, W_lft);
diskname = ['DISK_', regexprep(idname, '^.*?_(\d{4}_.*)$', '$1')];   % e.g. DISK_1028_etfe

%% Nominal model and disk uncertainty (ucover, Hindi et al. 2002)
clear frdmodels
for i = 1:nOP
    frdmodels(1, 1, i) = GFRF_w(i); %#ok<SAGROW>   1 x 1 x nOP model array
end
Z = squeeze(frdmodels.ResponseData);           % nw x nOP
if isvector(Z), Z = Z(:); end
if size(Z, 1) ~= numel(W_lft), Z = Z.'; end
switch nominal_type
    case "chebyshev"
        c = zeros(numel(W_lft), 1);
        for k = 1:numel(W_lft)
            c(k) = min_disk(Z(k, :));
        end
        G_hat = frd(reshape(c, 1, 1, []), W_lft, GFRF_w.Ts);
    case "mean"
        G_hat = frd(reshape(mean(Z, 2), 1, 1, []), W_lft, GFRF_w.Ts);
    case "elliptical"
        model_elliptical = datadriven.find_elliptical_uncertainty(frdmodels, W_lft, 'verbosity', 1);
        G_hat = model_elliptical.NominalModel;
end
[ny, nu] = size(G_hat);

r_opt = max(abs(Z - squeeze(G_hat.ResponseData)), [], 2);   % pointwise bound
W_disk = [];
try
    [~, info_disk] = ucover(frdmodels, G_hat, nWdisk, [], 'Additive');
    W_disk = info_disk.W1;
    r_opt  = abs(squeeze(info_disk.W1opt.ResponseData));
catch ME
    % ucover's rational fit (fitmagfrd) can fail; the pointwise bound does not need it
    warning('c:ucover', 'ucover failed (%s): pointwise radius computed directly.', ME.message);
end
switch disk_radius
    case "pointwise"
        r = r_opt(:);
    case "fitted"
        if isempty(W_disk), error('disk_radius = "fitted" needs the ucover fit.'); end
        r = max(abs(squeeze(freqresp(W_disk, W_lft))), r_opt(:));
end
r = max(r, 1e-6*max(r));                        % no zero radius (PI = diag(1, -1/r^2))
bad = ~isfinite(r) | ~isfinite(squeeze(G_hat.ResponseData));
if any(bad)
    error('Disk set / nominal model not finite at %s Hz.', mat2str(W_lft(bad).'/2/pi, 4));
end
model_disk = struct('NominalModel', G_hat, 'Radius', frd(reshape(r, 1, 1, []), W_lft, Ts), ...
    'RadiusPointwise', r_opt(:), 'Weight', W_disk, 'nominal_type', nominal_type, ...
    'disk_radius', disk_radius);

% FRFs, nominal model and |G_hat| +/- r
Gmag = abs(squeeze(G_hat.ResponseData));
figure; hold on
for i = 1:nOP
    utils.magplot(GFRF_w(i), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
utils.magplot(G_hat, 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Nominal G\_hat (' + nominal_type + ')');
plot(W_lft/2/pi, mag2db(Gmag + r), 'r--', 'DisplayName', '|G\_hat| \pm r')
plot(W_lft/2/pi, mag2db(max(Gmag - r, eps)), 'r--', 'HandleVisibility', 'off')
plot(W_lft/2/pi, mag2db(r), 'm', 'DisplayName', 'radius r')
title('Measured FRFs, nominal model and disk set')

% Nyquist: operating points and disks at a few frequencies
figure; hold on; axis equal; grid on
plot(real(squeeze(G_hat.ResponseData)), imag(squeeze(G_hat.ResponseData)), 'b', 'DisplayName', 'G\_hat')
th = linspace(0, 2*pi, 100);
for k = round(linspace(1, numel(W_lft), 8))
    plot(real(Z(k, :)), imag(Z(k, :)), 'k.', 'MarkerSize', 12, 'HandleVisibility', 'off')
    plot(real(G_hat.ResponseData(1, 1, k)) + r(k)*cos(th), ...
         imag(G_hat.ResponseData(1, 1, k)) + r(k)*sin(th), 'r', 'HandleVisibility', 'off')
end
xlabel('Re'); ylabel('Im'); title('Disk set (8 frequencies)')

%% Performance weight W1 (prewarped at omega_d so the discrete peak stays there)
W1_ct  = floor_gain + tf(peak_extra*[2*zeta_n*omega_d, 0], [1, 2*zeta_d*omega_d, omega_d^2]);
W1_frd = frd(c2d(W1_ct, Ts, c2dOptions('Method', 'tustin', 'PrewarpFrequency', omega_d)), W_lft);

%% Generalised plants
% performance: z = W1*(r - G u), y = r - G u
G_performance = mktito([W1_frd, -W1_frd*G_hat;
                        eye(ny), -G_hat], ny, nu);
% robustness: G = G_hat + delta, p = u, q = delta p (|delta| <= r), y = -q - G_hat u
G_robustness = mktito([zeros(ny*nu), eye(nu);
                       -eye(ny),     -G_hat], ny, nu);

%% IQC multiplier PI = diag(1, -1/r^2)  (elliptical set with A = I/r)
PI = frd(cat(3, reshape([ones(1, numel(r)); zeros(2, numel(r)); -1./r(:).'.^2], 2, 2, [])), W_lft, Ts);
PI = PI/norm(PI, inf);

stage_tags = {'rob', 'nperfhinf'};
for dir_ = {fullfile('Data', 'diaries'), fullfile('Data', 'controller')}
    if ~exist(dir_{1}, 'dir'), mkdir(dir_{1}); end
end

%% Synthesis
for orderK = orderK_list
    fprintf('\n===== Controller order %d =====\n', orderK);
    yalmip('clear');     % YALMIP numbering restarts, as in a fresh session
    rng(seed);

    K       = datadriven.Controller.SS(orderK, ny, nu, Ts, "right");
    K_stage = cell(2, 1);     % controller at the end of each stage
    K_itr   = {};

    % history of the solved iterations (to fall back on OL-stable controllers)
    hist_state  = cell(max_iter, 1);   % K.getstate() after the iteration
    hist_obj    = nan(max_iter, 1);
    hist_stable = false(max_iter, 1);  % controller open-loop stable
    hist_stage  = zeros(max_iter, 1);  % 0 = solver failed

    STAGE = 1;
    converged = false;
    nextstage = false;
    prev_obj  = Inf;          % for the convergence test
    iter_in_stage = 0;        % solved iterations in the current stage
    K_good_state  = K.getstate();   % last solved controller (rollback on failure)
    K_stage_state = cell(2, 1);     % parameters of K_stage{s} (for fallbacks)

    diaryfile = fullfile('Data', 'diaries', 'diary_temp.txt');
    if exist(diaryfile, 'file'), delete(diaryfile); end
    diary(diaryfile)

    for iter = 1:max_iter
        fprintf('Iter %3d | Stage %d\n', iter, STAGE);
        tic;
        switch STAGE
            case 1
                [LMIs, OBJ] = datadriven.Objective.IQC_reduced( ...
                    G_robustness, K, PI, alpha_PI, beta_PI, W_lft, ...
                    "addSlack", true, "epsNormalisation", 1e-6, "pinorm", true);
            case 2
                [LMIs, OBJ] = datadriven.Objective.Hinf(G_performance, K, W_lft);
                LMIs_rob = datadriven.Objective.IQC_reduced( ...
                    G_robustness, K, PI, alpha_PI, beta_PI, W_lft, ...
                    "addSlack", false, "pinorm", true);
                LMIs = [LMIs, LMIs_rob];
        end

        opts_sdp = sdpsettings('solver', solver, 'verbose', 0, ...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_NEAR_REL', mosek_near_rel);
        JOB = optimize(LMIs, OBJ, opts_sdp);
        elapsed = toc;

        K_itr{iter} = ss(K); %#ok<SAGROW>   raw solver output

        % ----- Solver failure: best OL-stable controller of the stage, then stop -----
        if JOB.problem ~= 0
            fprintf('  !! Solver failed (code %d: %s)\n', JOB.problem, JOB.info);
            [ok, idx] = select_stable(hist_stage, hist_stable, hist_obj, STAGE, 'best');
            if STAGE == 2 && ok
                K.setstate(hist_state{idx});
                fprintf('  --> back to the best OL stable controller of stage 2 (iter %d, obj = %.6f)\n', ...
                    idx, hist_obj(idx));
            else
                K.setstate(K_good_state);
                fprintf('  --> back to the last solved controller\n');
            end
            break
        end

        % ----- Solved iteration (oscillations of the objective are allowed) -----
        obj_val      = value(OBJ);
        K_good_state = K.getstate();
        K_is_stable  = isstable(ss(K));
        hist_state{iter}  = K_good_state;
        hist_obj(iter)    = obj_val;
        hist_stable(iter) = K_is_stable;
        hist_stage(iter)  = STAGE;
        iter_in_stage     = iter_in_stage + 1;
        if K_is_stable, stab = 'stable'; else, stab = 'unstable'; end
        fprintf('  obj = %.6f | K OL %s | time = %.2f s\n\n', obj_val, stab, elapsed);

        err = abs(prev_obj - obj_val);
        if STAGE == 2 && iter_in_stage >= max_iter_stage2
            fprintf('  --> Stage 2: %d iterations reached (max_iter_stage2)\n', iter_in_stage);
            nextstage = true;
        end
        if err < obj_err(STAGE) || err/prev_obj < obj_err(STAGE) || nextstage
            nextstage = false;
            K_stage{STAGE}       = ss(K);
            K_stage_state{STAGE} = K.getstate();
            fprintf('  --> Saved K_stage{%d}\n', STAGE);
            if STAGE < 2
                STAGE    = STAGE + 1;
                prev_obj = Inf;
                iter_in_stage = 0;
                fprintf('  --> Advancing to Stage %d\n', STAGE);
            else
                converged = true;
                fprintf('  --> Converged.\n');
                break
            end
        else
            prev_obj = obj_val;
        end
    end

    % ----- Final controller: best OL-stable stage-2 one if the last is unstable -----
    if STAGE == 2 && ~isstable(ss(K))
        [ok, idx] = select_stable(hist_stage, hist_stable, hist_obj, 2, 'best');
        if ok
            K.setstate(hist_state{idx});
            fprintf('Final controller OL unstable: best OL stable stage-2 controller (iter %d, obj = %.6f)\n', ...
                idx, hist_obj(idx));
        elseif ~isempty(K_stage_state{1}) && isstable(K_stage{1})
            K.setstate(K_stage_state{1});
            fprintf('No OL stable stage-2 controller: back to the stage-1 controller\n');
        else
            warning('No open-loop stable controller found: the final controller is OL unstable.');
        end
    end
    K_final = ss(K);
    K_stage{STAGE} = K_final;

    %% Nominal closed loop and save
    norminf_W1S = norm(W1_frd*frd(feedback(1, G_hat*K_final), W_lft), inf);
    fprintf('Order %d: ||W1 S||_inf = %.3f, K OL stable: %d\n', orderK, norminf_W1S, isstable(K_final));

    tag = stage_tags{STAGE};
    if ~converged, tag = [tag, '_notconv']; end %#ok<AGROW>
    expname = sprintf('%s_controller_order%d_%.0fhz_%.0fdb_%s', ...
        diskname, orderK, omega_d/2/pi, mag2db(min(abs(1./W1_frd.ResponseData))), tag);
    save(fullfile('Data', 'controller', sprintf('%s.mat', expname)), ...
        'K_final', 'K_stage', 'K_itr', 'hist_obj', 'hist_stable', 'hist_stage', ...
        'model_disk', 'W1_frd', 'orderK', 'idname', 'STAGE', 'converged', 'norminf_W1S')
    diary off
    movefile(diaryfile, fullfile('Data', 'diaries', sprintf('%s.txt', expname)))
end

%% Analysis (last synthesised order)
colors = lines(max(3, nOP));

% Sensitivity of the stage controllers on the nominal plant
figure; hold on
stage_names = {'Stage 1 (robust stab.)', 'Stage 2 (nom. H_\infty + robust stab.)'};
for s = 1:2
    if ~isempty(K_stage{s})
        utils.magplot(frd(feedback(1, G_hat*K_stage{s}), W_lft), 'Color', colors(s, :), ...
            'DisplayName', stage_names{s});
    end
end
utils.magplot(inv(W1_frd), 'Color', [0 0 1 .3], 'DisplayName', 'W_1^{-1}');
title(sprintf('Nominal sensitivity (disk set), order %d', orderK))

% Final controller on every measured operating point
figure; hold on
for i = 1:nOP
    utils.magplot(inv(1 + K_final*GFRF(i)), 'Color', colors(i, :), 'DisplayName', 'S OP'+string(i));
end
utils.magplot(inv(1 + K_final*G_hat), 'Color', 'k', 'DisplayName', 'S nominal');
utils.magplot(inv(W1_frd), 'Color', 'c', 'DisplayName', 'W_1^{-1}');
title(sprintf('Disk set: sensitivity on the measured FRFs, order %d', orderK))

figure; utils.magplot(K_final, 'DisplayName', 'K final'); title('Controller')
figure; pzplot(K_final); title('Controller poles/zeros')

%% Local functions
function [ok, idx] = select_stable(hist_stage, hist_stable, hist_obj, stage, mode)
%SELECT_STABLE OL-stable solved iteration of a stage: 'best' = lowest objective,
%   'last' = most recent.
    cand = find(hist_stage == stage & hist_stable);
    ok   = ~isempty(cand);
    idx  = [];
    if ~ok, return, end
    if mode == "best"
        [~, k] = min(hist_obj(cand));
        idx = cand(k);
    else
        idx = cand(end);
    end
end

function c = min_disk(z)
%MIN_DISK Centre of the smallest disk holding the points z (complex), exact:
%   the optimal disk passes through 2 (diameter) or 3 (circumcircle) of them.
    z = z(:);   n = numel(z);
    c = z(1);   rbest = Inf;   tol = 1e-12;
    if n == 1, return, end
    for i = 1:n-1
        for j = i+1:n
            cc = (z(i) + z(j))/2;   rr = abs(z(i) - cc);
            if rr < rbest && all(abs(z - cc) <= rr*(1 + 1e-9) + tol)
                c = cc;   rbest = rr;
            end
            for k = j+1:n
                [cc, rr] = circum(z(i), z(j), z(k));
                if isfinite(rr) && rr < rbest && all(abs(z - cc) <= rr*(1 + 1e-9) + tol)
                    c = cc;   rbest = rr;
                end
            end
        end
    end
end

function [c, r] = circum(a, b, d)
%CIRCUM Circumcircle of three points of the complex plane (Inf if collinear).
    b = b - a;   d = d - a;
    den = 2*imag(conj(b)*d);
    if abs(den) < eps*max(abs([b; d]))^2
        c = NaN;   r = Inf;   return
    end
    c = -1i*(abs(b)^2*d - abs(d)^2*b)/den;
    r = abs(c);   c = c + a;
end
