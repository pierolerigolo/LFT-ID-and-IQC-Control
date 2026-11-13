%% c_IQCcontrolSynElli  -  Step 3 (alternative): robust controller with an elliptical uncertainty set
%
%   Uncertainty directly on the measured FRFs (no mass block, no W2), as in
%   exp_controller.m of Gupta, Klauser, Karimi (2024)
%   (gitlab.epfl.ch/ddmac/paper-codes/2024/non-parametric-iqc-for-generalised-plant):
%       G(jw) in { G_hat(jw) + delta : || A(jw) [Re delta; Im delta] || <= 1 }
%   The nominal model G_hat and the ellipses A(jw) are computed together here
%   (minimum-volume ellipse around the operating points at every frequency).
%   Same frequency grid, measured FRFs and W1 as c_IQCcontrolSynW2.
%
%   Stages (frequency-domain convex-concave iterations, fixed-order SS controller):
%     1  robust stability only (slack minimised)
%     2  min nominal H-inf  s.t. robust stability
%   The objective may oscillate (no step is rejected). Open-loop (OL) stability
%   of the controller is tracked at every iteration:
%     - solver failure in stage 2 -> best (lowest objective) OL-stable stage-2 controller
%     - end of stage 2 (convergence, max_iter_stage2 or max_iter) -> best (lowest
%       objective) OL-stable stage-2 controller, even when the last one is stable
%       (stage-1 one as fallback)
%
%   Output: Data/controller/ELLI_<expname>_controller_... .mat and the run diary.
%
close all; clc; clear; clear("yalmip");

%% Config
idname        = 'LFT_GFRF_GdmGLTIW2opt_1028_etfe';  % only its frequency grid and measured FRFs are used
orderK_list   = 10:2:20;     % controller orders to synthesise
Ts            = 0.001;
seed          = 22;

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
elliname = ['ELLI_', regexprep(idname, '^.*?_(\d{4}_.*)$', '$1')];   % e.g. ELLI_1028_etfe

%% Elliptical uncertainty set and nominal model
clear frdmodels
for i = 1:nOP
    frdmodels(1, 1, i) = GFRF_w(i); %#ok<SAGROW>   1 x 1 x nOP model array
end
model_elliptical = datadriven.find_elliptical_uncertainty(frdmodels, W_lft, 'verbosity', 1);
G_hat = model_elliptical.NominalModel;
A_set = model_elliptical.UncertaintySet.A{1};          % one 2x2 matrix per frequency
[ny, nu] = size(G_hat);

bad = cellfun(@(x) any(~isfinite(x(:))), A_set) | ~isfinite(squeeze(G_hat.ResponseData));
if any(bad)
    error('Elliptical set / nominal model not finite at %s Hz.', mat2str(W_lft(bad).'/2/pi, 4));
end

% FRFs, nominal model and |G_hat| +/- largest semi-axis of the ellipse
r    = cellfun(@(x) 1/min(svd(x)), A_set);
Gmag = abs(squeeze(G_hat.ResponseData));
figure; hold on
for i = 1:nOP
    utils.magplot(GFRF_w(i), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
utils.magplot(G_hat, 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Nominal G\_hat');
plot(W_lft/2/pi, mag2db(Gmag + r), 'r--', 'DisplayName', 'Elliptical bound')
plot(W_lft/2/pi, mag2db(max(Gmag - r, eps)), 'r--', 'HandleVisibility', 'off')
title('Measured FRFs, nominal model and elliptical set')

%% Performance weight W1 (prewarped at omega_d so the discrete peak stays there)
W1_ct  = floor_gain + tf(peak_extra*[2*zeta_n*omega_d, 0], [1, 2*zeta_d*omega_d, omega_d^2]);
W1_frd = frd(c2d(W1_ct, Ts, c2dOptions('Method', 'tustin', 'PrewarpFrequency', omega_d)), W_lft);

%% Generalised plants
% performance: z = W1*(r - G u), y = r - G u
G_performance = mktito([W1_frd, -W1_frd*G_hat;
                        eye(ny), -G_hat], ny, nu);
% robustness (reference code): p = R u, y = -sum_block*L*q - G_hat u
L = [kron(kron(eye(ny), ones(1, nu)), [1, 0]);
     kron(kron(eye(ny), ones(1, nu)), [0, 1])];
R = kron(ones(ny, 1), eye(nu));
sum_block = kron(ones(1, 2), eye(ny));
G_robustness = mktito([zeros(ny*nu, ny*nu*2), R;
                       -sum_block*L,          -G_hat], ny, nu);

%% IQC multiplier PI = diag(I, -A_bar'*A_bar), A_bar = A*J' (reference code)
J   = [1, 0; 0, 1i];
tmp = model_elliptical.UncertaintySet.A';    % row-wise expansion
tmp = [tmp{:}];
tmp = cellfun(@(x) x*J', tmp, "UniformOutput", 0);
A_bar = splitapply(@(x) {blkdiag(x{:})}, tmp, (1:length(W_lft))');

PI_cell = cellfun(@(x) blkdiag(eye(ny*nu), -x'*x), A_bar, "UniformOutput", 0);
PI = frd(cat(3, PI_cell{:}), W_lft, Ts);
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

    % ----- Final controller: best (lowest objective) OL-stable stage-2 one,
    %       also when the last one is stable -----
    if STAGE == 2
        [ok, idx] = select_stable(hist_stage, hist_stable, hist_obj, 2, 'best');
        if ok
            K.setstate(hist_state{idx});
            fprintf('Final controller: best OL stable stage-2 controller (iter %d, obj = %.6f)\n', ...
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
        elliname, orderK, omega_d/2/pi, mag2db(min(abs(1./W1_frd.ResponseData))), tag);
    save(fullfile('Data', 'controller', sprintf('%s.mat', expname)), ...
        'K_final', 'K_stage', 'K_itr', 'hist_obj', 'hist_stable', 'hist_stage', ...
        'model_elliptical', 'W1_frd', 'orderK', 'idname', 'STAGE', 'converged', 'norminf_W1S')
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
title(sprintf('Nominal sensitivity (elliptical set), order %d', orderK))

% Final controller on every measured operating point
figure; hold on
for i = 1:nOP
    utils.magplot(inv(1 + K_final*GFRF(i)), 'Color', colors(i, :), 'DisplayName', 'S OP'+string(i));
end
utils.magplot(inv(1 + K_final*G_hat), 'Color', 'k', 'DisplayName', 'S nominal');
utils.magplot(inv(W1_frd), 'Color', 'c', 'DisplayName', 'W_1^{-1}');
title(sprintf('Elliptical set: sensitivity on the measured FRFs, order %d', orderK))

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
