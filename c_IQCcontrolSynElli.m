%% c_IQCcontrolSynElli  -  Step 3 (alternative): robust controller with an elliptical uncertainty set
%
%   Same approach as exp_controller.m of
%     V. Gupta, E. Klauser, A. Karimi, "Non-parametric IQC multipliers in
%     data-driven robust controller synthesis", 2024
%     (gitlab.epfl.ch/ddmac/paper-codes/2024/non-parametric-iqc-for-generalised-plant)
%
%   Uncertainty: directly on the measured FRFs (no mass block, no W2)
%       G(jw) in { G_hat(jw) + delta : || A(jw) [Re delta; Im delta] || <= 1 }
%   The nominal model G_hat and the ellipses A(jw) are computed together here
%   (minimum-volume ellipse around the operating points at every frequency).
%   Performance: nominal ||W1 S||_inf, same W1 and frequency grid as
%   c_IQCcontrolSynW2 (for comparison).
%
%   Stages (as in the reference code):
%     0  performance only                      (H-inf on the nominal plant)
%     1  robust stability only                 (IQC slack, until slack < tol)
%     2  nominal performance + robust stability (H-inf s.t. IQC)
%   OL-stability fallback: if the run ends on an OL-unstable controller, the
%   best (lowest objective) OL-stable stage-2 controller is used.
%
%   Output: Data/controller/ELLI_<expname>_controller_... .mat and the run diary.
%
close all; clc; clear; clear("yalmip");

%% Config
idname      = 'LFT_GFRF_GdmGLTIW2opt_1028_etfe';  % only its frequency grid and measured FRFs are used
orderK_list = 10:2:20;       % controller orders to synthesise
Ts          = 0.001;
seed        = 22;

% Performance weight W1 = floor_gain + peak_extra * band-pass at omega_d
omega_d     = 80*2*pi;       % disturbance frequency (rad/s)
floor_gain  = 0.5;
peak_extra  = 500;
zeta_d      = 0.001;
zeta_n      = 0.0001;

% IQC (values of the reference code)
alpha_PI    = 0;
beta_PI     = 0;

% Optimisation
solver          = 'mosek';
mosek_near_rel  = 1e3;               % MOSEK default
max_iter        = 100;
max_iter_stage  = [20, 30, 50];      % iterations allowed per stage 0, 1, 2
tol_stage       = [1e-4, 1e-4, 1e-4];% stage 0/2: rel. change of sqrt(obj); stage 1: slack

% Solver paths (machine specific)
addpath(genpath("C:\Users\admin\mosek"))
addpath(genpath("C:\Program Files\Mosek\10.1"))
addpath(genpath("C:\Users\admin\Documents\MATLAB\tbxmanager"))
% addpath(genpath("/Users/pierre/mosek"))
% addpath(genpath("/Users/pierre/Documents/MATLAB/tbxmanager"))

%% Data and frequency grid (same as the W2 design)
load(fullfile('Data', sprintf('%s.mat', idname)), 'GLFT', 'GFRF')
W_lft  = GLFT.Frequency;
nOP    = size(GFRF, 2);
GFRF_w = interp(GFRF, W_lft);
expname = regexprep(idname, '^.*?_(\d{4}_.*)$', '$1');   % e.g. 1028_etfe

%% Elliptical uncertainty set and nominal model
clear frdmodels
for i = 1:nOP
    frdmodels(1, 1, i) = GFRF_w(i); %#ok<SAGROW>   1 x 1 x nOP model array
end
model_elliptical = datadriven.find_elliptical_uncertainty(frdmodels, W_lft, 'verbosity', 1);
G_hat = model_elliptical.NominalModel;
A_set = model_elliptical.UncertaintySet.A{1};          % one 2x2 matrix per frequency

bad = cellfun(@(x) any(~isfinite(x(:))), A_set) | ~isfinite(squeeze(G_hat.ResponseData));
if any(bad)
    error('Elliptical set / nominal model not finite at %s Hz.', mat2str(W_lft(bad).'/2/pi, 4));
end
[ny, nu] = size(G_hat);

% Plot: FRFs, nominal model and |G_hat| +/- largest semi-axis of the ellipse
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

%% Performance weight W1 (prewarped at omega_d)
W1_ct  = floor_gain + tf(peak_extra*[2*zeta_n*omega_d, 0], [1, 2*zeta_d*omega_d, omega_d^2]);
W1_frd = frd(c2d(W1_ct, Ts, c2dOptions('Method', 'tustin', 'PrewarpFrequency', omega_d)), W_lft);

%% Generalised plants (reference code)
% performance: z = W1*(r - G u), y = r - G u
G_performance = mktito([W1_frd,  -W1_frd*G_hat;
                        eye(ny), -G_hat], ny, nu);

% robustness: p = R u, y = -sum_block*L*q - G_hat u, q = [Re delta; Im delta]-channel
L = [kron(kron(eye(ny), ones(1, nu)), [1, 0]);
     kron(kron(eye(ny), ones(1, nu)), [0, 1])];
R = kron(ones(ny, 1), eye(nu));
sum_block = kron(ones(1, 2), eye(ny));
G_robustness = mktito([zeros(ny*nu, ny*nu*2), R;
                       -sum_block*L,          -G_hat], ny, nu);

%% IQC multiplier PI = diag(I, -A_bar'*A_bar), A_bar = A*J'
J   = [1, 0; 0, 1i];
tmp = model_elliptical.UncertaintySet.A';    % row-wise expansion
tmp = [tmp{:}];
tmp = cellfun(@(x) x*J', tmp, "UniformOutput", 0);
A_bar = splitapply(@(x) {blkdiag(x{:})}, tmp, (1:length(W_lft))');

PI_cell = cellfun(@(x) blkdiag(eye(ny*nu), -x'*x), A_bar, "UniformOutput", 0);
PI = frd(cat(3, PI_cell{:}), W_lft, Ts);
PI = PI/norm(PI, inf);

stage_names = {'perf', 'rob', 'nperf'};
for dir_ = {fullfile('Data', 'diaries'), fullfile('Data', 'controller')}
    if ~exist(dir_{1}, 'dir'), mkdir(dir_{1}); end
end

%% Synthesis
for orderK = orderK_list
    fprintf('\n===== Controller order %d =====\n', orderK);
    yalmip('clear');
    rng(seed);

    K = datadriven.Controller.SS(orderK, ny, nu, Ts, "right");
    K_stage = cell(3, 1);     % controller at the end of stages 0, 1, 2
    K_itr   = {};
    hist_state  = cell(max_iter, 1);
    hist_obj    = nan(max_iter, 1);
    hist_stable = false(max_iter, 1);
    hist_stage  = -ones(max_iter, 1);   % -1 = solver failed

    STAGE = 0;
    iter_in_stage = 0;
    prev_obj  = Inf;
    converged = false;
    K_good_state = K.getstate();

    diaryfile = fullfile('Data', 'diaries', 'diary_temp.txt');
    if exist(diaryfile, 'file'), delete(diaryfile); end
    diary(diaryfile)

    for iter = 1:max_iter
        fprintf('Iter %3d | Stage %d\n', iter, STAGE);
        tic;
        switch STAGE
            case 0   % performance only
                [LMIs, OBJ] = datadriven.Objective.Hinf(G_performance, K, W_lft);
            case 1   % robustness only (slack)
                [LMIs, OBJ] = datadriven.Objective.IQC_reduced( ...
                    G_robustness, K, PI, alpha_PI, beta_PI, W_lft, "epsNormalisation", 1e-4);
            case 2   % nominal performance with robustness
                [LMIs, OBJ] = datadriven.Objective.Hinf(G_performance, K, W_lft);
                LMIs_rob = datadriven.Objective.IQC_reduced( ...
                    G_robustness, K, PI, alpha_PI, beta_PI, W_lft, ...
                    "addSlack", false, "epsNormalisation", 1e-4);
                LMIs = [LMIs, LMIs_rob];
        end

        opts_sdp = sdpsettings('solver', solver, 'verbose', 0, ...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_NEAR_REL', mosek_near_rel);
        JOB = optimize(LMIs, OBJ, opts_sdp);
        elapsed = toc;

        if JOB.problem ~= 0
            fprintf('  !! Solver failed (code %d: %s) -> last solved controller kept\n', JOB.problem, JOB.info);
            K.setstate(K_good_state);
            break
        end

        K.normalise();            % balanced realisation, as in the reference code
        K_good_state = K.getstate();
        K_itr{iter}  = ss(K); %#ok<SAGROW>
        if STAGE == 1
            obj_val = double(OBJ);          % slack
        else
            obj_val = sqrt(double(OBJ));    % ||W1 S||_inf on the grid
        end
        K_is_stable = isstable(K_itr{iter});
        hist_state{iter}  = K_good_state;
        hist_obj(iter)    = obj_val;
        hist_stable(iter) = K_is_stable;
        hist_stage(iter)  = STAGE;
        iter_in_stage     = iter_in_stage + 1;
        if K_is_stable, stab = 'stable'; else, stab = 'unstable'; end
        fprintf('  obj = %.6f | K OL %s | time = %.2f s\n\n', obj_val, stab, elapsed);

        % ----- Stage transitions -----
        err = abs(obj_val - prev_obj);
        switch STAGE
            case {0, 2}, done = err < tol_stage(STAGE+1) || err/prev_obj < tol_stage(STAGE+1);
            case 1,      done = abs(obj_val) < tol_stage(2);
        end
        if ~done && iter_in_stage >= max_iter_stage(STAGE+1)
            fprintf('  --> Stage %d: %d iterations reached\n', STAGE, iter_in_stage);
            if STAGE == 1
                warning('Robust stability not reached in stage 1 (slack = %.3g): stopping.', obj_val);
                K_stage{2} = ss(K);
                break
            end
            done = true;
        end
        if done
            K_stage{STAGE+1} = ss(K);
            fprintf('  --> Saved K_stage{%d} (%s)\n', STAGE+1, stage_names{STAGE+1});
            if STAGE < 2
                STAGE = STAGE + 1;
                iter_in_stage = 0;
                prev_obj = Inf;
                fprintf('  --> Advancing to Stage %d\n', STAGE);
                continue
            end
            converged = true;
            fprintf('  --> Converged.\n');
            break
        end
        prev_obj = obj_val;
    end

    % ----- Final controller: best OL-stable stage-2 one if the last is unstable -----
    if STAGE == 2 && ~isstable(ss(K))
        cand = find(hist_stage == 2 & hist_stable);
        if ~isempty(cand)
            [~, k] = min(hist_obj(cand));
            K.setstate(hist_state{cand(k)});
            fprintf('Final controller OL unstable: best OL stable stage-2 controller (iter %d, obj = %.4f)\n', ...
                cand(k), hist_obj(cand(k)));
        else
            warning('No open-loop stable stage-2 controller: the final controller is OL unstable.');
        end
    end
    K_final = ss(K);
    if STAGE == 2, K_stage{3} = K_final; end

    norminf_W1S = norm(W1_frd*frd(feedback(1, G_hat*K_final), W_lft), inf);
    fprintf('Order %d: ||W1 S||_inf (nominal) = %.3f, K OL stable: %d, robust (stage 2 reached): %d\n', ...
        orderK, norminf_W1S, isstable(K_final), STAGE == 2);

    tag = stage_names{STAGE+1};
    if ~converged, tag = [tag, '_notconv']; end %#ok<AGROW>
    outname = sprintf('ELLI_%s_controller_order%d_%.0fhz_%.0fdb_%s', expname, orderK, ...
        omega_d/2/pi, mag2db(min(abs(1./W1_frd.ResponseData))), tag);
    save(fullfile('Data', 'controller', sprintf('%s.mat', outname)), ...
        'K_final', 'K_stage', 'K_itr', 'hist_obj', 'hist_stable', 'hist_stage', ...
        'model_elliptical', 'W1_frd', 'orderK', 'idname', 'STAGE', 'converged', 'norminf_W1S')
    diary off
    movefile(diaryfile, fullfile('Data', 'diaries', sprintf('%s.txt', outname)))
end

%% Analysis (last synthesised order)
colors = lines(max(3, nOP));

figure; hold on
for i = 1:nOP
    utils.magplot(inv(1 + K_final*GFRF(i)), 'Color', colors(i, :), 'DisplayName', 'S OP'+string(i));
end
utils.magplot(inv(1 + K_final*G_hat), 'Color', 'k', 'DisplayName', 'S nominal');
utils.magplot(inv(W1_frd), 'Color', 'c', 'DisplayName', 'W_1^{-1}');
title(sprintf('Elliptical set: sensitivity on the measured FRFs, order %d', orderK))

figure; utils.magplot(K_final, 'DisplayName', 'K final'); title('Controller')
figure; pzplot(K_final); title('Controller poles/zeros')
