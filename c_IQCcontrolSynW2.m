%% c_IQCcontrolSynW2  -  Step 3: data-driven robust controller synthesis (IQC)
%
%   Plant: LFT from b_LFT2ndorderID,  G = G_dm(dm) * G_LTI * (1 + W2*Delta)
%     channel 1: Delta (W2 block, 1x1 LTI, ||Delta|| <= 1)  -> multiplier diag(1,-1)
%     channel 2: dm    (real, centred, NOT normalised)       -> multiplier
%                [eta^2*D, E; E', -D], eta in mass units (eta = dm_bound
%                covers the whole identified mass range)
%   Performance: nominal ||W1 S||_inf, W1 = floor + narrow peak at omega_d.
%
%   Stages (frequency-domain convex-concave iterations, fixed-order SS controller):
%     1  robust stability only, constant multiplier at eta0 (slack minimised)
%     2  min nominal H-inf  s.t. robust stability with the multiplier of stage 1
%     3  per iteration: (a) largest eta and optimised multiplier D(w), E(w) for
%        the current controller (bisection), (b) as stage 2 with that multiplier.
%   A step that increases the objective is rejected (rollback) and ends the stage.
%   Stage 2 keeps iterating through open-loop unstable controllers, but stage 3
%   starts from the last open-loop stable controller of stage 2. In stage 3 the
%   bisection is warm-started from the previous D(w), E(w).
%
%   Output: Data/controller/<idname>_controller_... .mat and the run diary.
%
close all; clc; clear; clear("yalmip");

%% Config
idname        = 'LFT_GFRF_GdmGLTIW2n40_1028_spa';   % file written by b_LFT2ndorderID
orderK_list   = 16:2:20;     % controller orders to synthesise
Ts            = 0.001;
seed          = 22;

% Performance weight W1 = floor_gain + peak_extra * band-pass at omega_d
omega_d       = 80*2*pi;   % disturbance frequency (rad/s)
floor_gain    = 0.5;
peak_extra    = 500;
zeta_d        = 0.001;
zeta_n        = 0.0001;

% IQC
eta0            = 0.05;      % initial bound on dm (mass units)
alpha_PI        = 0;         % split of the multiplier (Gupta et al. 2024, Cor. 1)
beta_PI         = 0.5;
multiplier_type = "frequency";   % "frequency" or "constant" D, E

% Optimisation
solver          = 'mosek';
mosek_near_rel  = 1e3;       % MOSEK default; 1e6 accepted very inaccurate solutions
max_iter        = 50;
obj_err         = [5e-3, 1e-3, 1e-3];   % convergence tolerance per stage
obj_rise_tol    = 1e-3;      % relative objective increase that rejects a step

% Solver paths (machine specific)
addpath(genpath("C:\Users\admin\mosek"))
addpath(genpath("C:\Program Files\Mosek\10.1"))
addpath(genpath("C:\Users\admin\Documents\MATLAB\tbxmanager"))
% addpath(genpath("/Users/pierre/mosek"))
% addpath(genpath("/Users/pierre/Documents/MATLAB/tbxmanager"))

%% Load the LFT plant
load(fullfile('Data', sprintf('%s.mat', idname)), 'GLFT', 'Delta', 'dm_bound', 'GFRF')
W_lft = GLFT.Frequency;
[nout_lft, nin_lft] = size(GLFT);
[nq, np] = size(Delta);          % Delta: p (plant outputs, np) -> q (plant inputs, nq)
ny    = nout_lft - np;
nu    = nin_lft - nq;
G_hat = GLFT(np+1:end, nq+1:end);   % nominal plant (dm = 0, Delta = 0)
nOP   = size(GFRF, 2);

figure; hold on
for i = 1:nOP
    utils.magplot(GFRF(i), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
utils.magplot(G_hat, 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Nominal G\_hat');
title('Measured FRFs and nominal model')

%% Performance weight W1 (prewarped at omega_d so the discrete peak stays there)
W1_ct  = floor_gain + tf(peak_extra*[2*zeta_n*omega_d, 0], [1, 2*zeta_d*omega_d, omega_d^2]);
W1_frd = frd(c2d(W1_ct, Ts, c2dOptions('Method', 'tustin', 'PrewarpFrequency', omega_d)), W_lft);

%% Generalised plants
% performance: z = W1*(r - G u), y = r - G u
G_performance = mktito([W1_frd, -W1_frd*G_hat;
                        eye(ny), -G_hat], ny, nu);
% robustness: [p; y] = GLFT [q; u]
G_robustness  = mktito(GLFT, ny, nu);

%% IQC multipliers, ordering (p1, p2, q1, q2)
N     = length(W_lft);
PILTI = frd(repmat(diag([1, -1]), [1 1 N]), W_lft, Ts);          % W2 block
PI0   = frd(repmat(diag([1, eta0^2, -1, -1]), [1 1 N]), W_lft, Ts); % W2 block + dm at eta0

stage_tags = {'rob', 'nperfhinf', 'nperfoptiIQC'};
for dir_ = {fullfile('Data', 'diaries'), fullfile('Data', 'controller')}
    if ~exist(dir_{1}, 'dir'), mkdir(dir_{1}); end
end

%% Synthesis
for orderK = orderK_list
    fprintf('\n===== Controller order %d =====\n', orderK);
    yalmip('clear');     % YALMIP numbering restarts, as in a fresh session
    rng(seed);

    K       = datadriven.Controller.SS(orderK, ny, nu, Ts, "right");
    K_stage = cell(3, 1);     % controller at the end of each stage
    K_itr = {}; D_itr = {}; E_itr = {}; eta_itr = {};
    accepted_itr = false(max_iter, 1);

    PI    = PI0;
    eta   = eta0;
    D_val = []; E_val = [];
    STAGE = 1;
    converged = false;
    nextstage = false;
    prev_obj  = Inf;          % for the convergence test
    obj_good  = Inf;          % objective of the last accepted controller
    K_good_state   = K.getstate();
    K_stable_state = [];      % last accepted open-loop stable controller (stage 2)
    obj_stable     = Inf;     % its objective

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
            case 3
                [LMIs, OBJ] = datadriven.Objective.Hinf(G_performance, K, W_lft);
                % (a) multiplier for the current (frozen) controller,
                %     warm-started from the previous D(w), E(w)
                warm = {};
                if ~isempty(D_val), warm = {"D0", D_val, "E0", E_val}; end
                [eta, D_val, E_val, PI_matrix] = ...
                    datadriven.Objective.IQC_reduced_combinedPI_par_fix_eta_bisection( ...
                    G_robustness, PILTI, datadriven.Controller.Fixed(K), alpha_PI, beta_PI, W_lft, ...
                    "etaUpperInit", eta, "etaLower", 1e-2, "etaHardCap", 100*dm_bound, ...
                    "multiplier", multiplier_type, warm{:});
                PI = frd(PI_matrix, W_lft, Ts);
                % (b) controller step with that multiplier
                LMIs_rob = datadriven.Objective.IQC_reduced( ...
                    G_robustness, K, PI, alpha_PI, beta_PI, W_lft, ...
                    "addSlack", false, "pinorm", true);
                LMIs = [LMIs, LMIs_rob];
        end

        opts_sdp = sdpsettings('solver', solver, 'verbose', 0, ...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_NEAR_REL', mosek_near_rel);
        JOB = optimize(LMIs, OBJ, opts_sdp);
        elapsed = toc;

        K_itr{iter}   = ss(K);    % raw solver output (also kept when rejected)
        D_itr{iter}   = D_val;
        E_itr{iter}   = E_val;
        eta_itr{iter} = eta;

        solver_failed = JOB.problem ~= 0;
        if solver_failed
            fprintf('  !! Solver failed (code %d: %s) -> rolled back to last accepted controller\n', ...
                JOB.problem, JOB.info);
            K.setstate(K_good_state);
            if STAGE ~= 2, break, end   % in stage 2: end the stage, go on with stage 3
        end

        % The previous controller is feasible for the convexified problem of
        % this iteration (in stage 3 the multiplier was just certified for it),
        % so the objective cannot increase in exact arithmetic: an increase is
        % an inaccurate solve -> reject, roll back and end the stage.
        if solver_failed
            obj_val   = obj_good;
            nextstage = true;
        elseif isfinite(obj_good) && value(OBJ) > obj_good*(1 + obj_rise_tol) + 1e-9
            obj_val = value(OBJ);
            fprintf('  !! Objective increased (%.6f -> %.6f): step rejected, rolled back\n', ...
                obj_good, obj_val);
            K.setstate(K_good_state);
            obj_val   = obj_good;
            nextstage = true;
        else
            obj_val      = value(OBJ);
            K_good_state = K.getstate();
            obj_good     = obj_val;
            accepted_itr(iter) = true;
        end

        K_is_stable = isstable(ss(K));
        if K_is_stable, stab = 'stable'; else, stab = 'unstable'; end
        if STAGE == 2 && K_is_stable && accepted_itr(iter)
            K_stable_state = K.getstate();
            obj_stable     = obj_val;
        end
        fprintf('  obj = %.6f | eta = %.4g | K OL %s | time = %.2f s\n\n', obj_val, eta, stab, elapsed);

        err = abs(prev_obj - obj_val);
        if err < obj_err(STAGE) || err/prev_obj < obj_err(STAGE) || nextstage
            nextstage = false;
            % Stage 2 -> 3: continue from the last open-loop stable controller
            if STAGE == 2 && ~K_is_stable
                if ~isempty(K_stable_state)
                    K.setstate(K_stable_state);
                    K_good_state = K_stable_state;   % rollback point for stage 3
                    obj_good     = obj_stable;
                    fprintf('  --> Stage 2 ended on an OL unstable controller: back to the last OL stable one (obj = %.6f)\n', obj_stable);
                else
                    warning('No open-loop stable controller found in stage 2: stage 3 starts from an unstable one.');
                end
            end
            K_stage{STAGE} = ss(K);
            fprintf('  --> Saved K_stage{%d}\n', STAGE);
            if STAGE < 3
                STAGE    = STAGE + 1;
                prev_obj = Inf;
                if STAGE == 2   % objective changes (slack -> H-inf)
                    obj_good = Inf;
                    if K_is_stable    % stage-1 controller as fallback (H-inf objective unknown)
                        K_stable_state = K.getstate();
                        obj_stable     = Inf;
                    end
                end
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

    % K always holds the last accepted controller (rolled back otherwise);
    % it is certified robustly stable for |dm| <= eta.
    K_final = ss(K);

    %% Nominal closed loop and save
    norminf_W1S = norm(W1_frd*frd(feedback(1, G_hat*K_final), W_lft), inf);
    fprintf('Order %d: ||W1 S||_inf = %.3f, certified eta = %.4g (dm_bound = %.4g)\n', ...
        orderK, norminf_W1S, eta, dm_bound);

    tag = stage_tags{STAGE};
    if ~converged, tag = [tag, '_notconv']; end %#ok<AGROW>
    expname = sprintf('%s_controller_order%d_%.0fhz_%.0fdb_%s_etai%.2f_etaf%.2f_multi_%s', ...
        idname, orderK, omega_d/2/pi, mag2db(min(abs(1./W1_frd.ResponseData))), ...
        tag, eta0, eta, extractBefore(multiplier_type, 5));
    save(fullfile('Data', 'controller', sprintf('%s.mat', expname)), ...
        'K_final', 'K_stage', 'K_itr', 'accepted_itr', 'D_itr', 'E_itr', 'eta_itr', ...
        'eta', 'eta0', 'dm_bound', 'W1_frd', 'orderK', 'idname', 'STAGE', 'converged', 'norminf_W1S')
    diary off
    movefile(diaryfile, fullfile('Data', 'diaries', sprintf('%s.txt', expname)))
end

%% Analysis (last synthesised order)
colors = lines(max(3, nOP));

% Sensitivity of the stage controllers on the nominal plant
figure; hold on
stage_names = {'Stage 1 (robust stab.)', 'Stage 2 (nom. H_\infty, \eta_0)', 'Stage 3 (nom. H_\infty, opt. IQC)'};
for s = 1:3
    if ~isempty(K_stage{s})
        utils.magplot(frd(feedback(1, G_hat*K_stage{s}), W_lft), 'Color', colors(s, :), ...
            'DisplayName', stage_names{s});
    end
end
utils.magplot(inv(W1_frd), 'Color', [0 0 1 .3], 'DisplayName', 'W_1^{-1}');
title(sprintf('Nominal sensitivity, order %d', orderK))

% Final controller on every measured operating point
figure; hold on
for i = 1:nOP
    utils.magplot(inv(1 + K_final*GFRF(i)), 'Color', colors(i, :), 'DisplayName', 'S OP'+string(i));
end
utils.magplot(inv(1 + K_final*G_hat), 'Color', 'k', 'DisplayName', 'S nominal');
utils.magplot(inv(W1_frd), 'Color', 'c', 'DisplayName', 'W_1^{-1}');
title(sprintf('Sensitivity on the measured FRFs, order %d, \\eta = %.3g', orderK, eta))

figure; utils.magplot(K_final, 'DisplayName', 'K final'); title('Controller')
figure; pzplot(K_final); title('Controller poles/zeros')
