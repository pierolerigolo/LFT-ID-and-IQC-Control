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
%   The objective may oscillate (no step is rejected). Open-loop (OL) stability
%   of the controller is tracked at every iteration:
%     - solver failure in stage 2/3 -> best (lowest objective) OL-stable
%       controller of that stage; stage 2 then goes on with stage 3, stage 3 stops
%     - end of stage 2 on an OL-unstable controller -> best OL-stable one of stage 2
%     - end of stage 3 on an OL-unstable controller -> LAST OL-stable one of
%       stage 3 (objective not compared)
%   The final eta is recomputed by bisection for the final controller.
%   In stage 3 the bisection is warm-started from the previous D(w), E(w).
%
%   Output: Data/controller/<idname>_controller_... .mat and the run diary.
%
close all; clc; clear; clear("yalmip");

%% Config
% idname        = 'LFT_GFRF_GdmGLTIW2opt_1028_etfe';          % file written by b_LFT2ndorderID (pointwise W2, no floor, no hpf)
% idname        = 'LFT_GFRF_GdmGLTIW2optfl40_1028_etfe';      % file written by b_LFT2ndorderID (pointwise W2, -40 dB floor)
% idname        = 'LFT_GFRF_GdmGLTIW2n25hpf100fl40_1028_etfe'; % file written by b_LFT2ndorderID (W2: high-pass + -40 dB floor)
% idname      = 'LFT_GFRF_GdmGLTIW2n25hpf100_1028_etfe';     % high-pass, no floor (w2_floor_db = -Inf in b)
% idname      = 'LFT_GFRF_GdmGLTIW2n25_1028_etfe';         % without high-pass (use_hpf = false in b)
% idname = 'LFT_GFRF_GdmGLTIW2opthpf75fl40_1028_etfe';
% idname ='LFT_GFRF_GdmGLTIW2opthpf75fl40_1028_etfe_50hz';
idname='LFT_GFRF_GdmGLTIW2opt_1028_etfe_80hz_fmin10';
orderK_list   = 18;     % controller orders to synthesise
Ts            = 0.001;
seed          = 22;

% Performance weight W1 = floor_gain + peak_extra * band-pass at omega_d
omega_d       = 80*2*pi;   % disturbance frequency (rad/s)
floor_gain    = 0.5;
peak_extra    = 500;
zeta_d        = 0.001;
zeta_n        = 0.0001;

% % Performance weight W1 = floor_gain + peak_extra * band-pass at omega_d
% omega_d       = 50*2*pi;   % disturbance frequency (rad/s)
% floor_gain    = 0.5;
% peak_extra    = 500;
% zeta_d        = 0.001;
% zeta_n        = 0.0001;

% IQC
alpha_PI        = 0;         % split of the multiplier (Gupta et al. 2024, Cor. 1)
beta_PI         = 0.5;
multiplier_type = "frequency";   % "frequency" or "constant" D, E

% Optimisation
solver          = 'mosek';
mosek_near_rel  = 1e6;       % controller steps only (MOSEK default 1e3); the eta
                             % bisection always runs at the MOSEK defaults
max_iter        = 30;
max_iter_stage2 = 20;        % stage 2 ends after this many iterations even if the
                             % objective still oscillates (Inf = only on convergence)
obj_err         = [5e-3, 1e-3, 1e-3];   % convergence tolerance per stage

% Solver paths (machine specific)
% addpath(genpath("C:\Users\admin\mosek"))
% addpath(genpath("C:\Program Files\Mosek\10.1"))
% addpath(genpath("C:\Users\admin\Documents\MATLAB\tbxmanager"))
addpath(genpath("/Users/pierre/mosek"))
addpath(genpath("/Users/pierre/Documents/MATLAB/tbxmanager"))

%% Load the LFT plant
load(fullfile('Data', sprintf('%s.mat', idname)), 'GLFT', 'Delta', 'dm_bound', 'GFRF')
eta0            = dm_bound/1.5;      % initial bound on dm (mass units)
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

    % history of the solved iterations (to fall back on OL-stable controllers)
    hist_state  = cell(max_iter, 1);   % K.getstate() after the iteration
    hist_obj    = nan(max_iter, 1);
    hist_stable = false(max_iter, 1);  % controller open-loop stable
    hist_stage  = zeros(max_iter, 1);  % 0 = solver failed

    PI    = PI0;
    eta   = eta0;
    D_val = []; E_val = [];
    STAGE = 1;
    converged = false;
    nextstage = false;
    prev_obj  = Inf;          % for the convergence test
    iter_in_stage = 0;        % solved iterations in the current stage
    K_good_state  = K.getstate();   % last solved controller (rollback on failure)
    K_stage_state = cell(3, 1);     % parameters of K_stage{s} (for fallbacks)

    diaryfile = fullfile('Data', 'diaries', 'diary_temp.txt');
    if exist(diaryfile, 'file'),diary off; delete(diaryfile); end
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
                % if iter==9, nextstage=true;end
            case 3
                [LMIs, OBJ] = datadriven.Objective.Hinf(G_performance, K, W_lft);
                % (a) multiplier for the current (frozen) controller,
                %     warm-started from the previous D(w), E(w)
                [eta, D_val, E_val, PI_matrix] = multiplier_bisection(G_robustness, PILTI, K, eta, D_val, E_val, ...
                    alpha_PI, beta_PI, W_lft, dm_bound, multiplier_type);
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

        K_itr{iter}   = ss(K);    % raw solver output
        D_itr{iter}   = D_val;
        E_itr{iter}   = E_val;
        eta_itr{iter} = eta;      % multiplier bound used in this K-step

        % ----- Solver failure: best OL-stable controller of the stage -----
        if JOB.problem ~= 0
            fprintf('  !! Solver failed (code %d: %s)\n', JOB.problem, JOB.info);
            if STAGE == 1
                K.setstate(K_good_state);
                break
            end
            [ok, idx] = select_stable(hist_stage, hist_stable, hist_obj, STAGE, 'best');
            if ok
                K.setstate(hist_state{idx});
                eta = eta_itr{idx};
                fprintf('  --> back to the best OL stable controller of stage %d (iter %d, obj = %.6f)\n', ...
                    STAGE, idx, hist_obj(idx));
            elseif STAGE == 2 && ~isempty(K_stage_state{1}) && isstable(K_stage{1})
                K.setstate(K_stage_state{1});
                fprintf('  --> no OL stable controller in stage 2: back to the stage-1 controller\n');
            else
                K.setstate(K_good_state);
                fprintf('  --> no OL stable controller in stage %d: back to the last solved one\n', STAGE);
            end
            K_good_state = K.getstate();
            if STAGE == 3, break, end
            K_stage{2}       = ss(K);   % stage 2: end the stage, go on with stage 3
            K_stage_state{2} = K.getstate();
            STAGE      = 3;
            prev_obj   = Inf;
            fprintf('  --> Advancing to Stage 3\n');
            continue
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
        fprintf('  obj = %.6f | eta = %.4g | K OL %s | time = %.2f s\n\n', obj_val, eta, stab, elapsed);

        err = abs(prev_obj - obj_val);
        if STAGE == 2 && iter_in_stage >= max_iter_stage2
            fprintf('  --> Stage 2: %d iterations reached (max_iter_stage2)\n', iter_in_stage);
            nextstage = true;
        end
        if err < obj_err(STAGE) || err/prev_obj < obj_err(STAGE) || nextstage
            nextstage = false;
            if STAGE == 2 && ~K_is_stable
                % Stage 2 -> 3 from the best (lowest objective) OL-stable controller
                [ok, idx] = select_stable(hist_stage, hist_stable, hist_obj, 2, 'best');
                if ok
                    K.setstate(hist_state{idx});
                    fprintf('  --> Stage 2 ended on an OL unstable controller: best OL stable one (iter %d, obj = %.6f)\n', ...
                        idx, hist_obj(idx));
                elseif ~isempty(K_stage_state{1}) && isstable(K_stage{1})
                    K.setstate(K_stage_state{1});
                    fprintf('  --> No OL stable controller in stage 2: back to the stage-1 controller\n');
                else
                    warning('No open-loop stable controller found in stage 2: stage 3 starts from an unstable one.');
                end
                K_good_state = K.getstate();
            end
            K_stage{STAGE}       = ss(K);
            K_stage_state{STAGE} = K.getstate();
            fprintf('  --> Saved K_stage{%d}\n', STAGE);
            if STAGE < 3
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

    % ----- Final controller -----
    % Run ended in stage 2 (max_iter reached): best OL-stable stage-2 controller.
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
    % Stage 3: the last OL-stable stage-3 controller (objective NOT compared,
    % it may be higher than earlier ones). Its eta is recomputed below.
    if STAGE == 3 && ~isstable(ss(K))
        [ok, idx] = select_stable(hist_stage, hist_stable, hist_obj, 3, 'last');
        if ok
            K.setstate(hist_state{idx});
            fprintf('Final controller OL unstable: last OL stable stage-3 controller (iter %d)\n', idx);
        elseif ~isempty(K_stage_state{2}) && isstable(K_stage{2})
            K.setstate(K_stage_state{2});
            fprintf('No OL stable stage-3 controller: back to the controller stage 3 started from\n');
        else
            warning('No open-loop stable controller found: the final controller is OL unstable.');
        end
    end
    % Certified eta of the final controller (bisection at default MOSEK tolerances)
    if STAGE >= 2
        try
            [eta, D_val, E_val] = multiplier_bisection(G_robustness, PILTI, K, eta, D_val, E_val, ...
                alpha_PI, beta_PI, W_lft, dm_bound, multiplier_type);
        catch ME
            warning('Final eta bisection failed (%s): eta kept at %.4g.', ME.message, eta);
        end
    end
    K_final = ss(K);
    K_stage{STAGE} = K_final;

    %% Nominal closed loop and save
    norminf_W1S = norm(W1_frd*frd(feedback(1, G_hat*K_final), W_lft), inf);
    fprintf('Order %d: ||W1 S||_inf = %.3f, certified eta = %.4g (dm_bound = %.4g), K OL stable: %d\n', ...
        orderK, norminf_W1S, eta, dm_bound, isstable(K_final));

    tag = stage_tags{STAGE};
    if ~converged, tag = [tag, '_notconv']; end %#ok<AGROW>
    expname = sprintf('%s_controller_order%d_%.0fhz_%.0fdb_%s_etai%.2f_etaf%.2f_multi_%s', ...
        idname, orderK, omega_d/2/pi, mag2db(min(abs(1./W1_frd.ResponseData))), ...
        tag, eta0, eta, extractBefore(multiplier_type, 5));
    save(fullfile('Data', 'controller', sprintf('%s.mat', expname)), ...
        'K_final', 'K_stage', 'K_itr', 'D_itr', 'E_itr', 'eta_itr', ...
        'hist_obj', 'hist_stable', 'hist_stage', 'D_val', 'E_val', ...
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


%% 80hz chosen controller on every measured operating point vs SOA
% opensave=load(fullfile('data','controller','LFT_GFRF_GdmGLTIW2n25hpf100fl40_1028_etfe_controller_order18_80hz_-34db_nperfhinf_notconv_etai0.18_etaf0.42_multi_freq.mat'))
% iteration=18;
% K_chosen=opensave.K_final;
% opensave=load(fullfile('data','controller','LFT_GFRF_GdmGLTIW2n25hpf100_1028_etfe_controller_order18_80hz_-34db_nperfoptiIQC_etai0.18_etaf0.88_multi_freq.mat'));
% iteration=26;
% K_chosen=opensave.K_itr{iteration};
% opensave=load(fullfile('data','controller','LFT_GFRF_GdmGLTIW2n25hpf100fl40_1028_etfe_controller_order16_80hz_-34db_nperfhinf_notconv_etai0.18_etaf0.33_multi_freq.mat'));
% K_chosen=opensave.K_final;
opensave=load(fullfile('data','controller','LFT_GFRF_GdmGLTIW2opthpf75fl40_1028_etfe_controller_order18_80hz_-34db_nperfoptiIQC_notconv_etai0.24_etaf0.37_multi_freq.mat'));
K_chosen=opensave.K_final;
SOAsave=load(fullfile('data','finalrepportcont','controller_ellips_sus_7nm_order20_82hz_-30db_nomperf.mat'));
KSOA=SOAsave.K_itr{end};
if isstable(K_chosen), fprintf('K chosen is OL stable\n'); else, fprintf('K chosen is OL unstable\n'); end
figure; hold on
GFRF=interp(GFRF,W_lft)
for i = 1:nOP
    utils.magplot(inv(1 + K_chosen*GFRF(i)), 'Color', colors(i, :), 'DisplayName', 'S OP'+string(i));
    utils.magplot(inv(1 + KSOA*GFRF(i)), 'Color', [colors(i, :) .3],'linestyle','--', 'DisplayName', 'S SOA OP'+string(i));
end
utils.magplot(inv(1 + K_chosen*G_hat), 'Color', 'k', 'DisplayName', 'S nominal');
utils.magplot(inv(1 + KSOA*G_hat), 'Color', [0 0 0 .3],'linestyle','--', 'DisplayName', 'S SOA nominal');
utils.magplot(inv(opensave.W1_frd), 'Color', 'c', 'DisplayName', 'W_1^{-1}');
title(sprintf('Sensitivity on the measured FRFs, order %d, \\eta = %.3g (dm_{bound}=%.3f), compared to ellipses SOA', opensave.orderK, opensave.eta,opensave.dm_bound))

figure; utils.magplot(K_chosen, 'DisplayName', 'K chosen')
utils.magplot(KSOA, 'DisplayName', 'K SOA');title('Controller')
%%
L0 = K_chosen*GFRF(1);                                 % 0 kg
[mm, k] = min(abs(squeeze(1 + L0.ResponseData)));
fprintf('min |1+L| = %.3f at %.1f Hz\n', mm, L0.Frequency(k)/2/pi);
figure; nyquist(L0)
%% 50hz chosen controller on every measured operating point
% % -30db
% opensave_A=load(fullfile('data','controller','LFT_GFRF_GdmGLTIW2opthpf75fl40_1028_etfe_50hz_controller_order16_50hz_-30db_nperfoptiIQC_etai0.24_etaf0.82_multi_freq.mat'));
% K_A=opensave_A.K_final;
% opensave_B=load(fullfile('data','controller','LFT_GFRF_GdmGLTIW2opthpf75fl40_1028_etfe_50hz_controller_order18_50hz_-30db_nperfoptiIQC_etai0.24_etaf1.63_multi_freq.mat'));
% K_B=opensave_B.K_final;

% -34db
fn_A='LFT_GFRF_GdmGLTIW2opthpf75fl40_1028_etfe_50hz_controller_order20_50hz_-34db_nperfoptiIQC_etai0.24_etaf1.72_multi_freq';
opensave_A=load(fullfile('data','controller',[fn_A,'.mat']));
K_A=opensave_A.K_final;
% fn_B='LFT_GFRF_GdmGLTIW2opthpf75fl40_1028_etfe_50hz_controller_order16_50hz_-34db_nperfoptiIQC_etai0.24_etaf1.50_multi_freq';
fn_B='LFT_GFRF_GdmGLTIW2opthpf75fl40_1028_etfe_50hz_controller_order22_50hz_-34db_nperfoptiIQC_notconv_etai0.24_etaf0.80_multi_freq'
opensave_B=load(fullfile('data','controller',[fn_B,'.mat']));
K_B=opensave_B.K_final;

if isstable(K_A), fprintf('K A is OL stable\n'); else, fprintf('K A is OL unstable\n'); end
if isstable(K_B), fprintf('K B is OL stable\n'); else, fprintf('K B is OL unstable\n'); end
figure; hold on
GFRF=interp(GFRF,W_lft);
for i = 1:nOP
    utils.magplot(inv(1 + K_A*GFRF(i)), 'Color', colors(i, :), 'DisplayName', ...
        'S A OP'+string(i));
    utils.magplot(inv(1 + K_B*GFRF(i)), 'Color', [colors(i, :) .3],'linestyle','--', ...
        'DisplayName', 'S B OP'+string(i));
end
utils.magplot(inv(1 + K_A*G_hat), 'Color', 'k', 'DisplayName', 'S B nominal');
utils.magplot(inv(1 + K_B*G_hat), 'Color', [0 0 0 .3],'linestyle','--', 'DisplayName', ...
    'S A nominal');
utils.magplot(inv(opensave_A.W1_frd), 'Color', 'c', 'DisplayName', 'W_1^{-1}');
title(sprintf(['Sensitivity on the measured FRFs, K A order %d, \\eta = %.3g' ...
    ' compared to K B order %d, \\eta = %.3g (dm_{bound}=%.3f)'], opensave_A.orderK, ...
    opensave_A.eta,opensave_B.orderK, opensave_B.eta,opensave_A.dm_bound))

figure; utils.magplot(K_A, 'DisplayName', 'K A')
utils.magplot(K_B, 'DisplayName', 'K B');title('Controller')


%% Local functions
function [eta, D_val, E_val, PI_matrix] = multiplier_bisection(G_robustness, PILTI, K, eta, D_val, E_val, ...
        alpha_PI, beta_PI, W_lft, dm_bound, multiplier_type)
%MULTIPLIER_BISECTION Largest eta (and D(w), E(w)) certified for the frozen controller K,
%   warm-started from the previous D(w), E(w) when available.
    warm = {};
    if ~isempty(D_val), warm = {"D0", D_val, "E0", E_val}; end
    [eta, D_val, E_val, PI_matrix] = ...
        datadriven.Objective.IQC_reduced_combinedPI_par_fix_eta_bisection( ...
        G_robustness, PILTI, datadriven.Controller.Fixed(K), alpha_PI, beta_PI, W_lft, ...
        "etaUpperInit", eta, "etaLower", 1e-2, "etaHardCap", 100*dm_bound, ...
        "multiplier", multiplier_type, warm{:});
end

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
