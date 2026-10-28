close all; clc; clear; clear("yalmip");
% warning('off')
%% Config
rng(22);
orderK   = 18 ;
Ts       = 0.001;

%% Load data
idname='LFT_GFRF_param_n7_2907';
load(fullfile('Data',sprintf('%s.mat',idname)))
% load(fullfile('data','Old id method/','idm_suspended_2907.mat'))
% load(fullfile('data','Old id method/','sus_data_2907_lft_full.mat'))
% dm_bound=dm_max;
% % For windows
% addpath(genpath("C:\Users\admin\mosek"))
% addpath(genpath("C:\Program Files\Mosek\10.1"))
% addpath(genpath("C:\Users\admin\Documents\MATLAB\tbxmanager"))

% For mac
addpath(genpath("/Users/pierre/mosek"))
addpath(genpath("/Users/pierre/Documents/MATLAB/tbxmanager"))

W_lft=GLFT.Frequency;
nomegas=length(W_lft);
G_hat=GLFT(2,2);
ny=1;nu=1;np=1;nq=1;nw=1;nz=1;

%% plot the plant
figure; hold on; grid on;
xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)');
set(gca, 'XScale', 'log');

% Plot raw FRFs
for i = 1:length(GFRF)
    utils.magplot(GFRF(i), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
end

% Plot nominal model
utils.magplot(GLFT(2,2), 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Nominal G\_hat');

legend('Location','best');
%%
% figure;bode(GLFT(2,2),GFRF(2))
%% Performance weight
omega_d    = 81.7*2*pi;
omega_n    = 1000;
floor_gain = 0.5;

zeta_d     = 0.001;
zeta_n     = 0.0001;
peak_extra = 500;

num_bp     = peak_extra * [2*zeta_n*omega_d, 0];
den_bp     = [1, 2*zeta_d*omega_d, omega_d^2];
W1_peak_ct = tf(num_bp, den_bp);

%W1_sin_ct=tf([0 2*110/1000 0],[1, 2*zeta_d*omega_d, omega_d^2]);
W1_ct      = floor_gain + W1_peak_ct;
% W1_ct      = floor_gain + W1_sin_ct;

W1_frd         = frd(c2d(W1_ct, G_hat.Ts, 'tustin'),W_lft);
% W1_frd     = frd((W1(1,1)+1)/4, W_lft);

utils.magplot(inv(W1_frd),'DisplayName','1/W1');
% exportgraphics(gcf,'plot/rp_input.svg')
%% Two-port performance plant  (Stage 0)
G_performance = mktito([ ...
    W1_frd      , -W1_frd*G_hat ; ...
    eye(ny)     , -G_hat        ; ...
    ], ny, nu);

%% Two-port robustness plant  (Stage 1)
G_robustness = mktito(GLFT,ny, nu);

%% Three-port combined plant  (Stage 3 — IQC_performance)
G11 = G_robustness(G_robustness.OutputGroup.Y1, G_robustness.InputGroup.U1);
G12 = G_robustness(G_robustness.OutputGroup.Y1, G_robustness.InputGroup.U2);
G21 = G_robustness(G_robustness.OutputGroup.Y2, G_robustness.InputGroup.U1);
G22 = G_robustness(G_robustness.OutputGroup.Y2, G_robustness.InputGroup.U2);

O_pw = frd(zeros(np, nw, length(W_lft)), W_lft);
I_y  = frd(repmat(eye(ny), [1,1,length(W_lft)]), W_lft);

G_mat = [  G11          ,  O_pw    ,  G12          ; ...  % p row
          -W1_frd*G21   ,  W1_frd  , -W1_frd*G22   ; ...  % z row
          -G21          ,  I_y     , -G22           ];     % y row

G_combined = frd(G_mat, W_lft);
G_combined.InputGroup  = struct('U1', 1:nq,        'U2', nq+(1:nw),    'U3', nq+nw+(1:nu));
G_combined.OutputGroup = struct('Y1', 1:np,        'Y2', np+(1:nz),    'Y3', np+nz+(1:ny));

disp(G_combined.InputGroup)
disp(G_combined.OutputGroup)

%% IQC multiplier PI  (bouded real scalar)
eta      = dm_bound*1;   % bound on delta = dm/m0
% eta=0.2;
PI_data  = repmat([eta^2, 0; 0, -1], [1, 1, length(W_lft)]);
PI       = frd(PI_data, W_lft, G_hat.Ts);

alpha_PI_soft = 0;
beta_PI_soft  = 0.5;
% alpha_PI_hard = 0;
% beta_PI_hard  = 1e-4;

%% Controller parametrisation
K = datadriven.Controller.SS(orderK, ny, nu, Ts, "right");

%% Storage for intermediate controllers
K_stage = cell(4, 1);   % one per stage 0..3

%% Synthesis
%
%  Stage 0 — Performance only        (H∞)
%  Stage 1 — Robustness only         (IQC_reduced slack)
%  Stage 2 — Robust nominal          (H∞ + IQC_reduced slack, joint obj)
%  Stage 3 — Robust performance      (IQC_performance)
%
solver       = 'mosek';
max_iter     = 1000;
prev_obj     = Inf;
STAGE        = 1;gamma=10;
has_slack    = false;
loss_tol = Inf;
eta_max = Inf;%2*dm_bound;
mutli_type='frequency';
nextstage=0;
last_hinfstb_itr=Inf;
D_val=repmat({eye(1)}, nomegas, 1);
E_val=repmat({zeros(1)}, nomegas, 1);

if ~exist(fullfile('data','diaries'),'dir')
    mkdir(fullfile('data','diaries'))
end
if exist(fullfile('data','diaries','diary_temp.txt'),'file')
    delete(fullfile('data','diaries','diary_temp.txt')) % clear the diary
end
diary(fullfile('data','diaries','diary_temp.txt'))
diary on


for iter = 1:max_iter
    fprintf('Iter %3d | Stage %d\n', iter, STAGE);
    tic;

    switch STAGE
        case 0
            % ── Stage 0: Performance only (H∞ warm-start) ───────────────
            [LMIs, OBJ,gamma] = datadriven.Objective.Hinf( ...
                G_performance, K, W_lft);
            % [LMIs_con, ~] = datadriven.Constraint.Hinf( ...
            %     G_performance, K, W_lft);
            LMIs     = [LMIs];
            has_slack = false; slack_val = [];


        case 1
            % ── Stage 1: Robustness only (IQC_reduced slack) ─────────────
            [LMIs, OBJ] = datadriven.Objective.IQC_reduced( ...
                G_robustness, K, PI, alpha_PI_soft, beta_PI_soft, W_lft, ...
                "addSlack", true, "epsNormalisation", 1e-6,"pinorm",true);
            slack_val = OBJ;   % OBJ is the slack her
            has_slack = true;
            % nextstage=1;

        case 2
            % STAGE=3;
            % ── Stage 2: Robust nominal (H∞ + IQC_reduced, joint obj) ────
            [LMIs, OBJ_perf,gamma] = datadriven.Objective.Hinf( ...
                G_performance, K, W_lft);
            % [LMIs_con, ~] = datadriven.Constraint.Hinf( ...
            %     G_performance, K, W_lft);
            [LMIs_rob, slack] = datadriven.Objective.IQC_reduced( ...
                G_robustness, K, PI, alpha_PI_soft, beta_PI_soft, W_lft, ...
                "addSlack", false,"pinorm",true);
            LMIs      = [LMIs, LMIs_rob];
            OBJ       = OBJ_perf;
            slack_val = slack;
            has_slack = false;     
            % nextstage=1;
            if iter==last_hinfstb_itr
                nextstage=1;
            end
        case 3
            % STAGE=3;
            % ── Stage 2: Robust nominal (H∞ + IQC_reduced, joint obj) ────
            [LMIs, OBJ_perf,gamma] = datadriven.Objective.Hinf( ...
                G_performance, K, W_lft);
            % [LMIs_con, ~] = datadriven.Constraint.Hinf( ...
            %     G_performance, K, W_lft);
            K_fixed = datadriven.Controller.Fixed(K);
            [eta, D_val, E_val, PI_matrix, diagEta] = datadriven.Objective.IQC_reduced_PI_par_fix_eta_bisection( ...
                G_robustness, K_fixed, alpha_PI_soft, beta_PI_soft, W_lft, ...
                "etaUpperInit", eta, "etaHardCap", 100*dm_bound, "verbose", true, ...
                "multiplier",mutli_type,"epsD",1e-4,"dropNullRows",false);%,"D0",D_val,"E0",E_val);
            
            % PI_data = repmat(PI_matrix, [1, 1, length(W_lft)]);
            PI = frd(PI_matrix, W_lft, G_hat.Ts); 

            [LMIs_rob, slack] = datadriven.Objective.IQC_reduced( ...
                G_robustness, K, PI, alpha_PI_soft, beta_PI_soft, W_lft, ...
                "addSlack", false,"pinorm",true);
            LMIs      = [LMIs, LMIs_rob];
            OBJ       = OBJ_perf;
            slack_val = slack;
            has_slack = false;
         
            if obj_val > (1+loss_tol)*obj_ref || eta >= eta_max
                nextstage=1;
            end
    end

    % Tolerances
    if STAGE >= 2
        tol     = 1e-5;
        verbose = 0;
        obj_err = 1e-3;
    else
        tol     = 1e-5;
        verbose = 0;
        obj_err = 5e-3;
    end

    opts_sdp = sdpsettings('solver', solver, 'verbose', verbose, ...
        'MSK_DPAR_INTPNT_CO_TOL_NEAR_REL', 1e6);    %, ...
        % 'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',   tol, ...
        % 'mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS',   tol, ...
        % 'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP', tol, ...
        % 'mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED',  tol);
    JOB = optimize(LMIs, OBJ, opts_sdp);

    K_itr{iter}=ss(K);
    D_itr{iter}=D_val;
    E_itr{iter}=E_val;
    eta_itr{iter}=eta;


    elapsed = toc;

    if JOB.problem ~= 0
        fprintf('  !! Solver failed (code %d)\n', JOB.problem);
        break
    end

    obj_val = value(OBJ);
    if has_slack
        fprintf('  obj = %.6f | gamma = %.6f | slack = %.6f | time = %.2f s\n', ...
            obj_val,value(gamma), value(slack_val), elapsed);
    else
        fprintf('  obj = %.6f | time = %.2f s\n\n', obj_val, elapsed);
    end
    if isstable(ss(K))
        disp('Controller is OL stable')
    else
        disp('Controller is OL unstable')
    end
    err=abs(prev_obj - obj_val);
    if (err < obj_err) ||nextstage  % || (err/prev_obj < obj_err)
        nextstage=0;
        % Save intermediate controller at end of each stage
        K_stage{STAGE+1} = ss(K);
        fprintf('  --> Saved K_stage{%d}\n', STAGE+1);
        obj_ref = obj_val; 
        if STAGE < 3
            STAGE    = STAGE + 1;
            prev_obj = Inf;
            fprintf('  --> Advancing to Stage %d\n', STAGE);
        else
            fprintf('  --> Converged.\n');
            break
        end
    else
        prev_obj = obj_val;
    end
end

% Use last valid saved stage controller if solver left K in bad state
K_final = [];
try
    K_candidate = ss(K);
    [A,~,~,~]   = ssdata(K_candidate);
    if ~any(isnan(A(:))) && ~any(isinf(A(:)))
        K_final = K_candidate;
    end
catch
end

if isempty(K_final)
    warning("K contains NaN/Inf -- falling back to last valid stage controller.");
    for s = 4:-1:1
        if ~isempty(K_stage{s})
            [At,~,~,~] = ssdata(K_stage{s});
            if ~any(isnan(At(:))) && ~any(isinf(At(:)))
                K_final = K_stage{s};
                fprintf("Using K_stage{%d} as K_final.\n", s);
                break
            end
        end
    end
end

if isempty(K_final)
    error("No valid controller found in any stage.");
end

%% debug checking the value of D and E in f of the freq
for i=1:length(D_val)
    D_array(i)=D_val{i};
    E_array(i)=E_val{i};
end
figure
semilogx(W_lft/2/pi,D_array,W_lft/2/pi,imag(E_array))
legend('real(D)','imag(E)')
%% load saved controller !!!!! remove the saving part
% load(fullfile('Selected controller/','controller_sus_7nm_order10_72hz_-8db_nomperfrobperf.mat'))

%% Analysis — nominal closed-loop
S_nom       = feedback(1, G_hat * K_final);
T_nom       = feedback(G_hat * K_final, 1);
norminf_W1S = norm(interp(W1_frd, W_lft) * frd(S_nom, W_lft), inf)
norminf_T   = norm(frd(T_nom, W_lft), inf)

%% Plot sensitivity for each saved intermediate controller
stage_names = {'Stage 0 (Perf)', 'Stage 1 (Rob)', ...
               'Stage 2 (Rob+Nomhinf)', 'Stage 3 (Rob+Nom opti IQC)'};
stage_names_red={'perf','rob','nperfhinf','nperfoptiIQC'};
stagefin=[];
for i =3:4
    if ~isempty(K_stage{i})
        stagefin=[stagefin, i];
    end
end

if ~exist(fullfile('data','controller'),'dir')
    mkdir(fullfile('data','controller'))
end



expname=sprintf('%s_controller_order%d_%.0fhz_%.0fdb_%s_lasthinfstabitr%d_eta%.2f_multi_%s', ...
    idname,orderK,omega_d/2/pi,mag2db(min(abs(inv(W1_frd).ResponseData))), ...
    strcat(stage_names_red{stagefin}),last_hinfstb_itr,eta,mutli_type(1:4))

save(fullfile('data','controller',sprintf('%s.mat',expname)))
diary off
movefile(fullfile('data','diaries','diary_temp.txt'),fullfile('data','diaries',sprintf('%s.txt',expname)))


colors      = lines(4);
colors = colors([1,2,4,3],:);
figure; hold on; grid on;
for s = 1:4
    if ~isempty(K_stage{s})
        K_s = ss(K_stage{s});
        S_s = feedback(1, G_hat * K_s);
        utils.magplot(frd(S_s, W_lft), 'Color', colors(s,:));
    end
end
utils.magplot(inv(W1_frd),'Color',[0 0 1 .3]);   % performance bound

legend([stage_names(~cellfun(@isempty, K_stage)),'1/W1'], 'Location', 'best');
title('Sensitivity — intermediate controllers');
xlabel('Frequency (rad/s)');
% exportgraphics(gcf,'plot/rpall_S.svg')
%% Final plots
figure;
bodemag(inv(W1_frd), frd(S_nom, W_lft));
legend('1/W1', 'S_{final}');
title('Sensitivity — final controller');
% exportgraphics(gcf,'plot/rp_S.svg')
figure; bodemag(frd(T_nom, W_lft));  title('Complementary Sensitivity T');
% exportgraphics(gcf,'plot/rp_T.svg')
figure; bodemag(K_final);            title('Controller');
% exportgraphics(gcf,'plot/rp_K.svg')
figure; pzplot(K_final);
% exportgraphics(gcf,'plot/rp_pzplot.svg')
%% nyquist plot for stability
% figure; hold on;title('Nyquist — 1 + L(FRFinterp)')
% for i=1:length(GFRF)
%     GFRF_inter(i)=interp(GFRF(i),W_lft);
%     nyquistplot(1 + K_final*GFRF_inter(i)); 
% end
% exportgraphics(gcf,'plot/rp_nyquist.svg')

% %% compelmentary sensitivity function T
% figure;hold on;
% for i =1:length(GFRF)
%     utils.magplot(feedback(GFRF_inter(i)*K_final,1))
% end
%% plots of the various FRF sensitivities

figure; hold on
for i=1:length(GFRF_dm)
    utils.magplot(inv(1+K_stage{4}*GFRF(i)),'Color',[colors(i,:) .3],'LineStyle','--','DisplayName','dmreal='+string(GFRF_dm{i}.dmreal))
    utils.magplot(inv(1+K_stage{4}*G(i)),'Color',colors(i,:),'DisplayName','G'+string(i))
end
utils.magplot(inv(1+K_stage{4}*GLFT(2,2)),'DisplayName','nominal (G22)')
utils.magplot(inv(W1_frd))
title('robust performance (final controller)')
xlim([10 1000*pi])
% xlim([400 500])
 ylim([-40 10])

figure; hold on
for i=1:length(GFRF_dm)
    utils.magplot(inv(1+K_stage{3}*GFRF(i)),'Color',[colors(i,:) .3],'LineStyle','--','DisplayName','dmreal='+string(GFRF_dm{i}.dmreal))
    utils.magplot(inv(1+K_stage{3}*G(i)),'Color',colors(i,:),'DisplayName','G'+string(i))
end
utils.magplot(inv(1+K_stage{3}*GLFT(2,2)),'DisplayName','nominal (G22)')
utils.magplot(inv(W1_frd))
title('robust nominal performance')
xlim([10 1000*pi])
% xlim([400 500])
 ylim([-40 10])
K_itr_current=K_itr;

%%
% load(fullfile('Data/','controller_sus_7nm_order18_82hz_-34db_nperfhinfnperfoptiIQC_eta9.139390e-01.mat'),'K_itr')
iteration=17
K_chosen=K_stage{3};
%load(fullfile('Data/','controller_sus_7nm_order18_82hz_-25db_nperfhinfnperfoptiIQC_eta4.681766e-01.mat'))
K_old=K_itr_current{17};


% load(fullfile('data','controller_sus_7nm_order20_81hz_-9db_nomperf_eta0.1.mat'),'K_itr')
% K_chosen=K_itr{end};

norminf_W1S = norm(interp(W1_frd, W_lft) * frd(feedback(K_chosen*GLFT(2,2),1), W_lft), inf)
norminf_W1S_old = norm(interp(W1_frd, W_lft) * frd(feedback(K_old*GLFT(2,2),1), W_lft), inf)

figure;pzplot(K_chosen)
figure;utils.magplot(frd(K_chosen,logspace(log10(10),log10(pi/Ts),1000)),'DisplayName','K chosen')
utils.magplot(frd(K_old,logspace(log10(10),log10(pi/Ts),1000)),'Displayname','K old')
% K_chosen=K_stage{1};
if isstable(K_chosen)
    disp('Controller is OL stable')
else
    disp('Controller is OL unstable')
end
if isstable(K_old)
    disp('Controller is OL stable')
else
    disp('Controller is OL unstable')
end
colors=lines(4);
% S
figure; hold on
% load("C:\Users\pgr\Github_local\pdm-pierre-gretillat\IQC LPV identification\idm_suspended.mat")
% for i=1:length(GFRF_dm)
    % utils.magplot(inv(1+K_chosen*G(i)),'Color',colors(i,:),'DisplayName','G'+string(i))
% end
 utils.magplot(inv(1+K_chosen*GFRF(2)),'Color',[colors(2,:) .3],'LineStyle','--','DisplayName','FRF chosen dmreal='+string(GFRF_dm{2}.dmreal))
utils.magplot(inv(1+K_chosen*GLFT(2,2)),'Color',colors(2,:),'DisplayName','nominal (G22) chosen')
utils.magplot(inv(1+K_old*GFRF(2)),'Color',[colors(3,:) .3],'LineStyle','--','DisplayName','FRF old dmreal='+string(GFRF_dm{2}.dmreal))
utils.magplot(inv(1+K_old*GLFT(2,2)),'Color',colors(3,:),'DisplayName','nominal (G22) old')
utils.magplot(inv(W1_frd))
title(sprintf('K_itr %d',iteration))
xlim([10 1000*pi])
% xlim([400 500])
 ylim([-40 10])

 figure; hold on
for i=1:length(GFRF_dm)
    utils.magplot(inv(1+K_chosen*GFRF(i)),'Color',colors(i,:),'DisplayName','chosen G'+string(i))
    utils.magplot(inv(1+K_old*GFRF(i)),'Color',colors(i,:),'DisplayName','old G'+string(i),'Linestyle','--')
end
utils.magplot(inv(W1_frd))

% 
%  % T
%  figure; hold on
% for i=1:length(GFRF_dm)
%     % utils.magplot(feedback(K_chosen*GFRF(i),1),'Color',[colors(i,:) .3],'LineStyle','--','DisplayName','dmreal='+string(GFRF_dm{i}.dmreal))
%     utils.magplot(feedback(K_chosen*G(i),1),'Color',colors(i,:),'DisplayName','G'+string(i))
% end
% utils.magplot(feedback(K_chosen*GLFT(2,2),1),'DisplayName','nominal (G22)')
% title(sprintf('K_itr %d',iteration))
% xlim([10 1000*pi])
% % xlim([400 500])
%  ylim([-40 10])
% 
% 
%   figure; hold on
% for i=1:length(GFRF_dm)
%     % utils.magplot(feedback(K_chosen*GFRF(i),1),'Color',[colors(i,:) .3],'LineStyle','--','DisplayName','dmreal='+string(GFRF_dm{i}.dmreal))
%     utils.magplot(G(i),'Color',colors(i,:),'DisplayName','G'+string(i))
% end