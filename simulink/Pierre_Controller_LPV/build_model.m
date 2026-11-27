%% build_model.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('C:\Projects\OSIP-VIB\Simulink\PRBS')
addpath('Controller_LTI\LFT_IQC\')
addpath('Controller_LTI\Screw_Torque_Uncertainty')

clearvars
close all
clc

%% System parameters
Switch_controller = 0;

% T_end   = 20;
Fs      = 1e3;     % Hz
Ts      = 1/Fs;    % s
V_sat   = 2/10; % V (factor of 10 for dSpace conversion)
 
V_sat_TMS = 5; % V (factor of 10 for dSpace conversion)
% Sensitivities Wilcoxon accelerometer
% SN #4546
S(1)            = 10.5;     % V/g
% SN #4786
S(2)            = 10.4;     % V/g
% SN #4787
S(3)            = 10.1;     % V/g
% Kistler 8786A5 204939
S_Kistler_39    = 1.002;    % V/g
% Kistler 8786A5 204940
S_Kistler_40    = 1.03;     % V/g

% Sensitivity Impulse Hammer 086E80
S_Hammer        = 0.02395;  % mV/N

% Sensitivity Force Sensor PCB208C01
S_Sensor    = 0.1134;

% Sensitivities Micromega accelerometer
% R(x) where x corresponds to MICROMEGA_ID
R(1)            = 0.898;     % V/g
R(2)            = 0.902;     % V/g
R(3)            = 0.901;     % V/g
R(4)            = 0.897;     % V/g
R(5)            = 0.901;     % V/g
R(6)            = 0.899;     % V/g

% Choose IMU type
TYPE = 'MF';

% Sensitivities ASC accelerometer
switch TYPE
    case 'MF'
        Sen_AX            = 0.5424;     % V/g
        Sen_AY            = 0.53845;    % V/g
        Sen_AZ            = 0.54161;    % V/g
    case 'LN'
        Sen_AX            = 2.0176;     % V/g
        Sen_AY            = 2.011;      % V/g
        Sen_AZ            = 2.0066;     % V/g
    otherwise
        error('Selected accelerometer type is not valid!')
end

% Sensitivity ASC gyroscope
switch TYPE
    case 'MF'
        Sen_GX          = 0.18811;      % V/rad/s
        Sen_GY          = 0.18991;      % mV/rad/s
        Sen_GZ          = 0.18764;      % mV/rad/s
    case 'LN'
        Sen_GX          = 0.75481;      % V/rad/s
        Sen_GY          = 0.75728;      % mV/rad/s
        Sen_GZ          = 0.74931;      % mV/rad/s
    otherwise
        error('Selected accelerometer type is not valid!')
end

% Generate PRBS excitation signal
amp     = 3;                        % V
N = 13;                             % Register length
Rep = 4;                            % Repetitions
P = 2^N-1;                          % Period
PRBS = genPRBS_freqDiv(N,Rep,Ts,1);   % Generate PRBS signal
%     signal = amp*[PRBS.u' PRBS.u' PRBS.u'];
time = 0:Ts:(PRBS.t(end)+20);
signal = amp*[zeros(1,10*Fs) PRBS.u' 10*zeros(1,10*Fs)];  % Add zeros in the beginning to account for DAQ time lag in ControlDesk
%     time = linspace(0, PRBS.t(end)+10, length(PRBS.t)+10*Fs); % Retrieve PRBS signal time and add time for zero signal in the beginning and end
ext_PRBS = timeseries(signal, time);
% ext_PRBS.signals.values = signal;
% ext_PRBS.time = [];
T_end = 20+PRBS.t(end);
% test signal
% [corr, lags] = xcorr(PRBS.u,PRBS.u);

% figure
% plot(time, signal)
% title('PRBS excitation signal')
%
%     figure
%     plot(lags, corr)
%     title('PRBS autocorrelation')

% FFT calculation
% [freq, Fz_fft] = FFT_Abs(signal, Fs);
% % FFT_db(Fz, Fs);
% figure
% plot(freq, Fz_fft)
% % title(['FFT of force, excitation frequency = ' filename(7:8) ' Hz'])
% xlim([1 500])
% xlabel('Frequency (Hz)')
% ylabel('Force (N)')


% Generate sinusoidal frequency sweep
f1 = 10;
f2 = 500;
time = 0:Ts:f2/10;
f= (f1+(f2-f1)/(2*time(end))*time);
signal = sin(2*pi*f.*time);
signal = amp*signal(:);

signal = [zeros(1,5*Fs) signal' zeros(1,5*Fs)];  % Add zeros in the beginning to account for DAQ time lag in ControlDesk
time = 0:Ts:f2/10+10;
%     
%     figure
%     plot(time, signal)
ext_SineSweep.signals.values = signal';
ext_SineSweep.time = [];
%     excitation = timeseries(signal,time);
% end

% Compensation filter
zeta = 0.1623;
omega_n = 24.44*2*pi;
zeta_c = 1.6;
omega_c = 15*2*pi;  % Alternatively 10, 15 Hz

Numerator = [1 2*zeta*omega_n omega_n^2];
Denominator = [1 2*zeta_c*omega_c omega_c^2];

C = tf(Numerator,Denominator);

Cz = c2d(C,Ts);
% Cz = tf2ss(Cz.Numerator{1}, Cz.Denominator{1});

options = bodeoptions;
options.FreqUnits = 'Hz';
% options.MagUnits = 'abs';
options.XLim = [1 500];
options.PhaseMatching = 'on';
options.PhaseMatchingFreq = 1;
% 
% figure
% bode(Cz, Gz, Cz*Gz, options)
% legend show
% title(['Compensation filter: zeta_c = ' num2str(zeta_critical) ', omega_c = ' num2str(omega_controller/2/pi) 'Hz'])

% Controller (feedback)
% addpath('Controller_LTI')
% load('K_mixsyn_CompFilt.mat')
% load('K_mixsyn_DistModel.mat')
% load('Kdd_G_10.mat')
% load('Kdd_FRF_3.mat')
% load('220816/Kdd_G_11.mat')
% load('Kdd_FRF_7.mat')
% load('Kdd_Gtesttt_mixsyn.mat')
% load('Kdd_G_order8_freq.mat')
% load('Kdd_G_freq.mat')
% load('Kdd_FRF_order8_freq.mat')
% load('DD_FRF_freq30_100_order5_unsure.mat');
% index_freq = 6;
% index_iter = 6;
% load('MM_Plate_Torque.mat')
% load('W2_Plate_Torque.mat')
% load('IQC_Payload_Torque_corrected.mat')
% load('230201_30dB.mat')
% load('Plate_Torque_2_comparison_6_7_8Nm.mat')
% load('Plate_Torque_2_comparison_6_7_8Nm_corrected.mat')
% load('Payload_Torque_comparison.mat')
% load('Plate_Torque_2_comparison_v2.mat')
% load('Plate_Torque_2_comparison_Vaibhav_final.mat')
% load('Plate_Torque_2_comparison_25dB.mat')
% load('K_IQC_LFR.mat')
% load('Controller_PlateTorque_30dB_v2.mat')
% load('K_IQC_LFR_RobustPerformance.mat')

% load('data.mat')
% K = K_hindi_nominal;
% K = K_hindi_robust;
% K = K_new_nominal;
% K = K_new_robust;

% K = Kdd{6};
% K = K_MM{3};
% K = K_W2{11};
% K = K_IQC{8};
% K_working = Kdd_freq{index_freq,index_iter};
% K = Kdd.*[1 0 0; 0 0 0; 0 0 0];
% K = K_IQC_V;
% K = K_IQC_V_disk;
% K = K_IQC_V;
% K = Kc;
% K = K_IQC_LFR;
% K = K_proposed;

% figure
% bodemag(K_working,options)
% title('Working controller')
% K = Kdd(1,1);
% K = [K 0 0; 0 0 0; 0 0 0];
% K = ss(1*K);

% load(fullfile('controller_data/','controller_sus_7nm_order10_72hz_-8db_nomperfrobperf_works_stbg0.5.mat'))
% K_old=K_itr{18};
% w1_old=W1_frd;
% load(fullfile('controller_data/','controller_sus_7nm_order10_72hz_-1db_nomperfrobperf.mat'))
% K_robperf=K_itr{20};
% w1_robperf=W1_frd;
% load(fullfile('controller_data/','controller_sus_7nm_order10_72hz_-2db_nomperf.mat'))
% K_nomperf=K_itr{15};
% w1_nomperf=W1_frd;
% 
% K=K_robperf;
% 
% W1_frd_current=w1_robperf;


% load(fullfile('controller_data/','controller_sus_7nm_order18_81hz_-26db_nomperfrobperf_eta2.515400e-01.mat'))
% K=K_itr{10}; %marche mais pas a vide (FRFmean /1.5 dans lidentificaiton)
% 
% load(fullfile('controller_data/','controller_ellips_sus_7nm_order21_81hz_-26db_nomperf.mat'))
% K=K_itr{end};% %marche mais perf nul

% load(fullfile('controller_data/','controller_sus_7nm_order16_81hz_-26db_nomperfrobperf_eta3.018480e-01.mat'))
% K=K_itr{5}; % itr 8 marche à vide mais pas avec 5kg itr5 marche pour les 2 extreme


% eta*1 but no tweaking on the frf mean
% 
% load(fullfile('controller_data/','controller_sus_7nm_order16_81hz_-26db_nomperfrobperf_eta2.515400e-01.mat'))
% K=K_itr{24}; marche et a besoin des rondelle ppur le disque

% load(fullfile('controller_data/','controller_2907_7nm_order18_82hz_-30db_nomperfrobperf_eta3.446890e-01.mat'))
% K=K_itr{11}; %nom perf works and has the best inf norm

% iqc opti qui marche en partant de lancien nom perf qui marchait 
% load(fullfile('controller_data/','controller_sus_7nm_order18_82hz_-30db_nperfhinfnperfoptiIQC_eta5.106693e-01.mat'))
% K=K_itr{20}; %20 focntio yay eta = 0.495194

load(fullfile('controller_data/','controller_sus_7nm_order18_82hz_-34db_nperfhinfnperfoptiIQC_eta9.139390e-01.mat'))
% K=K_stage{3}; % eta initial de dm_max
K=K_itr{26}; %work till end at 3.8kg work till 26 with 0kg

% load(fullfile('controller_data/','controller_sus_7nm_order18_82hz_-40db_nperfhinfnperfoptiIQC_eta6.511209e-01.mat'))
% K=K_itr{end}; %works with 3.8kg a gauche


% load(fullfile('controller_data/','controller_sus_7nm_order18_82hz_-25db_nperfhinfnperfoptiIQC_eta4.681766e-01.mat'))
% K=K_itr{end} work till end with 0kg but less perf thatn the -34 at itr 26

% load(fullfile("controller_data/",'controller_sus_7nm_order18_82hz_-30db_nperfhinfnperfoptiIQC_eta3.118336e-01.mat'))
% K=K_itr{11}; % de lamerde marche pas l'eta est jamais assez grand

% load(fullfile('controller_data/','controller_sus_7nm_order18_82hz_-32db_nperfhinfnperfoptiIQC_eta4.441903e-01.mat'))
% % K=K_stage{3};
% K=K_itr{22}; % denouveau full pas stable mme avec eta plus grand

% load(fullfile('controller_data/','controller_sus_7nm_order20_82hz_-32db_nperfhinf_eta7.568131e-01.mat'))
% K=K_stage{3}; %eta 0.2, stable avec masse de 2kg instable a vide donc ce quon attend 
% K=K_itr{18}; % denouvveau stable juste a cet iter a vide mais le eta est
% vraiment plus grand et genre a17 ca marche pas a 19 ca marche plus

% load(fullfile('controller_data/','controller_sus_7nm_order18_61hz_-34db_nperfhinfnperfoptiIQC_eta4.931591e-01.mat'))
% K=K_itr{end}; %59.94 hz jcrois la 2e decimal pas sur

%49.76hz
% load(fullfile('controller_data/','controller_sus_7nm_order18_50hz_-33db_nperfhinfnperfoptiIQC_eta5.851874e-01.mat'))
% K=K_itr{12}; %marche po

% load(fullfile('controller_data/','controller_sus_7nm_order18_50hz_-29db_nperfhinfnperfoptiIQC_eta5.421448e-01.mat'))
% K=K_itr{10}; %marche pour itr=10 % ORX WAS HERE ? 

% load(fullfile('controller_data','controller_ellips_sus_7nm_order20_82hz_-30db_nomperf.mat'))
% K=K_itr{5};

% load(fullfile('controller_data','controller_ellips_sus_7nm_order16_82hz_-30db_nomperf.mat'))
% K=K_itr{7};

% load(fullfile('LFT-ID-and-IQC-Control/','Data/','controller/','LFT_GFRF_param_n7_2907_controller_order18_82hz_-34db_nperfhinfnperfoptiIQC_lasthinfstabitrInf_eta1.02_multi_freq.mat'),'K_itr','W1_frd')
% K=K_itr{3};

fdname=fullfile('LFT-ID-and-IQC-Control/','Data/','controller/');
% load(fullfile(fdname,'LFT_GFRF_param_n7_1028_controller_order16_82hz_-34db_nperfhinfnperfoptiIQC_lasthinfstabitrInf_eta0.77_multi_freq.mat'),'K_itr','K_stage')
% K=K_itr{5};

% load(fullfile(fdname,'LFT_GFRF_param_n7_1028_controller_order18_82hz_-34db_nperfhinfnperfoptiIQC_lasthinfstabitr10_dminit0.41_etaf0.62_multi_freq.mat'),'K_itr','K_stage')
% K=K_itr{3};

% opensave=load(fullfile(fdname,'LFT_GFRF_GdmGLTIW2opt_1028_etfe_controller_order18_80hz_-34db_nperfoptiIQC_notconv_etai0.36_etaf0.82_multi_freq.mat'));
% K=opensave.K_stage{2};
% % K=opensave.K_stage{2};
% opensave=load(fullfile(fdname,'ELLI_1028_etfe_controller_order18_80hz_-34db_nperfhinf.mat'));
% K=opensave.K_itr{21};
% % K=opensave.K_stage{2};

%
fn='LFT_GFRF_GdmGLTIW2opt_1110_etfe_80hz_fmin10_hyb_controller_order18_80hz_-34db_nperfoptiIQC_notconv_etai0.35_etaf0.36_multi_freq'
% fn='DISK_1028_etfe_controller_order18_80hz_-34db_nperfhinf';
opensave=load(fullfile(fdname,[fn,'.mat']));
% K=opensave.K_stage{2};
K=opensave.K_final;


% Nominal second order model
modelsave=load('../../data/LFT_GFRF_GdmGLTIW2opt_1110gc_etfe_80hz_fmin10_hyb.mat');
[Gnomnum,Gnomden]=tfdata(modelsave.Gdmnom);


% opensave=load(fullfile('LFT-ID-and-IQC-Control/Data/finalrepportcont/', ...
%    'controller_ellips_sus_7nm_order20_82hz_-30db_nomperf.mat' ));
% K=opensave.K_itr{end};



% figure; hold on;
% for i =25:5:50
%     bode(K_itr{i})
% end
W1_frd_current=opensave.W1_frd;
figure
% bode(K,K_working,options)
pzplot(K)
figure
bode(K,options)
legend
% legend('Robust K', 'Classic K')

if isstable(K)
    disp('Controller is OL stable')
else
    disp('Controller is OL unstable')
end
K(2,2)=c2d(tf(0),Ts);
K(3,3)=c2d(tf(0),Ts);
% figure
% bodemag(K, options)

% checking s and t with the controller and the last measured frf
load("C:\Users\pgr\Github_local\pdm-pierre-gretillat\controller syn\Data\sus_data_2907_lft_full.mat")
W_lft=GLFT.frequency;
norminf_W1S = norm(interp(W1_frd_current, W_lft) * frd(feedback(K(1,1)*GLFT(2,2),1), W_lft), inf)
% %%
% % figure;hold on;bode(inv(W1_frd_current),options)
% for i =1:length(GFRF_dm)
%     GFRF = GFRF_dm{i}.FRF;
%     S = feedback(1, GFRF * K(1,1));
%     bode(S,options)
% end
% figure;hold on;
% for i =1:length(GFRF_dm)
%     GFRF = GFRF_dm{i}.FRF;
%     T = feedback(GFRF * K(1,1), 1);
%     bode(T,options)
% end    


%%
% Proportional gains for x,y,z-axes
Px = 1;
Py = 1;
Pz = 1;

% Butterworth lowpass filter desing for SDFT signal preconditioning
n = 6;
Wn = 200*Ts;
[A,B,C,D] = butter(n,Wn,'low');

BUTTER = ss(A,B,C,D,Ts);

% options = bodeoptions;
% options.FreqUnits = 'Hz';
% options.Xlim = [1 1000];
% figure
% bodemag(BUTTER,options)

% Parameter variables for commands in dSpace
TMS_amp=2;
% freq = 30:10:100;
% TMS_freq1=freq(index_freq);
TMS_freq1=80;
TMS_freq2=0;
TMS_freq3=0;
Sweep_amp = 0.5;

currentswitch=2;

Switch_excit=4;


%% Build 
rti_build2('Pierre_Controller_LTI')

