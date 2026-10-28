clear;close all;clc

% Folder with the sinesweep IO data pairs
idfolder_old=fullfile('Data','MIVIDA 260729 sswp');
idfolder=fullfile('Data','MIVIDA 261028/');

%% nouvelle id avec les rondelles sous la masse de 3.09 kg
GFRF0_2907_sswp_stft=utils.expdata2frf_stft(fullfile(idfolder_old,'260729_7nm_p10_rdm_0kg_x_sswp.mat'), 'Fs',1000);
GFRF1_2907=utils.expdata2frf_stft(fullfile(idfolder_old,'260729_7nm_p10_rdm_2.12kg_x_sswp.mat'), 'Fs',1000);
%% id avec les 2 masse de 3.8 kg
GFRF2g_2907=utils.expdata2frf_stft(fullfile(idfolder_old,'260729_rdm_3.8kg_x_sswp.mat'), 'Fs',1000);
GFRF4_2907=utils.expdata2frf_stft(fullfile(idfolder_old,'260729_rdm_7.6kg_x_sswp.mat'),'Fs',1000);

%% PRBS vs Sinesweep
% GFRF0_sinesweepDFT=utils.expdata2frf(fullfile(idfolder,'261028_gfrf0kg_xsswp.mat'),65000);
GFRF0_sinesweepDFT=[];
GFRF0_sinesweepSTFT_a=utils.expdata2frf_stft(fullfile(idfolder,"261028_gfrf0kg_xsswp.mat"),'itmax',60000,'Fs',1000);
GFRF0_sinesweepSTFT_b=utils.expdata2frf_stft(fullfile(idfolder,"261028_gfrf0kg_xsswp.mat"),'itmin',60000,'itmax',120000,'Fs',1000);
GFRF1_sinesweepSTFT=utils.expdata2frf_stft(fullfile(idfolder,"261028_gfrf2.1kg_xsswp.mat"),'Fs',1000);
GFRF2_sinesweepSTFT=utils.expdata2frf_stft(fullfile(idfolder,"261028_gfrf5.9kg_xsswp.mat"),'Fs',1000);
GFRF3_sinesweepSTFT=utils.expdata2frf_stft(fullfile(idfolder,"261028_gfrf7.6kg_xsswp.mat"),'itmax',65000,'Fs',1000);


% prbs frf average on 5 periods
[u,y]=utils.openIOtemp(fullfile(idfolder,"261028_gfrf0kg_xprbs.mat"),8000,60000);

[Ru,Lag]=xcorr(u);
figure;plot(Lag,Ru)
T_period=8191; %=2^13-1 13bit shift register
N=length(u)
m = floor(N / T_period);                 
fprintf('Number of periods: m = %d\n', m);
Ts=0.001;
u_periods = reshape(u(1:T_period*m), T_period, m);
y_periods = reshape(y(1:T_period*m), T_period, m);

% Average spectra first, then take ratio
U_avg = mean(fft(u_periods), 2);   
Y_avg = mean(fft(y_periods), 2);   
G_avg = Y_avg ./ U_avg;

freq_hz = (0:T_period-1) / (T_period * Ts);
idx = 1:floor(T_period/2);

GFRF0_prbsfourier = frd(G_avg(idx), freq_hz(idx),Ts,'FrequencyUnit','Hz');

fig=figure;
hold on
bode_opts = bodeoptions; bode_opts.FreqUnits = 'Hz';
% bodemag(GFRF0_prbsfourier,GFRF0_sinesweepSTFT_a,GFRF1_sinesweepSTFT,GFRF2_sinesweepSTFT,GFRF3_sinesweepSTFT,GFRF0_2907_sswp_stft,bode_opts); grid on;
bodemag(GFRF0_prbsfourier,GFRF0_sinesweepSTFT_a,GFRF0_2907_sswp_stft,bode_opts); grid on;
legend
xlim([10 500])
title(sprintf('Fourier FRF estimate (prbs m=%d periods averaged)', m));
% exportgraphics(gcf,fullfile('plots','prbsvssineseeep.jpg'),'Resolution','600')
%% temporal id data
% [u{1}, y{1}]=utils.openIOtemp(fullfile(idfolder_old,'260729_7nm_p10_rdm_0kg_x_sswp.mat'));
% [u{2}, y{2}]=utils.openIOtemp(fullfile(idfolder_old,'260729_7nm_p10_rdm_2.12kg_x_sswp.mat'));
% [u{3}, y{3}]=utils.openIOtemp(fullfile(idfolder_old,'260729_rdm_3.8kg_x_sswp.mat'));
% [u{4}, y{4}]=utils.openIOtemp(fullfile(idfolder_old,'260729_rdm_7.6kg_x_sswp.mat'));

[u{1}, y{1}]=utils.openIOtemp(fullfile(idfolder,'261028_gfrf0kg_xsswp.mat'),1,60000);
[u{2}, y{2}]=utils.openIOtemp(fullfile(idfolder,'261028_gfrf2.1kg_xsswp.mat'));
[u{3}, y{3}]=utils.openIOtemp(fullfile(idfolder,'261028_gfrf5.9kg_xsswp.mat'));
[u{4}, y{4}]=utils.openIOtemp(fullfile(idfolder,'261028_gfrf7.6kg_xsswp.mat'),1,65000);

%%
data=cell(1,4);
for i=1:4
    data{i}=iddata(y{i}',u{i}',0.001);
end
%% Create the various variable used
% GFRF_dm{1}.FRF=GFRF0_2907;
% GFRF_dm{1}.dmreal=0;
% GFRF_dm{2}.FRF=GFRF1_2907;
% GFRF_dm{2}.dmreal=2.12;
% GFRF_dm{3}.FRF=GFRF2g_2907;
% GFRF_dm{3}.dmreal=3.8;
% GFRF_dm{4}.FRF=GFRF4_2907;
% GFRF_dm{4}.dmreal=7.6;

GFRF_dm{1}.FRF=GFRF0_sinesweepSTFT_a;
GFRF_dm{1}.dmreal=0;
GFRF_dm{2}.FRF=GFRF1_sinesweepSTFT;
GFRF_dm{2}.dmreal=2.12;
GFRF_dm{3}.FRF=GFRF2_sinesweepSTFT;
GFRF_dm{3}.dmreal=5.9;
GFRF_dm{4}.FRF=GFRF3_sinesweepSTFT;
GFRF_dm{4}.dmreal=7.6;




% Multimodel as a MIMO of 1Xn
freq=GFRF_dm{1}.FRF.Frequency;
for i=1:length(GFRF_dm)
    GFRF(1,i)=interp(GFRF_dm{i}.FRF,freq);
end
save(fullfile('Data','GFRF_iddata_1028_stft.mat'),'GFRF','GFRF_dm','data')