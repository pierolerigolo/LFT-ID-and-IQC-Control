clear;close all;clc

% Folder with the sinesweep IO data pairs
idfolder=fullfile('Data','MIVIDA 260729 sswp');

%% nouvelle id avec les rondelles sous la masse de 3.09 kg
GFRF0_2907=utils.expdata2frf_stft(fullfile(idfolder,'260729_7nm_p10_rdm_0kg_x_sswp.mat'), 'Fs',1000);
GFRF1_2907=utils.expdata2frf_stft(fullfile(idfolder,'260729_7nm_p10_rdm_2.12kg_x_sswp.mat'), 'Fs',1000);
%% id avec les 2 masse de 3.8 kg
GFRF2g_2907=utils.expdata2frf_stft(fullfile(idfolder,'260729_rdm_3.8kg_x_sswp.mat'), 'Fs',1000);
GFRF4_2907=utils.expdata2frf_stft(fullfile(idfolder,'260729_rdm_7.6kg_x_sswp.mat'),'Fs',1000);

%% temporal id data
[u{1}, y{1}]=utils.openIOtemp(fullfile(idfolder,'260729_7nm_p10_rdm_0kg_x_sswp.mat'));
[u{2}, y{2}]=utils.openIOtemp(fullfile(idfolder,'260729_7nm_p10_rdm_2.12kg_x_sswp.mat'));
[u{3}, y{3}]=utils.openIOtemp(fullfile(idfolder,'260729_rdm_3.8kg_x_sswp.mat'));
[u{4}, y{4}]=utils.openIOtemp(fullfile(idfolder,'260729_rdm_7.6kg_x_sswp.mat'));

%%
data=cell(1,4);
for i=1:4
    data{i}=iddata(y{i}',u{i}',0.001);
end
%% Create the various variable used
GFRF_dm{1}.FRF=GFRF0_2907;
GFRF_dm{1}.dmreal=0;
GFRF_dm{2}.FRF=GFRF1_2907;
GFRF_dm{2}.dmreal=2.12;
GFRF_dm{3}.FRF=GFRF2g_2907;
GFRF_dm{3}.dmreal=3.8;
GFRF_dm{4}.FRF=GFRF4_2907;
GFRF_dm{4}.dmreal=7.6;

% Multimodel as a MIMO of 1Xn
freq=GFRF_dm{1}.FRF.Frequency;
for i=1:length(GFRF_dm)
    GFRF(1,i)=interp(GFRF_dm{i}.FRF,freq);
end
save(fullfile('Data','GFRF_iddata_2907_stft.mat'),'GFRF','GFRF_dm','data')