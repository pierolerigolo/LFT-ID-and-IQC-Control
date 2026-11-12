close all;clc;clear
%% open the experimental sinesweep
expforlder=fullfile('data','Experiments/','Old Experiments/','controller_sus_7nm_order18_82hz_-34db_nperfhinfnperfoptiIQC_eta9.139390e-01_itr26/');

[S,Uclfrd]=utils.sweepS(fullfile(expforlder,'gfrf0_sswep.mat'));

load(fullfile('data','GFRF_iddata_0729_spa.mat'))

figure
hold on 
bode(Uclfrd)
bode(S)


