clear;close all;clc;

[S, Uclfrd, Yclfrd, Yolfrd] = utils.sweepS(fullfile('data','Experiments/','gfrf0_K_final_''controller_ellips_sus_7nm_order20_82hz_-30db_nomperf.mat'));
[S_old, Uclfrd, Yclfrd, Yolfrd] = utils.sweepS(fullfile('data','Experiments/','Old Experiments/','gfrf0_ellips_80hz_sswep_order20.mat'));
[S_stbg5_6, Uclfrd, Yclfrd, Yolfrd] = utils.sweepS(fullfile('data','Experiments/','gfrf0_K_final_''controller_ellips_sus_7nm_order20_82hz_-30db_nomperf_stbgain5.6.mat'));

contsave=load(fullfile('data','finalrepportcont/','controller_ellips_sus_7nm_order20_82hz_-30db_nomperf.mat'));

K_save=contsave.K_itr{end};

load(fullfile('data','GFRF_iddata_1028_etfe.mat'));

bode(S,S_old,S_stbg5_6,inv(1+GFRF(1)*K_save),inv(1+GFRF(1)*K_save/5.6))
legend