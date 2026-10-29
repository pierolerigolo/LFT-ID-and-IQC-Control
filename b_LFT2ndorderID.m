clear;close all;clc

if ~exist('plots', 'dir')
   mkdir('plots')
end
%% Multimodel non-parametric
expname='1028_spa';
load(fullfile('data',sprintf('GFRF_iddata_%s.mat',expname)));
bodemag(GFRF_cell{:})
legend

%% Identification of Gd(m) on a freq range of the FRF
[Gdm,Gdmnom,sys_est_Gdm,dm_bound,m0,m]=utils.Gdm_greyestid(GFRF,25,75);
exportgraphics(gcf,fullfile('plots','gdmcompare.jpg'),'Resolution',300)

%%
% s=tf('s');
% m0=0.1;
% c_est=1;
% k_est=50000;
% test=(1*s^2+2*491*s+k_est)/(m0*s^2+c_est*s+k_est);
% test2=(1*s^2)/(m0*s^2+c_est*s+k_est);
% figure
% bode(test,test2)
%%
figure;hold on
% for i=1:length(Gdm)
%     bode(Gdm(i))
% end

% figure;hold on
for i=1:length(Gdm)
    Gresid(1,i)=GFRF(i)-Gdm(i);
    bode(GFRF(i)-Gdm(i))
end
%%
Gresidmeandata=mean(Gresid.ResponseData,2);
Gresidmean=frd(Gresidmeandata,Gresid.Frequency,Gresid.Ts)
hold on
bode(Gresidmean)

%%

Parray=stack(1,Gresid(1,1),Gresid(1,2),Gresid(1,3),Gresid(1,4));
Pnom=Gresidmean;
[usys, info]=ucover(Parray,Pnom,7);
relerr = Pnom\(Pnom-Parray);
W=info.W1;
bo.FreqUnits='Hz';
figure
bodemag(relerr,'b--',W,'r',bo);




%% Identification of GLTI on the whole range
filename=sprintf('LFT_estimated_%s.mat',expname)
if exist(fullfile('data',filename),'file')
    load(fullfile('data',filename),"Gdm","m","m0","dm_bound", ...
        "Gdmnom","sys_est_Gdm","Gest","GLTI","sys_est_GLTI","GFRF_dm")
end
% n=find(cellfun(@isempty,sys_est_GLTI))
for n=1:15
    fprintf('order of GLTI: %d\n',n)
    [Gest{n},GLTI{n},sys_est_GLTI{n}] = utils.GLTI_greyestid(GFRF,n,Gdm);
    sys_est_GLTI{n}
end

save(fullfile('data',filename),"Gdm","m","m0","dm_bound", ...
   "Gdmnom","sys_est_Gdm","Gest","GLTI","sys_est_GLTI","GFRF_dm")
%% Recovering the fit percent in function of the ord of GLTI
for n=find(~cellfun(@isempty,sys_est_GLTI))
    fitpercent(n,:)=sys_est_GLTI{n}.Report.Fit.FitPercent;
end
fitpercent
mean(fitpercent,2)
% Ploting the fit percent in f of the order of GLTI
figure
plot(fitpercent,'o-')
hold on
plot(mean(fitpercent,2),'x-','Linewidth',2)
xlim([1 length(sys_est_GLTI)])
legend('OP 1','OP 2','OP 3','OP 4','Mean')
xlabel('n')
ylabel('fit%')

%% Comparing with another setting
% load(fullfile('data','full_est_sys_data_stbenf_fullfrf.mat'))
% for n=find(~cellfun(@isempty,sys_est_GLTI))
%     fitpercent(n,:)=sys_est_GLTI{n}.Report.Fit.FitPercent;
% end
% fitpercent
% mean(fitpercent,2)
%%
% figure
% plot(fitpercent,'o-')
% hold on
% plot(mean(fitpercent,2),'x-','Linewidth',2)
% xlim([1 15])
% legend('OP 1','OP 2','OP 3','OP 4','Mean')
% 
%% Comparaison and validation of the chosen model order
% compare with NRMSE fit %
figure
compare(Gest{6},GFRF)

%% compare the db difference in function of the frequency
for i=1:4
    dberr(i)=frd(abs(freqresp(Gest{10}(1),GFRF.Frequency))./abs(GFRF(1).ResponseData),GFRF.Frequency);
end
bo=bodeoptions;
bo.Grid='on';

figure;bodemag(dberr,bo)
exportgraphics(gcf,fullfile('plots','dberr.jpg'),'Resolution',300)
%% poles zeros cancellation
figure
h = iopzplot(sys_est_GLTI{7});
showConfidence(h, 2);
%%
% Residual analysis (not significant due to the low fit of the model and
% the fact that it is an infinite order system
G_est_current=Gest{7};

for i=1:4
    figure
    resid(data{i},G_est_current(i))
    title(sprintf('OP %d',i))
    exportgraphics(gcf,fullfile('plots',sprintf('residop%d.jpg',i)),'Resolution',300)
end

% Temporal validation of the select model 
for i=1:4
    figure
    compare(G_est_current(i),data{i})
    title(sprintf('OP %d',i))
    exportgraphics(gcf,fullfile('plots',sprintf('tempcompop%d.jpg',i)),'Resolution',300)
end

%% Magnitude plot of only the LTI part fitted
bodemag(GLTI{8})


%% Fit percent only on the magnitude for comparision with the old method
for i=1:4
    Gest_temp=Gest{7};
    disp(utils.NRMSE(Gest_temp(i),GFRF(i),'magnitudeonly',true))
end

%% Various order of the system compared to the FRF in one graph (only 1 OP)
load(fullfile('data','full_est_sys_data.mat'),'Gest',"GLTI","sys_est_GLTI")
figure
for i=7:12
    Gest{i}.Name=sprintf('n=%d',i);
end
compare(GFRF(1),Gest{7}(1),Gest{8}(1),Gest{9}(1),Gest{10}(1),Gest{11}(1),Gest{12}(1))

%% Coherance plot 
G_est_chosen=Gest{7};
[ymod,fit,ic]=compare(data{1},G_est_chosen(1));
figure
ymeas=data{1}.OutputData;
ymodel=ymod.OutputData;
mscohere(ymeas,ymodel)
figure
plot(ymeas)
hold on
plot(ymodel)
%-> oke jcrois ca sert à rien la coherance

%% Abs diff in f of freq of bt the model and the measurements
figure
title('dB absolute error between model vs measurements')
hold on
for i=1:length(GFRF)
    [absdbfit,nonnanid]=utils.NRMSE(G_est_chosen(i),GFRF(i), ...
        'inFofFreq',true,'absDB',true);
    plot(GFRF.Frequency(nonnanid)/2/pi,absdbfit,'DisplayName',sprintf('OP %d',i))
end
xscale('log')
xline(25,'Label','Freq. min.','HandleVisibility','off')
xline(75,'Label','Freq. max.','HandleVisibility','off')
xlabel('Frequency (Hz)')
ylabel('Error (dB)')
legend('Location','best')
xlim([10 500])
%% Prepare the LFT for the control synthesis

% Frequency grid
% Main coarse grid
omega_global = logspace(log10(10), log10(pi/GFRF.Ts), 150);
% Tiny dense patch: ±30 Hz around the peak (25–100 Hz)
omega_peak_dense = logspace(log10(25*2*pi), log10(100*2*pi), 150);
% omega_peak_dense=[];
w_lft = unique([omega_global(:); omega_peak_dense(:)]);

% Selected order of GLTI
ord=7;

GLTI_selected=frd(GLTI{ord},w_lft);

GLFT= [  -Gdmnom, Gdmnom*GLTI_selected;...
         -Gdmnom, Gdmnom*GLTI_selected];

G=Gest{ord};

save(fullfile('data',sprintf('LFT_GFRF_param_n%d_%s.mat',ord,expname)), ...
    "GLFT","dm_bound","m0","m","GFRF","G","GFRF_dm")
%% Table for the obtained parameter and their std 
function T = param_table(sys_est, names)
% Parameter estimates and standard deviations of a greyest model.
% names: string/cellstr, one entry per estimated parameter (in ParVector order)
    P   = sys_est.Report.Parameters;
    val = P.ParVector(:);
    C   = P.FreeParCovariance;

    if isempty(C) || all(isnan(C(:)))
        sd = nan(size(val));
        warning('param_table:noCov','No covariance available for this model.')
    else
        sd = sqrt(diag(C));
    end

    names = string(names(:));
    assert(numel(names) == numel(val), ...
        'Expected %d names, got %d.', numel(val), numel(names));

    T = table(names, val, sd, 100*sd./abs(val), ...
        'VariableNames',{'Parameter','Value','StdDev','RelStd_pct'});
end

% G_Delta(m): parameters are m_1..m_n, c, k
n     = numel(GFRF);
names = [compose("m_%d",1:n), "c", "k"];
T_Gdm = param_table(sys_est_Gdm, names)

% G_LTI (part of G_est): b_0..b_ord, a_1..a_ord
ord    = 7;                                   % the order you keep
names  = [compose("b_%d",0:ord), compose("a_%d",1:ord)];
T_GLTI = param_table(sys_est_GLTI{ord}, names)

% Formatted "value ± std" strings for the slides
fmt = @(T) T.Parameter + " = " + compose("%.3g",T.Value) + " ± " + compose("%.2g",T.StdDev);
disp(fmt(T_Gdm))
disp(fmt(T_GLTI))

% Export (paste into PowerPoint as a table)
writetable(T_Gdm,  fullfile('plots','params_Gdm.xlsx'))
writetable(T_GLTI, fullfile('plots','params_GLTI.xlsx'))

%% using multiplicative uncertain system to replace GLTI optional needs to be reviewed
% combine uncertain system G=Gdm*Gw2
for i = 1:4
    P(i)=GFRF(i)/Gdm(i);
end
Parray=stack(1,P(1),P(2),P(3),P(4));
Pnom=GLTI{7};
[usys, info]=ucover(Parray,GLTI{7},7);
relerr = Pnom\(Pnom-Parray);
W=info.W1;
bo.FreqUnits='Hz';
figure
bodemag(relerr,'b--',W,'r',bo);
%%
% only G=Gw2
Parray=stack(1,GFRF(1),GFRF(2),GFRF(3),GFRF(4));
Pnom=Gdmnom%*GLTI{7};
[usys, info]=ucover(Parray,Pnom,7);
relerr = Pnom\(Pnom-Parray);
W=info.W1;
bo.FreqUnits='Hz';
figure
bodemag(relerr,'b--',W,'r',bo);
%%
% high pass filter
wbd=2*pi*100;
hpf=tf([1 0 0],[1 2*wbd wbd^2]);
hpf=c2d(hpf,0.001);
figure;bode(hpf)
figure
bodemag(relerr,'b--',hpf*W,'r',bo);
