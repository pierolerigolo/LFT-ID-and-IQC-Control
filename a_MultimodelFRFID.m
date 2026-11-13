%% a_MultimodelFRFID  -  Step 1: non-parametric multimodel FRF identification
%
%   One sine-sweep experiment per operating point (added mass) -> one FRF
%   per operating point, estimated with spa.
%
%   Output (when save_output = true): Data/GFRF_iddata_<date>_etfe.mat (outfile)
%     GFRF      1 x nOP frd (channel i = operating point i)
%     GFRF_cell 1 x nOP cell of idfrd (same FRFs)
%     d         1 x nOP cell of iddata (detrended time data)
%     realmass  added mass of each operating point (kg)
%
clear; close all; clc

%% Config
% one sine-sweep file per operating point (0 kg re-identified on 2026-11-10)
files    = {fullfile('Data', 'MIVIDA 261110', '261110_gfrf0kg_xsswp.mat'), ...
            fullfile('Data', 'MIVIDA 261028', '261028_gfrf2.1kg_xsswp.mat'), ...
            fullfile('Data', 'MIVIDA 261028', '261028_gfrf5.9kg_xsswp.mat'), ...
            fullfile('Data', 'MIVIDA 261028', '261028_gfrf7.6kg_xsswp.mat')};
realmass = [0 2.1 5.9 7.6];          % added mass of each file (kg)
itmax    = {[], [], [], 65000};      % last sample kept per file ([] = whole record)

Ts  = 0.001;
% FRF frequency vector (rad/s): log spacing f_lo..f_split, linear spacing above
% (fine enough to resolve the high-frequency modes)
f_lo    = 10;          % lowest frequency (Hz)
f_split = 20;          % log below, linear above (Hz)
n_log   = 50;          % points f_lo..f_split
n_lin   = 1000;        % points f_split..Nyquist (~0.5 Hz spacing)
n_etfe  = 2000;        % etfe frequency points (0..Nyquist, ~0.25 Hz)
w_lin = linspace(2*pi*f_split, pi/Ts, n_lin + 1);
w   = [logspace(log10(2*pi*f_lo), log10(2*pi*f_split), n_log), w_lin(2:end)];
lag = 1000;      % spa lag window; [] = record length (no smoothing), smaller = smoother

save_output = true;   % true: (over)write the output file below
outfile     = fullfile('Data', 'GFRF_iddata_1110_etfe.mat');

%% FRF of each operating point
nOP       = numel(files);
GFRF_cell = cell(1, nOP);
d         = cell(1, nOP);
for i = 1:nOP
    [u, y] = utils.openIOtemp(files{i}, 1, itmax{i});

    d{i} = iddata(y(:), u(:), Ts, 'InputName', 'u'+string(i), 'OutputName', 'y');
    d{i} = detrend(d{i});                          % remove the operating-point offset

    if isempty(lag), lag_i = length(d{i}.OutputData); else, lag_i = lag; end
    GFRF_cell{i} = interp(etfe(d{i}, lag_i, n_etfe), w);
    GFRF_cell{i}.Name      = sprintf('OP%d', i);
    GFRF_cell{i}.InputName = sprintf('u%d', i);

    GFRF(1, i) = frd(GFRF_cell{i}.ResponseData, GFRF_cell{1}.Frequency, Ts); %#ok<SAGROW>
end
bo=bodeoptions;
bo.FreqUnits='Hz';
figure; bode(GFRF_cell{:},bo); grid on; legend('location','best')
xlim([10 500])
exportgraphics(gcf,fullfile('plots','GFRF.jpg'),'Resolution',300)
%% Save
if save_output
    save(outfile, 'GFRF_cell', 'GFRF', 'd', 'realmass')
end
