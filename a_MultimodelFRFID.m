%% a_MultimodelFRFID  -  Step 1: non-parametric multimodel FRF identification
%
%   One sine-sweep experiment per operating point (added mass) -> one FRF
%   per operating point, estimated with spa.
%
%   Output (when save_output = true): Data/GFRF_iddata_1028_spa.mat
%     GFRF      1 x nOP frd (channel i = operating point i)
%     GFRF_cell 1 x nOP cell of idfrd (same FRFs)
%     d         1 x nOP cell of iddata (detrended time data)
%     realmass  added mass of each operating point (kg)
%
clear; close all; clc

%% Config
idfolder = fullfile('Data', 'MIVIDA 261028');
files    = {'261028_gfrf0kg_xsswp.mat', '261028_gfrf2.1kg_xsswp.mat', ...
            '261028_gfrf5.9kg_xsswp.mat', '261028_gfrf7.6kg_xsswp.mat'};
realmass = [0 2.1 5.9 7.6];          % added mass of each file (kg)
itmax    = {60000, [], [], 65000};   % last sample kept per file ([] = whole record)

Ts  = 0.001;
w   = logspace(log10(2*pi*10), log10(pi/Ts), 500);   % rad/s, 10 Hz .. Nyquist
lag = [];      % spa lag window; [] = record length (no smoothing), smaller = smoother

save_output = false;   % true: (over)write the output file below
outfile     = fullfile('Data', 'GFRF_iddata_1028_spa.mat');

%% FRF of each operating point
nOP       = numel(files);
GFRF_cell = cell(1, nOP);
d         = cell(1, nOP);
for i = 1:nOP
    [u, y] = utils.openIOtemp(fullfile(idfolder, files{i}), 1, itmax{i});

    d{i} = iddata(y(:), u(:), Ts, 'InputName', 'u'+string(i), 'OutputName', 'y');
    d{i} = detrend(d{i});                          % remove the operating-point offset

    if isempty(lag), lag_i = length(d{i}.OutputData); else, lag_i = lag; end
    GFRF_cell{i} = spa(d{i}, lag_i, w);
    GFRF_cell{i}.Name      = sprintf('OP%d', i);
    GFRF_cell{i}.InputName = sprintf('u%d', i);

    GFRF(1, i) = frd(GFRF_cell{i}.ResponseData, GFRF_cell{1}.Frequency, Ts); %#ok<SAGROW>
end

figure; bode(GFRF_cell{:}); grid on; legend show

%% Save
if save_output
    save(outfile, 'GFRF_cell', 'GFRF', 'd', 'realmass')
end
