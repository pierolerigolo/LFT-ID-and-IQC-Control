%% d_ExperimentalResults  -  Step 4: closed-loop experimental results
%
%   For every case (controller) and every operating point (added mass):
%     sine  : attenuation of the single tone and SNR        (utils.sineatt)
%     sweep : experimental sensitivity S = Y_cl/Y_ol         (utils.sweepS)
%     theory: S = 1/(1 + K G_i) on the measured FRF of that mass, if it is
%             one of the identified operating points
%   all read at the disturbance frequency of the case.
%
%   Data layout (segments are found automatically from the recorded signals):
%     Data/experimental/<case folder>/sine/<mass>kg*.mat
%     Data/experimental/<case folder>/sweep/<mass>kg*.mat
%   e.g. Data/experimental/80hz_gdmw2/sine/2.1kg.mat and .../sweep/2.1kg.mat
%
%   Output: one table per case (command window), plots/d_<case>_S.jpg and
%           plots/d_experimental_attenuation.tex
%
clear; close all; clc

%% Config
expfolder = fullfile('Data', 'experimental');
trans     = 3000;      % samples skipped after each switch in the sine experiments
Ts        = 0.001;

% folder      : case folder in expfolder
% controller  : file in Data/controller (K_final, or K_itr{itr} if itr is set)
% f_dist      : disturbance frequency (Hz)
cases = struct( ...
    'folder',     {'80hz_gdmw2',                    '80hz_elli',                     '50hz'}, ...
    'label',      {'80 Hz $G_{dm}$-$W_2$',          '80 Hz ellipse',                 '50 Hz $G_{dm}$-$W_2$'}, ...
    'controller', {'<80 Hz W2 controller>.mat',     '<80 Hz ellipse controller>.mat', '<50 Hz W2 controller>.mat'}, ...
    'itr',        {[],                              [],                              []}, ...
    'f_dist',     {80,                              80,                              50});

%% Analysis
if ~exist('plots', 'dir'), mkdir('plots'); end
nC  = numel(cases);
res = cell(nC, 1);
for c = 1:nC
    C    = cases(c);
    dirC = fullfile(expfolder, C.folder);
    ctrl = load(fullfile('Data', 'controller', C.controller));
    if isempty(C.itr), K = ctrl.K_final; else, K = ctrl.K_itr{C.itr}; end
    frf  = load(fullfile('Data', [ctrl.idname, '.mat']), 'GFRF', 'realmass');
    w0   = 2*pi*C.f_dist;

    [mass, fsine] = list_mass(fullfile(dirC, 'sine'));
    [mass_sw, fsweep] = list_mass(fullfile(dirC, 'sweep'));
    n = numel(mass);
    [sine, snr, sweep, theory] = deal(nan(n, 1));

    figure; hold on
    col = lines(n);
    for i = 1:n
        [sine(i), snr(i)] = utils.sineatt(fsine{i}, trans);

        j = find(abs(mass_sw - mass(i)) < 1e-6, 1);
        if ~isempty(j)
            S = utils.sweepS(fsweep{j}, Ts);
            sweep(i) = mag2db(abs(squeeze(interp(S, w0).ResponseData)));
            utils.magplot(S, 'Color', col(i, :), 'DisplayName', sprintf('%.2f kg exp.', mass(i)));
        end

        k = find(abs(frf.realmass - mass(i)) < 0.05, 1);
        if ~isempty(k)
            Sth = feedback(1, K*frf.GFRF(k));
            theory(i) = mag2db(abs(squeeze(interp(Sth, w0).ResponseData)));
            utils.magplot(Sth, 'Color', col(i, :), 'LineStyle', '--', 'DisplayName', sprintf('%.2f kg th.', mass(i)));
        end
    end
    utils.magplot(inv(ctrl.W1_frd), 'Color', 'k', 'DisplayName', 'W_1^{-1}');
    xline(C.f_dist, ':', 'HandleVisibility', 'off');
    title(sprintf('%s: sensitivity, experiment (-) and theory (--)', C.folder), 'Interpreter', 'none')
    exportgraphics(gcf, fullfile('plots', sprintf('d_%s_S.jpg', C.folder)), 'Resolution', 300)

    res{c} = table(mass, sine, sweep, theory, snr, ...
        'VariableNames', {'Mass_kg', 'Sine_dB', 'Sweep_dB', 'Theory_dB', 'SNR_dB'});
    fprintf('\n== %s (%g Hz) ==\n', C.folder, C.f_dist);
    disp(res{c})
end

%% LaTeX table: one row per mass, Sine / Sweep / Th. / SNR per case
allmass = unique(cell2mat(cellfun(@(T) T.Mass_kg, res, 'UniformOutput', false)));
fid = fopen(fullfile('plots', 'd_experimental_attenuation.tex'), 'w');
fprintf(fid, '\\begin{table}[h]\n\\centering\n\\begin{tabular}{l%s}\n\\toprule\n', repmat('cccc', 1, nC));
fprintf(fid, '%s \\\\\n', strjoin([{''}, arrayfun(@(C) sprintf('\\multicolumn{4}{c}{%s}', C.label), cases, 'UniformOutput', false)], ' & '));
fprintf(fid, '$\\Delta m$ [kg]%s \\\\\n\\midrule\n', repmat(' & Sine & Sweep & Th. ($S$) & SNR', 1, nC));
for m = allmass.'
    fprintf(fid, '%.2f', m);
    for c = 1:nC
        r = res{c}(abs(res{c}.Mass_kg - m) < 1e-6, :);
        v = nan(1, 4);
        if ~isempty(r), v = [r.Sine_dB, r.Sweep_dB, r.Theory_dB, r.SNR_dB]; end
        for x = v
            if isnan(x), fprintf(fid, ' & --'); else, fprintf(fid, ' & %.1f', x); end
        end
    end
    fprintf(fid, ' \\\\\n');
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fprintf(fid, ['\\caption{Experimental attenuation [dB] at the disturbance frequency for each added mass: ' ...
    'single tone (Sine), swept sine (Sweep), sensitivity on the measured FRF (Th., identified masses only) ' ...
    'and SNR of the closed-loop tone.}\n\\label{tab:exp_attenuation}\n\\end{table}\n']);
fclose(fid);

%% Local functions
function [mass, files] = list_mass(folder)
%LIST_MASS .mat files of a folder named <mass>kg*.mat, sorted by mass
    d = dir(fullfile(folder, '*.mat'));
    tok  = regexp({d.name}, '^(\d+(?:\.\d+)?)kg', 'tokens', 'once');
    keep = ~cellfun(@isempty, tok);
    d = d(keep);
    mass = cellfun(@(t) str2double(t{1}), tok(keep)).';
    [mass, o] = sort(mass);
    files = fullfile({d(o).folder}, {d(o).name}).';
end
