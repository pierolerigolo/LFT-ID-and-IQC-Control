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
%     Data/Experiments/<case folder>/sine/<mass>kg*.mat
%     Data/Experiments/<case folder>/sweep/<mass>kg*.mat
%   e.g. Data/Experiments/80hz_gdmw2/sine/2.1kg.mat and .../sweep/2.1kg.mat
%
%   Output: one table per case (command window), plots/d_<case>_S.jpg and
%           plots/d_experimental_attenuation.tex
%
%   Unstable runs (instab_folder, top-level .mat files): growth frequency and
%   rate, decay, plateau clipping (utils.instability). With the controller of
%   the run (contfn, or k_final_<controller file> in the run's file name):
%   its poles near the growth frequency, the recorded controller u/y compared
%   with -K_final at that frequency, and the lsim check (K_final applied to the
%   recorded acceleration vs the recorded controller output). The designed
%   controller is evaluated at s = sigma + j*omega (the oscillation grows).
%
%   Sweep runs (sweep_runs): implemented vs designed controller over the whole
%   band (-u/y in the closed-loop sweep vs K), estimated implementation gain,
%   and measured S = Y_cl/Y_ol vs 1/(1+K Gu) and 1/(1+gain K Gu), with Gu the
%   controller-output -> acceleration path measured from the sweeps.
%
%   Nyquist (nyq_K): loop on every measured FRF as designed and with the
%   implementation gain (gain_impl, or the gain estimated from the sweep runs).
%
clear; close all; clc

%% Config
run_cases     = false;                         % attenuation analysis of the cases below
run_instab    = true;                          % analysis of unstable runs
run_sweepchk  = true;                          % implemented vs designed controller on sweep runs
run_nyquist   = true;                          % Nyquist of the loop with the implementation gain
instab_folder = fullfile('Data', 'Experiments');
contfn        = 'LFT_GFRF_GdmGLTIW2opt_1028_etfe_80hz_fmin10_hyb_controller_order18_80hz_-34db_nperfoptiIQC_notconv_etai0.24_etaf0.60_multi_freq.mat';                            % controller file (Data/controller) used in the unstable
                                               % runs; '' = from the run's file name after k_final_
expfolder = fullfile('Data', 'Experiments');
trans     = 3000;      % samples skipped after each switch in the sine experiments
Ts        = 0.001;

% Sweep runs (open-loop then closed-loop disturbance sweep) and their controller
% (file, and itr: [] = K_final, n = K_itr{n}, Inf = last K_itr)
sweep_runs = struct( ...
    'file',       {fullfile('Data', 'Experiments', 'Old Experiments', ...
                   'controller_sus_7nm_order18_82hz_-34db_nperfhinfnperfoptiIQC_eta9.139390e-01_itr26', 'gfrf0_sswep.mat')}, ...
    'controller', {fullfile('Data', 'finalrepportcont', 'controller_sus_7nm_order18_82hz_-34db_nperfhinfnperfoptiIQC_eta9.139390e-01.mat')}, ...
    'itr',        {26});
f_chk = (15:5:450).';   % frequencies of the controller / sensitivity comparison (Hz)

% Nyquist: controllers to compare (same fields), plant = measured FRFs of nyq_frf
nyq_K = struct( ...
    'name',       {'K chosen (Gdm-W2)',                         'KSOA (ellipse)'}, ...
    'controller', {fullfile('Data', 'controller', contfn),     fullfile('Data', 'finalrepportcont', 'controller_ellips_sus_7nm_order20_82hz_-30db_nomperf.mat')}, ...
    'itr',        {[],                                          Inf});
nyq_frf   = 'LFT_GFRF_GdmGLTIW2opt_1028_etfe_80hz_fmin10_hyb';   % Data/<nyq_frf>.mat holds GFRF
gain_impl = [];         % implementation gain; [] = estimated from the sweep runs

% folder      : case folder in expfolder
% controller  : file in Data/controller (K_final, or K_itr{itr} if itr is set)
% f_dist      : disturbance frequency (Hz)
cases = struct( ...
    'folder',     {'80hz_gdmw2',                    '80hz_elli',                     '50hz'}, ...
    'label',      {'80 Hz $G_{dm}$-$W_2$',          '80 Hz ellipse',                 '50 Hz $G_{dm}$-$W_2$'}, ...
    'controller', {'<80 Hz W2 controller>.mat',     '<80 Hz ellipse controller>.mat', '<50 Hz W2 controller>.mat'}, ...
    'itr',        {[],                              [],                              []}, ...
    'f_dist',     {80,                              80,                              50});

if ~exist('plots', 'dir'), mkdir('plots'); end

%% Unstable runs
if run_instab
    files = dir(fullfile(instab_folder, '*.mat'));
    for i = 1:numel(files)
        fprintf('\n== %s ==\n', files(i).name);
        r = utils.instability(fullfile(files(i).folder, files(i).name));
        fprintf(['controller on at %.2f s\n' ...
            'growth : %.2f Hz, sigma = %+.2f 1/s, closed-loop damping %.4f\n' ...
            'decay  : %.2f Hz, sigma = %+.2f 1/s, damping %.4f\n' ...
            'plateau: |u| = %.3g (crest %.2f, 3rd harm. %.0f dB, clipped %d), |y| = %.3g (crest %.2f)\n'], ...
            r.t_on, r.f_growth, r.sigma_growth, r.zeta_cl, r.f_decay, r.sigma_decay, r.zeta_decay, ...
            r.u_amp, r.u_crest, r.u_harm3_dB, r.u_clipped, r.y_amp, r.y_crest);
        cfile = contfn;
        if isempty(cfile)
            tok = regexp(files(i).name, 'k_final_(.*)$', 'tokens', 'once');
            if ~isempty(tok), cfile = tok{1}; end
        end
        if isempty(cfile) || ~exist(fullfile('Data', 'controller', cfile), 'file')
            fprintf('controller file not found (set contfn): no controller check\n');
            continue
        end
        Kf = load(fullfile('Data', 'controller', cfile), 'K_final');   K = Kf.K_final;
        [wn, z] = damp(K);
        sel = abs(wn/2/pi - r.f_growth) < 0.2*r.f_growth;
        fprintf('controller poles near %.0f Hz: %s Hz, damping %s\n', r.f_growth, ...
            mat2str(wn(sel).'/2/pi, 5), mat2str(z(sel).', 3));

        % --- implemented vs designed controller (u = -K y) ---
        if isempty(r.idx_growth)
            fprintf('no growth window found: no controller check\n');
            continue
        end
        [y, ~, u] = utils.loadexp(fullfile(files(i).folder, files(i).name));
        t   = (0:numel(y)-1).'*Ts;
        idx = r.idx_growth(1):r.idx_growth(2);
        Hrec = band_ratio(u(idx), y(idx), r.f_growth, 1/Ts);          % recorded u/y
        % the oscillation grows as exp((sigma + j w) t): compare with K at s = sigma + j w
        % (K(j w) gives a fake phase shift next to a lightly damped controller pole)
        s_g  = r.sigma_growth + 1j*2*pi*r.f_growth;
        Hdes = -evalfr(K, exp(s_g*Ts));                               % designed  -K
        fprintf(['at %.2f Hz (sigma %+.2f): recorded u/y %6.1f dB %7.1f deg | designed -K %6.1f dB %7.1f deg' ...
            ' | difference %+.2f dB %+.1f deg\n'], r.f_growth, r.sigma_growth, ...
            mag2db(abs(Hrec)), rad2deg(angle(Hrec)), mag2db(abs(Hdes)), rad2deg(angle(Hdes)), ...
            mag2db(abs(Hrec/Hdes)), rad2deg(angle(Hrec/Hdes)));

        usim = lsim(K, -y, t);                                         % K_final on the recorded y
        % (with an implementation gain g, the recorded output is g times this curve)
        figure('Name', ['lsim check: ' files(i).name])
        plot(t, u, t, usim, '--'); grid on
        xlim(t(idx([1 end])).' + [-0.05 0.05])
        legend('recorded controller output', 'K\_final applied to the recorded y', 'Location', 'best')
        xlabel('time (s)'); ylabel('controller output')
        title(sprintf('growth window, %.1f Hz', r.f_growth))
    end
end

%% Sweep runs: implemented vs designed controller, gain estimate, sensitivity
g_est = [];
for i = 1:numel(sweep_runs)*run_sweepchk
    R = sweep_runs(i);
    K = load_K(R.controller, R.itr);
    [y, d, u] = utils.loadexp(R.file);
    [iol, icl, L] = sweep_windows(d, u, 1000);
    F  = @(x) fft(x - mean(x));
    Yo = F(y(iol:iol+L-1));  Do = F(d(iol:iol+L-1));
    Yc = F(y(icl:icl+L-1));  Dc = F(d(icl:icl+L-1));  Uc = F(u(icl:icl+L-1));
    fv = (0:L-1).'/(L*Ts);
    n  = numel(f_chk);
    [ratio, Gd, Gu, Sexp] = deal(complex(zeros(n, 1)));
    for k = 1:n
        m = abs(fv - f_chk(k)) < 1;                                  % +/- 1 Hz band
        Gd(k)    = sum(Yo(m).*conj(Do(m)))/sum(abs(Do(m)).^2);      % disturbance -> y (open loop)
        rr       = Yc(m) - Gd(k)*Dc(m);                              % closed loop: y = Gd d + Gu u
        Gu(k)    = sum(rr.*conj(Uc(m)))/sum(abs(Uc(m)).^2);         % controller output -> y
        Himpl    = -sum(Uc(m).*conj(Yc(m)))/sum(abs(Yc(m)).^2);     % implemented controller -u/y
        ratio(k) = Himpl/evalfr(K, exp(1j*2*pi*f_chk(k)*Ts));        % implemented / designed
        Sexp(k)  = (sum(Yc(m).*conj(Dc(m)))/sum(abs(Dc(m)).^2))/Gd(k);
    end
    g = median(abs(ratio))*sign(median(cos(angle(ratio))));          % robust to the controller peak
    g_est(end+1) = g; %#ok<AGROW>
    [~, nm] = fileparts(R.file);
    fprintf(['\n== sweep run %s ==\nimplemented / designed controller: median %.4f (%+.2f dB), ' ...
        'phase median %+.1f deg (1 sample of delay would be -360 f Ts)\n'], ...
        nm, g, mag2db(abs(g)), rad2deg(median(angle(ratio))));
    Kw   = arrayfun(@(f) evalfr(K, exp(1j*2*pi*f*Ts)), f_chk);
    Sth  = 1./(1 + Kw.*Gu);
    Sthg = 1./(1 + g*Kw.*Gu);
    disp(table(f_chk, mag2db(abs(ratio)), rad2deg(angle(ratio)), mag2db(abs(Sexp)), mag2db(abs(Sth)), mag2db(abs(Sthg)), ...
        'VariableNames', {'f_Hz', 'impl_over_K_dB', 'impl_over_K_deg', 'S_meas_dB', 'S_th_K_dB', 'S_th_gK_dB'}))

    figure('Name', ['sweep check: ' nm])
    subplot(3, 1, 1); semilogx(f_chk, mag2db(abs(ratio)), 'o-'); grid on
    yline(mag2db(abs(g)), '--', sprintf('%.2f dB', mag2db(abs(g)))); ylabel('|impl / K| (dB)')
    title('implemented controller (-u/y, closed-loop sweep) / designed K')
    subplot(3, 1, 2); semilogx(f_chk, rad2deg(angle(ratio)), 'o-', f_chk, -360*f_chk*Ts, '--'); grid on
    ylabel('phase (deg)'); legend('impl / K', '1 sample delay', 'Location', 'best')
    subplot(3, 1, 3); semilogx(f_chk, mag2db(abs(Sexp)), 'k-', f_chk, mag2db(abs(Sth)), 'b--', f_chk, mag2db(abs(Sthg)), 'r--'); grid on
    ylabel('|S| (dB)'); xlabel('Frequency (Hz)')
    legend('measured Y_{cl}/Y_{ol}', '1/(1+K G_u)', sprintf('1/(1+%.2f K G_u)', g), 'Location', 'best')
end

%% Nyquist: loop as designed vs with the implementation gain
if run_nyquist
    if isempty(gain_impl)
        if isempty(g_est), error('Set gain_impl or enable run_sweepchk to estimate it.'); end
        gain_impl = median(g_est);
    end
    FRF = load(fullfile('Data', [nyq_frf '.mat']), 'GFRF');
    w   = unique([FRF.GFRF.Frequency(:); 2*pi*linspace(60, 100, 20001).']);   % + fine grid at 60-100 Hz
    w   = w(w >= min(FRF.GFRF.Frequency) & w <= max(FRF.GFRF.Frequency));
    nOP = size(FRF.GFRF, 2);
    col = lines(nOP);
    gains = [1, gain_impl];
    figure('Name', 'Nyquist: designed vs implemented gain')
    for c = 1:numel(nyq_K)
        Kw = squeeze(freqresp(load_K(nyq_K(c).controller, nyq_K(c).itr), w));
        for v = 1:2
            subplot(numel(nyq_K), 2, 2*(c - 1) + v); hold on; grid on
            for i = 1:nOP
                Lw = gains(v)*Kw.*squeeze(interp(FRF.GFRF(i), w).ResponseData);
                [mm, k] = min(abs(1 + Lw));
                fprintf('%-18s gain %6.3f OP%d: min|1+L| = %.3f at %.2f Hz\n', nyq_K(c).name, gains(v), i, mm, w(k)/2/pi);
                plot(real(Lw), imag(Lw), 'Color', col(i, :), 'DisplayName', 'OP'+string(i))
            end
            plot(-1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2, 'HandleVisibility', 'off')
            axis equal; xlim([-6 4]); ylim([-5 5])
            title(sprintf('%s, gain %.3g', nyq_K(c).name, gains(v)))
        end
    end
    legend('Location', 'best')
    exportgraphics(gcf, fullfile('plots', 'd_nyquist_gain.png'), 'Resolution', 300)
end

%% Analysis of the cases
nC  = numel(cases);
res = cell(nC, 1);
for c = 1:nC*run_cases
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
if ~run_cases, return, end
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
function K = load_K(file, itr)
%LOAD_K Controller from a saved run: [] = K_final, n = K_itr{n}, Inf = last K_itr
    if isempty(itr)
        S = load(file, 'K_final');   K = S.K_final;
    else
        S = load(file, 'K_itr');
        if isinf(itr), itr = numel(S.K_itr); end
        K = S.K_itr{itr};
    end
end

function [iol, icl, L] = sweep_windows(d, u, pre)
%SWEEP_WINDOWS Start of the open-loop and closed-loop sweeps (same length L),
%   from the disturbance command d and the controller output u.
    m  = sqrt(movmean(d.^2, 500));   on = m > 0.05*max(m);
    st = find(diff([0; on]) == 1);   en = find(diff([on; 0]) == -1);
    full = (en - st) >= 0.8*max(en - st);   st = st(full);   en = en(full);
    ru  = arrayfun(@(k) rms(u(st(k):en(k))), 1:numel(st));
    isc = ru > 0.1*max(ru);
    kol = find(~isc, 1);   kcl = find(isc, 1);
    if isempty(kol) || isempty(kcl), error('Need one open-loop and one closed-loop sweep.'); end
    iol = max(st(kol) - pre, 1);   icl = max(st(kcl) - pre, 1);
    L   = min([en(kol) - st(kol), en(kcl) - st(kcl)]) + pre;
    L   = min([L, numel(d) - iol + 1, numel(d) - icl + 1]);
end

function H = band_ratio(num, den, f0, Fs)
%BAND_RATIO num/den at f0 (Hann window, cross-spectrum over +/- 3 Hz)
    w = hann(numel(num));
    N = fft((num(:) - mean(num)).*w);   D = fft((den(:) - mean(den)).*w);
    f = (0:numel(num)-1).'*Fs/numel(num);
    m = abs(f - f0) < 3;
    H = sum(N(m).*conj(D(m)))/sum(abs(D(m)).^2);
end

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
