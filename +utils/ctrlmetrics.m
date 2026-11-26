function T = ctrlmetrics(Ks, names, GFRF, W1, omega_d, opts)
%CTRLMETRICS Comparison metrics of controllers on the measured operating points
%
%   T = utils.ctrlmetrics(Ks, names, GFRF, W1, omega_d)
%   T = utils.ctrlmetrics(..., 'file', 'plots/x.csv', 'S_notch_dB', -20)
%
%   Ks, names : controllers and their labels;  GFRF: one channel per operating point
%   W1        : performance weight (frd or lti);  omega_d: disturbance frequency (rad/s)
%   One row per controller, over the operating points (worst = largest value):
%     S80_worst_dB / S80_spread_dB : |S(omega_d)| worst and max - min across operating points
%     W1S_worst, W1S_worst_OP, W1S_worst_Hz : max_i ||W1 S_i||_inf (what the robust
%                                    performance gamma bounds), where and at which frequency
%     Speak_dB, min_1pL           : peak |S| over frequency and operating points, and min |1+L|
%     notch_width_Hz              : width of the band around omega_d where |S| < S_notch_dB
%                                    on EVERY operating point (tolerance to a frequency drift)
%     KS80_max_dB, KSpeak_dB      : control effort |K S| at omega_d (worst) and its peak
%     N_enc_total                 : total encirclements of -1 (0 = stable on all FRFs)
%   The FRFs are interpolated on a fine grid (dense around omega_d) between their first
%   and last measured frequency. Printed, written to 'file' (csv) and copied to the
%   clipboard (tab separated).
%
arguments
    Ks cell
    names string
    GFRF
    W1
    omega_d (1,1) double
    opts.file string = ""
    opts.S_notch_dB (1,1) double = -20
    opts.n_grid (1,1) double = 6000
end
    nOP = size(GFRF, 2);
    w0  = min(GFRF.Frequency);   w1 = max(GFRF.Frequency);
    w   = unique([linspace(w0, w1, opts.n_grid).'; ...
                  omega_d + 2*pi*(-10:0.005:10).']);              % dense around omega_d
    w   = w(w >= w0 & w <= w1);
    Gw  = zeros(numel(w), nOP);
    for i = 1:nOP, Gw(:, i) = squeeze(interp(GFRF(i), w).ResponseData); end
    if isa(W1, 'frd'), W1w = squeeze(interp(W1, w).ResponseData); else, W1w = squeeze(freqresp(W1, w)); end
    [~, k80] = min(abs(w - omega_d));

    nK = numel(Ks);
    V  = nan(nK, 11);
    for j = 1:nK
        Kw = squeeze(freqresp(Ks{j}, w));
        L  = Kw.*Gw;   S = 1./(1 + L);   KS = Kw.*S;
        S80 = mag2db(abs(S(k80, :)));
        W1S = abs(W1w.*S);   [w1s, kk] = max(W1S(:));   [kw, iop] = ind2sub(size(W1S), kk);
        below = all(mag2db(abs(S)) < opts.S_notch_dB, 2);       % band where every OP is below
        width = 0;
        if below(k80)
            a = find(~below(1:k80), 1, 'last');   b = k80 - 1 + find(~below(k80:end), 1);
            if isempty(a), a = 0; end
            if isempty(b), b = numel(w) + 1; end
            width = (w(min(b - 1, end)) - w(max(a + 1, 1)))/2/pi;
        end
        V(j, :) = [max(S80), max(S80) - min(S80), w1s, iop, w(kw)/2/pi, ...
                   mag2db(max(abs(S(:)))), min(abs(1 + L(:))), width, ...
                   mag2db(max(abs(KS(k80, :)))), mag2db(max(abs(KS(:)))), ...
                   sum(abs(utils.encirclements(ss(Ks{j}), GFRF)))];
    end
    vn = {'S80_worst_dB', 'S80_spread_dB', 'W1S_worst', 'W1S_worst_OP', 'W1S_worst_Hz', ...
          'Speak_dB', 'min_1pL', 'notch_width_Hz', 'KS80_max_dB', 'KSpeak_dB', 'N_enc_total'};
    T = array2table(round(V, 3), 'VariableNames', vn, 'RowNames', cellstr(names(:)));
    fprintf('Controller comparison over %d operating points (notch width: |S| < %g dB on all of them)\n', ...
        nOP, opts.S_notch_dB);
    disp(T)
    if strlength(opts.file) > 0
        writetable(T, opts.file, 'WriteRowNames', true);
    end
    lines = ["Controller", string(vn); names(:), compose("%.3g", V)];
    try
        clipboard('copy', char(strjoin(join(lines, char(9), 2), newline)));
        fprintf('Metrics copied to the clipboard (tab separated).\n');
    catch
    end
end
