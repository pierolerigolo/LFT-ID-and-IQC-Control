function T = nyquistcheck(K, GFRF, opts)
%NYQUISTCHECK Nyquist plot of L = G K on every operating point, closest point to -1
%
%   T = utils.nyquistcheck(K, GFRF)
%   T = utils.nyquistcheck(K, GFRF, 'G_hat', G_hat, 'realmass', m, 'name', "W2 order 18")
%
%   GFRF  measured FRFs, one channel per operating point (as in the c scripts)
%   For each operating point (and the nominal model G_hat if given) the
%   frequency where L(jw) = G(jw) K(jw) is closest to -1 is marked on the
%   Nyquist plot (full plot and zoom around -1), and returned in T:
%     f_Hz        frequency of the closest point
%     min_1pL     min |1 + L|  (modulus margin; 1/min_1pL = peak |S|)
%     peak_S_dB   -20 log10(min_1pL)
%     angL_deg    phase of L at that point
%     N_enc       encirclements of -1 by L (counter-clockwise > 0), from the winding
%                 of 1 + L over -w_max..w_max (positive frequencies and their mirror)
%     CL_stable   Nyquist: with K and G stable, closed loop stable iff N_enc = 0
%                 (K unstable: |N_enc| must equal its number of unstable poles)
%   The frequency grid is refined around the lightly damped poles of K (zeta < 0.05),
%   so narrow controller peaks between FRF points are not missed.
%   Below the first FRF frequency L is assumed small (warning if |L| > 0.5 there).
%   'file', 'plots/x' saves plots/x.pdf and plots/x.csv.
%   Only positive frequencies are drawn. The FRFs are interpolated (interp, linear)
%   on n_interp points between their first and last measured frequency (no
%   extrapolation); K is evaluated exactly on that grid. 'n_interp', 0 = measured grid.
%
arguments
    K
    GFRF
    opts.G_hat = []
    opts.realmass double = []
    opts.name string = ""
    opts.zoom (1,1) double = 1.5      % half-width of the zoom window around -1
    opts.n_interp (1,1) double = 5000  % FRFs interpolated on this many linear points (0 = measured grid)
    opts.file string = ""            % if given: <file>.pdf (figure) and <file>.csv (table)
end
    nOP = size(GFRF, 2);
    labels = "OP" + string((1:nOP).');
    if numel(opts.realmass) == nOP
        labels = labels + " (" + compose("%.2f", opts.realmass(:)) + " kg)";
    end
    Gs = cell(nOP, 1);
    for i = 1:nOP, Gs{i} = GFRF(i); end
    if ~isempty(opts.G_hat)
        Gs{end+1} = opts.G_hat;   labels(end+1) = "Nominal";
    end
    n = numel(Gs);
    col = lines(max(n, 3));   col(n, :) = 0;     % nominal (last) in black when present
    if isempty(opts.G_hat), col = lines(max(n, 3)); end

    f_c = zeros(n, 1);   m_c = f_c;   a_c = f_c;   N_enc = f_c;

    % dense points around the lightly damped poles of K (narrow peaks)
    pK  = pole(K);
    if K.Ts > 0, sK = log(pK)/K.Ts; else, sK = pK; end
    wn  = abs(sK);   zK = -real(sK)./max(wn, eps);
    sel = imag(sK) > 0 & zK < 0.05;
    w_K = [];
    for q = find(sel).'
        w_K = [w_K; wn(q)*linspace(1 - 30*max(zK(q), 1e-4), 1 + 30*max(zK(q), 1e-4), 2000).']; %#ok<AGROW>
    end
    if K.Ts > 0, P = sum(abs(pK) > 1); else, P = sum(real(pK) > 0); end
    fig = figure('Name', 'Nyquist L = G K');
    tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
    ax  = [nexttile(tl), nexttile(tl)];
    for a = ax, hold(a, 'on'); grid(a, 'on'); axis(a, 'equal'); end
    for i = 1:n
        Gi = Gs{i};
        if opts.n_interp > 0
            wg = linspace(min(Gi.Frequency), max(Gi.Frequency), opts.n_interp).';
            wg = unique([wg; w_K(w_K >= min(Gi.Frequency) & w_K <= max(Gi.Frequency))]);
            Gi = interp(Gi, wg);
        end
        w = Gi.Frequency;
        L = squeeze(freqresp(K, w)).*squeeze(Gi.ResponseData);
        ok = isfinite(L);   w = w(ok);   L = L(ok);
        [m_c(i), k] = min(abs(1 + L));
        f_c(i) = w(k)/2/pi;   a_c(i) = rad2deg(angle(L(k)));
        % encirclements of -1: winding of 1 + L over -w_max..w_max
        ph = unwrap(angle(1 + [conj(flipud(L)); L]));
        N_enc(i) = round((ph(end) - ph(1))/(2*pi));
        if abs(L(1)) > 0.5
            warning('nyquistcheck:lowfreq', '%s: |L| = %.2f at the first FRF frequency (%.1f Hz): the encirclement count assumes L small below it.', ...
                labels(i), abs(L(1)), w(1)/2/pi);
        end
        for a = ax
            plot(a, real(L), imag(L), 'Color', col(i, :), 'DisplayName', labels(i));
            plot(a, [-1 real(L(k))], [0 imag(L(k))], '--', 'Color', col(i, :), 'HandleVisibility', 'off');
            plot(a, real(L(k)), imag(L(k)), 'o', 'Color', col(i, :), 'MarkerFaceColor', col(i, :), ...
                'HandleVisibility', 'off');
        end
        text(ax(2), real(L(k)), imag(L(k)), sprintf('  %.1f Hz, |1+L| = %.2f', f_c(i), m_c(i)), ...
            'Color', col(i, :), 'FontSize', 8, 'Clipping', 'on');
    end
    th = linspace(0, 2*pi, 200);
    for a = ax
        plot(a, -1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2, 'HandleVisibility', 'off');
        plot(a, -1 + min(m_c)*cos(th), min(m_c)*sin(th), 'r:', 'HandleVisibility', 'off');
        xlabel(a, 'Re L'); ylabel(a, 'Im L');
    end
    xlim(ax(2), -1 + opts.zoom*[-1 1]);   ylim(ax(2), opts.zoom*[-1 1]);
    title(ax(1), 'L = G K');   title(ax(2), 'zoom around -1 (circle: smallest |1+L|)');
    legend(ax(1), 'Location', 'best');
    if strlength(opts.name) > 0, title(tl, opts.name, 'Interpreter', 'none'); end

    % stable K (P = 0): stable closed loop iff no net encirclement of -1.
    % unstable K: |N_enc| = P required (direction not checked here)
    CL_stable = abs(N_enc) == P;
    if P > 0, warning('nyquistcheck:Kunstable', 'K has %d unstable poles: check the direction of the encirclements.', P); end
    T = table(f_c, m_c, -mag2db(m_c), a_c, N_enc, CL_stable, ...
        'VariableNames', {'f_Hz', 'min_1pL', 'peak_S_dB', 'angL_deg', 'N_enc', 'CL_stable'}, 'RowNames', cellstr(labels));
    disp(T)
    if strlength(opts.file) > 0
        [d, ~] = fileparts(opts.file);
        if strlength(d) > 0 && ~exist(d, 'dir'), mkdir(d); end
        fig.Position(3:4) = [1100 500];
        exportgraphics(fig, opts.file + ".pdf", 'ContentType', 'vector');
        writetable(T, opts.file + ".csv", 'WriteRowNames', true);
        fprintf('Nyquist check saved: %s.pdf / .csv\n', opts.file);
    end
end
