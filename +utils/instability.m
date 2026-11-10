function r = instability(path, Fs, do_plot)
%INSTABILITY Analyse a closed-loop run that diverges and saturates
%
%   r = utils.instability(path)            Fs = 1000 Hz, with plots
%   r = utils.instability(path, Fs, false) without plots
%
%   Finds when the controller output becomes active, fits the exponential
%   growth of the acceleration envelope (from 10x the noise level to 30 % of
%   the plateau) and the decay after the peak, and checks the plateau for
%   clipping. Fields of r:
%     t_on            controller output becomes active (s)
%     f_growth        oscillation frequency while it grows (Hz)
%     sigma_growth    envelope growth rate (1/s), > 0 = unstable closed loop
%     zeta_cl         equivalent closed-loop damping  = -sigma/(2 pi f)
%     f_decay, sigma_decay, zeta_decay   same for the decay after the peak
%     u_amp, y_amp    plateau amplitudes of controller output / acceleration
%     u_crest, y_crest  peak/rms on the plateau (pure sine 1.414, clipped < 1.414)
%     u_harm3_dB      3rd harmonic of u relative to the fundamental (plateau)
%     u_clipped       true if u looks hard-clipped (crest < 1.38)
%
arguments
    path
    Fs (1,1) double = 1000
    do_plot (1,1) logical = true
end
    [y, ~, u] = utils.loadexp(path);
    t   = (0:numel(y)-1).'/Fs;
    env = movmax(abs(y), round(Fs/20));               % 50 ms envelope
    lvl = movmedian(env, Fs);                          % 1 s smoothed level

    k_on  = find(abs(u) > 1e-3*max(abs(u)), 1);
    noise = median(env(1:max(k_on - 1, 1)));
    [pl, k_pk] = max(lvl);

    % ----- growth -----
    g = find(env > 10*noise & env < 0.3*pl);   g = g(g > k_on & g < k_pk);
    [sg, fg] = fit_env(g);
    % ----- decay after the peak -----
    d = find(env > 10*noise & env < 0.3*pl);   d = d(d > k_pk);
    [sd, fd] = fit_env(d);
    % ----- plateau -----
    p = find(lvl > 0.8*pl);
    up = u(p) - mean(u(p));   yp = y(p) - mean(y(p));
    [~, Xu, f] = dominant(up);
    [~, i1] = max(Xu(2:end));   i1 = i1 + 1;
    [~, i3] = min(abs(f - 3*f(i1)));

    r = struct('t_on', t(k_on), 'f_growth', fg, 'sigma_growth', sg, 'zeta_cl', -sg/(2*pi*fg), ...
        'f_decay', fd, 'sigma_decay', sd, 'zeta_decay', -sd/(2*pi*fd), ...
        'u_amp', max(abs(up)), 'y_amp', max(abs(yp)), ...
        'u_crest', max(abs(up))/rms(up), 'y_crest', max(abs(yp))/rms(yp), ...
        'u_harm3_dB', 20*log10(max(Xu(max(i3-2, 1):min(i3+2, end)))/Xu(i1)));
    r.u_clipped = r.u_crest < 1.38;

    if do_plot
        figure('Name', path)
        subplot(3, 1, 1); plot(t, y); hold on; plot(t, env, 'k');
        if ~isempty(g), plot(t(g), exp(polyval(polyfit(t(g), log(env(g)), 1), t(g))), 'r', 'LineWidth', 2); end
        xline(t(k_on), '--', 'controller on'); ylabel('acc. (m/s^2)'); grid on
        title(sprintf('growth %.1f Hz, \\sigma = %.2f 1/s (\\zeta_{cl} = %.4f)', fg, sg, r.zeta_cl))
        subplot(3, 1, 2); plot(t, u); xline(t(k_on), '--'); ylabel('controller output'); grid on
        title(sprintf('plateau |u| = %.3g, crest %.2f (sine 1.41), 3rd harmonic %.0f dB', ...
            r.u_amp, r.u_crest, r.u_harm3_dB))
        subplot(3, 1, 3);
        if ~isempty(g), [~, Xg, fgv] = dominant(y(g(1):g(end))); semilogy(fgv, Xg); end
        xlim([1 Fs/2]); grid on; xlabel('Frequency (Hz)'); ylabel('|Y| growth window')
    end

    function [sigma, fdom] = fit_env(idx)
        sigma = NaN;   fdom = NaN;
        if numel(idx) < Fs/20, return, end
        c = polyfit(t(idx), log(env(idx)), 1);   sigma = c(1);
        fdom = dominant(y(idx(1):idx(end)));
    end

    function [fdom, X, fv] = dominant(x)
        % dominant frequency of x (Hann window, parabolic interpolation)
        x = x(:) - mean(x);   N = numel(x);
        X  = abs(fft(x.*hann(N)));   X = X(1:floor(N/2) + 1);
        fv = (0:floor(N/2)).'*Fs/N;
        [~, k] = max(X(2:end));   k = k + 1;
        if k > 1 && k < numel(X)
            a = log(X(k-1));  b = log(X(k));  cc = log(X(k+1));
            fdom = (k - 1 + 0.5*(a - cc)/(a - 2*b + cc))*Fs/N;
        else
            fdom = fv(k);
        end
    end
end
