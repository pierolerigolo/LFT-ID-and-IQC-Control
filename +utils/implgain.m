function r = implgain(path, K, lags, skip)
%IMPLGAIN Gain (and delay) between the designed controller and the one that ran
%
%   r = utils.implgain(path, K)          lags -3..3 samples
%   r = utils.implgain(path, K, lags)
%   r = utils.implgain(path, K, lags, skip)   skip samples after switch-on in the fit
%                                             (default: 4 time constants of the slowest
%                                             pole of K, at most half of the segment)
%
%   Simulates the designed controller on the recorded acceleration of the
%   closed-loop part of the run (controller output active: moving RMS above 5 % of
%   its maximum), u_th = lsim(K, -y) from zero state, and fits the recorded
%   controller output:  u_rec(k) = g * u_th(k - lag).
%   Fields of r: gain (best g), lag (samples), residual (relative, best lag),
%   table (all lags), k_on (sample where the controller output becomes active),
%   t, u_rec, u_th (closed-loop part, for plots).
%   residual << 1 at lag 0 means u_rec = gain * K(-y): the controller that ran
%   is gain * K (pure gain, no delay).
%
arguments
    path
    K
    lags (1,:) double = -3:3
    skip double = []
end
    [y, ~, u] = utils.loadexp(path);
    m    = sqrt(movmean(u.^2, 500));
    on   = m > 0.05*max(m);                       % controller output active
    k_on = find(on, 1);   k_off = find(on, 1, 'last');
    idx  = (k_on:k_off).';
    t    = (idx - idx(1))*K.Ts;
    u_th = lsim(K, -y(idx), t);                   % from zero state: transient at the start
    if isempty(skip)
        tau  = -K.Ts/log(max(abs(pole(K))));      % slowest time constant (samples*Ts)
        skip = min(round(4*tau/K.Ts), floor(numel(idx)/2));
    end
    keep  = (skip + 1):numel(idx);
    t = t(keep);   u_th = u_th(keep);
    u_rec = u(idx(keep));

    g = zeros(size(lags));   res = g;
    for j = 1:numel(lags)
        a = u_rec(max(1, 1 + lags(j)):end + min(0, lags(j)));
        b = u_th(max(1, 1 - lags(j)):end + min(0, -lags(j)));
        g(j)   = (b.'*a)/(b.'*b);
        res(j) = norm(a - g(j)*b)/norm(a);
    end
    [~, jb] = min(res);
    r = struct('gain', g(jb), 'lag', lags(jb), 'residual', res(jb), 'k_on', k_on, 'skip', skip, ...
        'table', table(lags.', g.', mag2db(abs(g.')), 100*res.', ...
        'VariableNames', {'lag_samples', 'gain', 'gain_dB', 'residual_pct'}), ...
        't', t, 'u_rec', u_rec, 'u_th', u_th);
end
