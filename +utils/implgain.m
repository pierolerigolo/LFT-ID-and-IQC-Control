function r = implgain(path, K, lags)
%IMPLGAIN Gain (and delay) between the designed controller and the one that ran
%
%   r = utils.implgain(path, K)          lags -3..3 samples
%   r = utils.implgain(path, K, lags)
%
%   Simulates the designed controller on the recorded acceleration of the
%   closed-loop part of the run, u_th = lsim(K, -y), and fits the recorded
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
end
    [y, ~, u] = utils.loadexp(path);
    k_on = find(abs(u) > 1e-3*max(abs(u)), 1);
    idx  = (k_on:numel(y)).';
    t    = (idx - idx(1))*K.Ts;
    u_th = lsim(K, -y(idx), t);
    u_rec = u(idx);

    g = zeros(size(lags));   res = g;
    for j = 1:numel(lags)
        a = u_rec(max(1, 1 + lags(j)):end + min(0, lags(j)));
        b = u_th(max(1, 1 - lags(j)):end + min(0, -lags(j)));
        g(j)   = (b.'*a)/(b.'*b);
        res(j) = norm(a - g(j)*b)/norm(a);
    end
    [~, jb] = min(res);
    r = struct('gain', g(jb), 'lag', lags(jb), 'residual', res(jb), 'k_on', k_on, ...
        'table', table(lags.', g.', mag2db(abs(g.')), 100*res.', ...
        'VariableNames', {'lag_samples', 'gain', 'gain_dB', 'residual_pct'}), ...
        't', t, 'u_rec', u_rec, 'u_th', u_th);
end
