function [S, Uclfrd, Yclfrd, Yolfrd] = sweepS(path, Ts, pre)
%SWEEPS Experimental sensitivity from an open-loop / closed-loop sine-sweep pair
%
%   S = utils.sweepS(path, Ts, pre)
%   [S, Uclfrd, Yclfrd, Yolfrd] = utils.sweepS(...) also returns the raw FFT
%   bins (frd) of the controller output and of the acceleration in the
%   closed-loop window, and of the acceleration in the open-loop window.
%
%   The recording holds identical sweeps of the disturbance, with the controller
%   off (open loop) and on (closed loop), in any order. Only full-length sweeps
%   are kept (>= 80 % of the longest); the first of each kind is used.
%   Both windows start pre samples before their sweep and have the same length:
%       S(jw) = FFT(y closed loop) / FFT(y open loop)      (frd, 0 < w <= pi/Ts)
%
arguments
    path
    Ts  (1,1) double = 0.001
    pre (1,1) double = 1000
end
    [y, d, u] = utils.loadexp(path);
    % figure;hold on;plot(y);plot(d);plot(u); %debugger
    m = sqrt(movmean(d.^2, 500));
    on = m > 0.05*max(m);
    starts = find(diff([0; on]) == 1);
    ends   = find(diff([on; 0]) == -1);
    full   = (ends - starts) >= 0.8*max(ends - starts);
    starts = starts(full);   ends = ends(full);

    % closed loop = sweep during which the controller output is active
    ru  = arrayfun(@(i) rms(u(starts(i):ends(i))), 1:numel(starts));
    isc = ru > 0.1*max(ru);
    iol = find(~isc, 1);   icl = find(isc, 1);
    if isempty(iol) || isempty(icl)
        error('sweepS:blocks', 'Need one open-loop and one closed-loop sweep in %s.', path);
    end
    starts = starts([iol, icl]);   ends = ends([iol, icl]);
    iol = 1;   icl = 2;

    L  = min(ends - starts + 1) + pre;
    s0 = max(starts - pre, 1);
    L  = min([L, numel(y) - s0(1) + 1, numel(y) - s0(2) + 1]);
    Yol = fft(y(s0(iol):s0(iol) + L - 1));
    Ycl = fft(y(s0(icl):s0(icl) + L - 1));
    Ucl = fft(u(s0(icl):s0(icl) + L - 1));

    k = (1:floor(L/2)).';                 % 0 < f <= Nyquist
    w = 2*pi*k/(L*Ts);
    Uclfrd =frd(reshape(Ucl(k + 1), 1, 1, []), w, Ts);
    Yclfrd = frd(reshape(Ycl(k + 1), 1, 1, []), w, Ts);
    Yolfrd = frd(reshape(Yol(k + 1), 1, 1, []), w, Ts);
    S = frd(reshape(Ycl(k + 1)./Yol(k + 1), 1, 1, []), w, Ts);
end
