function [att_dB, snr_dB, seg] = sineatt(path, trans)
%SINEATT Attenuation and SNR of a single-tone closed-loop experiment
%
%   [att_dB, snr_dB, seg] = utils.sineatt(path, trans)
%
%   The recording is: disturbance on with the controller off (open loop),
%   controller switched on (closed loop), disturbance off (noise). The three
%   segments are found from the disturbance command and the controller output
%   (the closed-loop segment ends when the controller or the disturbance stops);
%   trans samples are skipped after each switch (transients), 500 before each
%   segment end, and the three segments get the same length.
%     att_dB = 20 log10( rms(closed loop) / rms(open loop) )
%     snr_dB = 20 log10( rms(closed loop) / rms(noise) )
%     seg    = [start; end] of the open-loop, closed-loop and noise segments
%
arguments
    path
    trans (1,1) double = 3000
end
    [y, d, u] = utils.loadexp(path);

    dist = active(d);
    k = find(dist, 1);   dist_end = find(dist, 1, 'last');
    if isempty(k), error('sineatt:noDist', 'No disturbance found in %s.', path); end
    ctrl = active(u) & dist;
    k_on = find(ctrl, 1);   ctrl_end = find(ctrl, 1, 'last');
    if isempty(k_on), error('sineatt:noCtrl', 'Controller never switched on in %s.', path); end

    starts = [k + trans,   k_on + trans,     dist_end + trans];
    ends   = [k_on - 500,  ctrl_end - 500,   numel(y)];
    L = min(ends - starts + 1);
    if L < 1000
        error('sineatt:short', 'Segments too short (%d samples) in %s: reduce trans.', L, path);
    end
    seg = [starts; starts + L - 1];

    r = arrayfun(@(i) rms(detrend(y(seg(1, i):seg(2, i)))), 1:3);
    att_dB = 20*log10(r(2)/r(1));
    snr_dB = 20*log10(r(2)/r(3));
end

function a = active(x)
% samples where the moving RMS of x exceeds 5 % of its maximum
    m = sqrt(movmean(x.^2, 500));
    a = m > 0.05*max(m);
end
