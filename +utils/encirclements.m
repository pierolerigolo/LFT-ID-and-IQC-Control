function N = encirclements(K, GFRF, n_interp)
%ENCIRCLEMENTS Net encirclements of -1 by L = G K on every measured FRF
%
%   N = utils.encirclements(K, GFRF)          (n_interp = 5000)
%   N = utils.encirclements(K, GFRF, n_interp)
%
%   N(i) = winding of 1 + K*GFRF(i) around 0 over -w_max..w_max (counter-
%   clockwise > 0). With K and the plant stable, the closed loop on FRF i is
%   stable iff N(i) = 0. Same computation as utils.nyquistcheck: the FRFs are
%   interpolated on n_interp linear points plus dense points around the
%   lightly damped poles of K (zeta < 0.05), between their first and last
%   measured frequency (L assumed small below the first one).
%
arguments
    K
    GFRF
    n_interp (1,1) double = 5000
end
    pK = pole(K);
    if K.Ts > 0, sK = log(pK)/K.Ts; else, sK = pK; end
    wn = abs(sK);   zK = -real(sK)./max(wn, eps);
    w_K = [];
    for q = find(imag(sK) > 0 & zK < 0.05).'
        w_K = [w_K; wn(q)*linspace(1 - 30*max(zK(q), 1e-4), 1 + 30*max(zK(q), 1e-4), 2000).']; %#ok<AGROW>
    end
    nOP = size(GFRF, 2);
    N = zeros(nOP, 1);
    for i = 1:nOP
        Gi = GFRF(i);
        w0 = min(Gi.Frequency);   w1 = max(Gi.Frequency);
        w  = unique([linspace(w0, w1, n_interp).'; w_K(w_K >= w0 & w_K <= w1)]);
        L  = squeeze(freqresp(K, w)).*squeeze(interp(Gi, w).ResponseData);
        L  = L(isfinite(L));
        ph = unwrap(angle(1 + [conj(flipud(L)); L]));
        N(i) = round((ph(end) - ph(1))/(2*pi));
    end
end
