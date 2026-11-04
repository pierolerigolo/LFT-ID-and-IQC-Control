function [Gdm, Gdmnom, sys_est, dm_bound, m0, m_est] = Gdm_greyestid(GFRF, freqmin, freqmax, options)
%GDM_GREYESTID Grey-box fit of the mass-dependent 2nd-order part G_dm
%
%   [Gdm, Gdmnom, sys_est, dm_bound, m0, m_est] = utils.Gdm_greyestid(GFRF, fmin, fmax)
%
%   Fits, on the band [fmin, fmax] Hz of the multimodel FRF GFRF (1 x n frd,
%   one channel per operating point), the models
%       G_i(s) = s^2 / (m_i s^2 + c s + k)        (Tustin-discretised)
%   with one mass m_i per operating point and common c, k (grey-box model
%   Gdm.m).
%
%   Outputs
%     Gdm      1 x n discrete tf, one model per operating point
%     Gdmnom   nominal model, mass m0
%     sys_est  idgrey estimate (parameters + covariance)
%     m0       centre of the identified mass interval
%     dm_bound half-width:  m_i in [m0 - dm_bound, m0 + dm_bound]
%     m_est    identified masses m_i
%
arguments
    GFRF
    freqmin
    freqmax
    options.phasesgn = 1
end
    Ts   = GFRF.Ts;
    freq = GFRF.Frequency;
    freq = freq((freq <= 2*pi*freqmax) & (freq >= 2*pi*freqmin));
    GFRF_crop = interp(GFRF, freq);
    n = length(GFRF);

    % Initial guess: 1 dB gain, damping ratio 1, resonance in the middle of the band
    m = 1/db2mag(1);
    w = 2*pi*(freqmax + freqmin)/2;
    c = 2*w*m;
    k = w^2*m;

    % Grey-box estimation
    Parameters = {'masses', ones(1, n)*m; 'damping', c; 'rigidity', k};
    sys0 = idgrey(@Gdm, Parameters, 'd', {options.phasesgn}, Ts);
    opt  = greyestOptions('Display', 'on');
    opt.EnforceStability = true;
    sys_est = greyest(GFRF_crop, sys0, opt);

    estim_param = sys_est.Report.Parameters.ParVector;
    m_est = estim_param(1:n)';
    c_est = estim_param(1+n);
    k_est = estim_param(2+n);

    % Centre/half-width of the identified interval (NOT the mean: with the
    % mean, m0 +/- dm_bound misses the extreme mass on the far side whenever
    % the masses are not symmetric about their mean).
    m0       = (max(m_est) + min(m_est))/2;
    dm_bound = (max(m_est) - min(m_est))/2;

    % Tustin is a linear fractional map s -> z, so Gdmnom/(1 + Gdmnom*dm)
    % (the LFT built in b_LFT2ndorderID) is EXACTLY the discretised model of
    % mass m0 - dm. Must match the discretisation in Gdm.m.
    s   = tf('s');
    Gdm = ones(size(m_est))*s;
    for i = 1:n
        Gdm(i) = s^2/(m_est(i)*s^2 + c_est*s + k_est);
    end
    Gdm    = -options.phasesgn*c2d(Gdm, Ts, 'tustin');
    Gdmnom = -options.phasesgn*c2d(s^2/(m0*s^2 + c_est*s + k_est), Ts, 'tustin');

    figure; compare(Gdm, GFRF_crop)
end
