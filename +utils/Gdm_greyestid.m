function [Gdm,Gdmnom,sys_est,dm_bound,m0,m_est] = Gdm_greyestid(GFRF,freqmin,freqmax,options)
% Estimation of the uncertain 2nd order system on a frequency band limited
% of the non parametric multimodel part of the uncertain pole
% 
arguments
    GFRF 
    freqmin 
    freqmax 
    options.phasesgn = 1 
end
    Ts=GFRF.Ts;
    freq=GFRF.Frequency;
    freq=freq((freq<=2*pi*freqmax)&(freq>=2*pi*freqmin));
    
    GFRF_crop=interp(GFRF,freq);
    
    n=length(GFRF);

    % initial guess from the frequency range and the gain of the system
    gain=1;
    damp=1;
    m=1/db2mag(gain);
    w=2*pi*(freqmax+freqmin)/2;
    c=2*damp*w*m;
    k=w^2*m;
    dm=ones(1,n)*m;
    
    % Greyestimator
    Parameters = {'masses',dm;'damping',c;'rigidity',k;};
    sys0 = idgrey(@Gdm, Parameters, 'd', {options.phasesgn}, Ts);
    opt = greyestOptions('Display','on');
    opt.EnforceStability=true;
    compare(sys0,GFRF_crop)
    sys_est = greyest(GFRF_crop, sys0, opt);
    
    % optimal parameter and estimated model
    estim_param=sys_est.Report.Parameters.ParVector;
    m_est=estim_param(1:n)';
    c_est = estim_param(1+n);
    k_est = estim_param(2+n);
    
    s=tf('s');
    Gdm=ones(size(m_est))*s;
    for i=1:length(m_est)
        Gdm(i)=s^2/((m_est(i))*s^2+c_est*s+k_est);
    end

    m0=mean(m_est);
    dm_bound=max(m_est)-m0;
    Gdm=-options.phasesgn*c2d(Gdm,Ts);
    Gdmnom=s^2/(m0*s^2+c_est*s+k_est);
    Gdmnom=-options.phasesgn*c2d(Gdmnom,Ts);
    compare(Gdm,GFRF_crop)
end