function [Gest,GLTI,sys_est] = GLTI_greyestid(GFRF,orderGLTI,Gdm)
%GLTI_GREYESTID Summary of this function goes here
%   Detailed explanation goes here
    Ts=GFRF.Ts;

    b  = [1,zeros(1,orderGLTI)];                  % 1 x (n+1)
    a=  zeros(1,orderGLTI);                    % 1 x n, the leading 1 is added in the function
    
    Parameters = {'numerator',b;'denominator',a};
    sys0 = idgrey('GdmGLTI', Parameters, 'd',{Gdm},0.001);
    compare(sys0,GFRF)
    
    % fixer le paramtre qui ne change pas 
    opt = greyestOptions('Display','on');
    opt.EnforceStability=true;
    sys_est = greyest(GFRF, sys0, opt);
    
    b_est=sys_est.Report.Parameters.ParVector(1:(1+orderGLTI))';
    a_est=sys_est.Report.Parameters.ParVector((2+orderGLTI):end)';
    
    GLTI=tf(real(b_est),[1,real(a_est)],Ts);
    Gest=Gdm*GLTI;
end

