function [LMI, OBJ, slack, gamma] = IQC_performance(G, K, PI, alpha, beta, omegas, opts)
    %IQC 
    % 
    %   Hinf only
    % 
    %   ----------------------------------------------------------------------------
    %   Output:
    %       LMI : LMI
    %       obj : OBJ
    % 
    %   ----------------------------------------------------------------------------
    %   Copyright
    % 
    
    arguments
        G lti
        K datadriven.Controller.Controller
        PI lti
        alpha (1, 1) double
        beta (1, 1) double
        omegas (:, 1) double
        opts.epsNormalisation (1, 1) double = 1e-4;   
        opts.epsSVD (1, 1) double = 1e-6;
        opts.beta_performance (1, 1) double = 1e-2;
        opts.last_gamma (1, 1) double = 0;
        opts.addSlack (1, 1) logical = false;
        opts.slack_max (1, 1) double = 1;
        opts.slack_weight (1, 1) double = 1000;   % OBJ = gamma + slack_weight*slack (must dominate gamma)
        opts.pinorm (1,1) logical = false;
    end

    % Ensure that frequencies are ordered!
    omegas = sort(omegas);
    
    nq = length(G.InputGroup.U1);
    nw = length(G.InputGroup.U2);
    % nu = length(G.InputGroup.U3);

    np = length(G.OutputGroup.Y1);
    nz = length(G.OutputGroup.Y2);
    ny = length(G.OutputGroup.Y3);
    
    % Check if base assumption is satisfied
    if K.factorisation == "Right"
        if ny > nq
            error("Base assumption is not satisfied!!!");
        end
    elseif K.factorisation == "Left"
        error("IQC not defined for Left factorisation! Yet!");
    end
     
    % ----- SDPVAR definitions -----
    gamma = sdpvar(1);
    gammawarm=norm(lft(G(np+1:end, np+1:end), ss(K)),'inf')^2 * 1e2;
    warmstart(gamma, gammawarm);
    fprintf('gamma warm %.6f\n',gammawarm)
    if opts.addSlack
        slack = sdpvar(1);
        warmstart(slack, .1);
    else
        slack = 0;
    end


    PI_11 = PI(1:np, 1:np);
    PI_12 = PI(1:np, np+1:end);
    PI_21 = PI(np+1:end, 1:np);
    PI_22 = PI(np+1:end, np+1:end);
    

    PI_plus = [ ...
        PI_11 + alpha * eye(np) , PI_12         ;
        PI_21                   , beta * eye(nq);
        ];
    
        % ----- PI factorisation -----
    PI_scaler = norm(PI_plus, inf);  % REMOVED: PI already normalised in controllersyn

    fprintf('Pi scaler: %.3e ',PI_scaler)
    % PI_scaler = 1e13;
    if opts.pinorm
        fprintf('used\n')
    else
        PI_scaler=1;
        fprintf('not used\n')
    end
    PI_plus=PI_plus/PI_scaler;
    PI_minus_22 = (PI_22 - beta * eye(nq))/PI_scaler;
    
    % show the range of Pi plus and Pi minus
    fprintf('range of Pi plus\nmin:')
    disp(min(min(min(PI_plus.ResponseData))))
    fprintf('max:')
    disp(max(max(max(PI_plus.ResponseData))))
    fprintf('range of Pi minus\nmin:')
    disp(min(min(min(PI_minus_22.ResponseData))))
    fprintf('max:')
    disp(max(max(max(PI_minus_22.ResponseData))))
    % ----- Frequency responses -----
    Gf = freqresp(G, omegas);
    PI_plus_f = freqresp(PI_plus, omegas);
    PI_minus_22_f = freqresp(PI_minus_22, omegas);
    
    % ===== Constraints =====    
    OBJ = gamma;
    LMI_collect = cell(length(omegas), 1);

    p_index = [G.OutputGroup.Y1, G.OutputGroup.Y2];
    q_index = [G.InputGroup.U1 , G.InputGroup.U2 ];

    for i_W = 1:length(omegas)
        % Prepare different values
        [X, Y] = freqresp(K, omegas(i_W));
    
        G11 = Gf(p_index, q_index, i_W);
        G12 = Gf(p_index, G.InputGroup.U3, i_W);
        G21 = Gf(G.OutputGroup.Y3, q_index, i_W);
        G22 = Gf(G.OutputGroup.Y3, G.InputGroup.U3, i_W);

        PI_plus_inv = inv(PI_plus_f(:, :, i_W));

        PIrp_plus_inv = [ ...
            blkdiag(PI_plus_inv(1:np, 1:np), gamma * eye(nz)), ...
                blkdiag(PI_plus_inv(1:np, np+1:end), zeros(nz, nw));
            blkdiag(PI_plus_inv(np+1:end, 1:np), zeros(nw, nz)), ...
                blkdiag(PI_plus_inv(np+1:end, np+1:end), (opts.beta_performance)^-1 * eye(nw));
            ];

        tmp = blkdiag(PI_minus_22_f(:, :, i_W), -(1 + opts.beta_performance) * eye(nw));   

        if K.factorisation == "Right"
            % Normalisation of plant
            if (opts.epsSVD >= 0)
                s = svds(G21, 1, "smallest"); % Smallest singular value
                G21 = G21 / s;
                G12 = G12 * s;
            end
            Phi = pinv(G21) * (Y - G22*X);
            Phi_c = double(Phi);
            
            Psi = eye(nq+nw) - pinv(G21)*G21;

            L  = [...
                G11 * Phi + G12 * X , G11 * Psi ;
                Phi                 , Psi       ;
                ];

            frakL11 = Phi'*tmp*Phi_c + Phi_c'*tmp*Phi - Phi_c'*tmp*Phi_c;
            frakL12 = Phi'*tmp*Psi - alpha * (G11 * Phi + G12 * X)' * (G11*Psi);
            frakL22 = Psi'*tmp*Psi - alpha * (G11 * Psi)'           * (G11*Psi);
            frakL = [ frakL11, frakL12; frakL12', frakL22];

            M = [PIrp_plus_inv, L; L', -frakL];
            
        elseif K.factorisation == "Left"
            error("IQC not defined for Left factorisation! Yet!");
        end
       
        % Normalisation for numerical stability
        if opts.epsNormalisation > 0
            [s,v,d] = svd(double(M));
            v = max(diag(v), opts.epsNormalisation);
            L = s * sqrt(diag(1./v)) * d';
            M = ((L' * M * L) + (L' * M * L)') / 2;
        else
            M = (M + M') / 2;
        end

        % Add slack to diagonal (mirrors IQC_reduced behaviour)
        M = M + slack * eye(size(M));

        Phi = (Y - G22*X);
        Phi_c = double(Phi);
        stability = Phi'*Phi_c + Phi_c'*Phi - Phi_c'*Phi_c;
        if opts.epsNormalisation > 0
            [s,v,d] = svd(double(stability));
            v = max(diag(v), opts.epsNormalisation);
            L = s * sqrt(diag(1./v)) * d';
            stability = ((L' * stability * L) + (L' * stability * L)') / 2;
        else
            stability = (stability + stability') / 2;
        end

        % Stability constraint: hard or slacked
        if opts.addSlack
            LMI_collect{i_W} = [
                M >= 0;
                stability >= 1e-6;
                ];
        else
            LMI_collect{i_W} = [
                M >= 0;
                stability >= 1e-6;
                ];
        end
    end
    if opts.addSlack
        LMI = [LMI_collect{:}, gamma >= 0, slack >= 0, slack <= opts.slack_max];
        OBJ = gamma + opts.slack_weight*slack;
    else
        LMI = [LMI_collect{:}, gamma >= 0];
        OBJ = gamma;
    end

end