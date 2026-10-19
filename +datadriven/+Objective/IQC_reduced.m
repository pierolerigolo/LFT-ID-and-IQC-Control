  function [constraints, slack] = IQC_reduced(plant, controller, PI, alpha, beta, omegas, opts)
    %IQC Sets up slaked, reduced IQC constraint for a LFR plant and given IQC multiplier
    %   
    %   This function sets up the slacked IQC constraint for the generalised
    %   plant represented in Linear Fraction Representation (LFR) form. The
    %   function returns the LMI contraints and the objective as specified in
    %   Gupta et. al. (2024).
    % 
    %   [constraints, slack] = IQC_reduced(plant, controller, PI, alpha, beta, omegas)
    %   gives the constraint and slack for the IQC defined by the IQC multiplier
    %   "PI" which can be split into positive and negative definite matrices
    %   using 'alpha' and 'beta' as specified in Corollary 1 of Gupta et. al.
    %   (2024).
    % 
    %   [constraints, slack] = IQC_reduced(..., "epsNormalisation", eps) sets
    %   the eps for normalisation of the constraint. Negative eps disables the
    %   normalisation.
    % 
    %   [constraints, slack] = IQC_reduced(..., "epsSVD", eps) sets the eps for
    %   normalisation of the plant. Negative eps disables the normalisation.
    % 
    %   [constraints, ~] = IQC_reduced(..., "addSlack", false) returns the
    %   strict version of the constraints.
    % 
    %   ------------------------------------------------------------------------
    %   Outputs
    %     * constraints : List of LMI constraints in YALMIP framework
    %     * slack       : Slack in YALMIP framework
    % 
    %   ------------------------------------------------------------------------
    %   References
    %     [1] V. Gupta, E. Klauser, and A. Karimi, “Non-parametric IQC
    %     multipliers in data-driven robust controller synthesis," 2024, 
    %     [Preprint] Available: http://infoscience.epfl.ch/record/307954.
    % 
    %   ------------------------------------------------------------------------
    %   Copyright 2024 Vaibhav Gupta, DDMAC, EPFL (MIT License)
    %
    
    arguments
        plant lti
        controller datadriven.Controller.Controller
        PI lti
        alpha (1, 1) double
        beta (1, 1) double
        omegas (:, 1) double
        opts.epsNormalisation (1, 1) double = 1e-6;
        opts.epsSVD (1, 1) double = 1e-6;
        opts.addSlack (1, 1) logical = true;
        opts.pinorm (1,1) logical = false;
    end

    % Ensure that frequencies are ordered!
    omegas = sort(omegas);
    nq = length(plant.InputGroup.U1);
    nu = length(plant.InputGroup.U2);
    np = length(plant.OutputGroup.Y1);
    ny = length(plant.OutputGroup.Y2);
    
    % Check if base assumption is satisfied
    if controller.factorisation == "Right"
        if ny > np
            error("Base assumption is not satisfied!!!");
        end
    elseif controller.factorisation == "Left"
        error("IQC not defined for Left factorisation! Yet!");
    end

    PI_11 = PI(1:np, 1:np);
    PI_12 = PI(1:np, np+1:end);
    PI_21 = PI(np+1:end, 1:np);
    PI_22 = PI(np+1:end, np+1:end);

    PI_plus = [ ...
        PI_11 + alpha * eye(np) , PI_12         ;
        PI_21                   , (beta + eps(1.0)) * eye(nq) ;
        ];
        % ----- PI factorisation -----
    PI_scaler = norm(PI_plus, inf);  

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
    % ----- SDPVAR definitions -----
    if opts.addSlack
        slack = sdpvar(1);
        warmstart(slack, 1);
    else
        slack = 0;
    end
    
    % ----- Frequency responses -----
    Gf = freqresp(plant, omegas);
    PI_plus_f = freqresp(PI_plus, omegas);
    PI_minus_22_f = freqresp(PI_minus_22, omegas);
    
    % ===== Constraints =====    
    constraints_tmp = cell(length(omegas), 1);

    for i_W = 1:length(omegas)
        % Prepare different values
        [X, Y] = freqresp(controller, omegas(i_W));
    
        G11 = Gf(plant.OutputGroup.Y1, plant.InputGroup.U1, i_W);
        G12 = Gf(plant.OutputGroup.Y1, plant.InputGroup.U2, i_W);
        G21 = Gf(plant.OutputGroup.Y2, plant.InputGroup.U1, i_W);
        G22 = Gf(plant.OutputGroup.Y2, plant.InputGroup.U2, i_W);

        if controller.factorisation == "Right"
            % Normalisation of plant
            if (opts.epsSVD >= 0)
                s = svds(G21, 1, "smallest"); % Smallest singular value
                G21 = G21 / s;
                G12 = G12 * s;
            end
    
            Phi = pinv(G21) * (Y - G22*X);
            Phi_c = double(Phi);
            
            Psi = eye(nq) - pinv(G21)*G21;

            L  = [...
                G11 * Phi + G12 * X , G11 * Psi ;
                Phi                 , Psi       ;
                ];

            PI_plus_inv = inv(PI_plus_f(:, :, i_W));

            tmp = PI_minus_22_f(:, :, i_W);
            frakL11 = Phi'*tmp*Phi_c + Phi_c'*tmp*Phi - Phi_c'*tmp*Phi_c;
            frakL12 = Phi'*tmp*Psi - alpha * (G11 * Phi + G12 * X)' * (G11*Psi);
            frakL22 = Psi'*tmp*Psi - alpha * (G11 * Psi)'           * (G11*Psi);
            frakL = [ frakL11, frakL12; frakL12', frakL22];

            M = [PI_plus_inv, L; L', -frakL];
            
        elseif controller.factorisation == "Left"
            error("IQC not defined for Left factorisation! Yet!");
        end
       
        M = M + slack * eye(size(M));
        % Normalisation for numerical stability
        if opts.epsNormalisation > 0
            [s,v,d] = svd(double(M));
            v = max(diag(v), opts.epsNormalisation);
            L = s * sqrt(diag(1./v)) * d';
            M = ((L' * M * L) + (L' * M * L)') / 2;
        else
            M = (M + M') / 2;
        end

        if opts.addSlack        
            % Stability condition for slacked variation
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

            constraints_tmp{i_W} = [
                M >= 0;
                stability >= 0;  % Ensure stability for non-zero gamma
                ];
        else
            constraints_tmp{i_W} = M >= 0;
        end
    end

    if opts.addSlack
        constraints = [constraints_tmp{:}, slack >= 0];
    else
        constraints = [constraints_tmp{:}];
    end

end