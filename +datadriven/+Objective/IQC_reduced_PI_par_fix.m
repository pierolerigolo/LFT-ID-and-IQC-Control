function [constraints, slack, D, E, PI] = IQC_reduced_PI_par_fix(plant, controller, eta, alpha, beta, omegas, opts)
    %IQC_REDUCED_PI Slacked, reduced IQC constraint with an OPTIMISED
    %(full-block, real-uncertainty) multiplier
    %
    %       PI = [ eta^2*D ,  E  ;
    %                E^T   , -D  ]
    %
    %   with D = D^T >= 0, E = -E^T, and 'eta' a FIXED scalar bound on the
    %   uncertainty (||delta|| <= eta). D and E are introduced here as
    %   YALMIP decision variables, so the multiplier is optimised jointly
    %   with the slack, instead of being supplied as fixed frequency-
    %   domain data (as in the original IQC_reduced_PI).
    %
    %   ------------------------------------------------------------------------
    %   WHY THE CONTROLLER MUST BE FIXED FOR THIS CALL
    %   ------------------------------------------------------------------------
    %   The original LMI is
    %       M = [PI_plus^-1, L ; L', -frakL] >= 0
    %   which is only an LMI (affine in the decision variables) because
    %   PI_plus was a KNOWN CONSTANT — inv() of a decision variable is not
    %   affine. Since PI_plus is now built from D, E (decision variables),
    %   we cannot keep this form. Instead we use the congruence transform
    %   T = blkdiag(PI_plus, I), which (given PI_plus > 0) gives the
    %   equivalent condition
    %       M' = T' * M * T = [PI_plus, PI_plus*L ; L'*PI_plus, -frakL] >= 0.
    %   M' is affine in D, E *provided L and frakL are numeric*, i.e. the
    %   controller (X, Y) is evaluated at FIXED values for this call. If
    %   the controller still contains free sdpvars, the PI_plus*L term is
    %   bilinear in (D,E,X,Y) and this is no longer a solvable SDP.
    %
    %   This function is therefore meant to be used as the "optimise the
    %   multiplier" half of an alternating scheme, mirrored against your
    %   existing controller-synthesis stages:
    %       1) fix D,E  -> optimise controller (X,Y)  [IQC_reduced_PI, original, with numeric PI]
    %       2) fix (X,Y) -> optimise D,E               [THIS function]
    %       repeat until convergence (e.g. objective / eta stops improving).
    %
    %   A warning is raised if the controller frequency response still
    %   contains free decision variables.
    %
    %   ------------------------------------------------------------------------
    %   [constraints, slack, D, E, PI] = IQC_reduced_PI(plant, controller,
    %   eta, alpha, beta, omegas) returns the LMI constraints, the slack
    %   variable, the optimised multiplier blocks D, E, and the assembled
    %   PI = [eta^2*D, E; E', -D] (constant across frequency).
    %
    %   Name-value options mirror the original function:
    %     "epsNormalisation", "epsSVD", "addSlack" - as before.
    %     "epsD"       - D >= epsD*I (avoids D, hence PI_plus, singularity). Default 1e-4.
    %     "epsPIplus"  - PI_plus >= epsPIplus*I, needed for the congruence
    %                    transform above to be a valid equivalence. Default 1e-6.
    %     "warmstartD", "warmstartE" - optional numeric warm starts for D, E.
    %
    %   ------------------------------------------------------------------------
    %   References
    %     [1] V. Gupta, E. Klauser, and A. Karimi, "Non-parametric IQC
    %     multipliers in data-driven robust controller synthesis," 2024,
    %     [Preprint] Available: http://infoscience.epfl.ch/record/307954.
    %     [2] Full-block S-procedure multiplier structure for norm-bounded
    %     real uncertainty ||delta|| <= eta: PI = [eta^2*D, E; E', -D],
    %     D = D' >= 0, E = -E' (cf. Scherer, "LPV control and full-block
    %     multipliers").
    %
    %   Copyright 2024 Vaibhav Gupta, DDMAC, EPFL (MIT License)
    %   Adapted 2026 to optimise the IQC multiplier itself.
    %

    arguments
        plant lti
        controller datadriven.Controller.Controller
        eta (1, 1) double
        alpha (1, 1) double
        beta (1, 1) double
        omegas (:, 1) double
        opts.epsNormalisation (1, 1) double = 1e-6;
        opts.epsSVD (1, 1) double = 1e-6;
        opts.addSlack (1, 1) logical = true;
        opts.epsD (1, 1) double = 1e-4;
        opts.epsPIplus (1, 1) double = 1e-6;
        opts.warmstartD =[]
        opts.warmstartE =[]
        opts.warmstartPI double = [];   % previous full PI = [eta^2*D,E;E',-D] (numeric, 2n x 2n)
        opts.verbose (1, 1) logical = false;
    end

    % Ensure that frequencies are ordered!
    omegas = sort(omegas);
    nomegas = length(omegas);
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

    % The [eta^2*D, E; E', -D] structure requires a square delta-channel
    if np ~= nq
        error("IQC_reduced_PI:dimensionMismatch", ...
            "The [eta^2*D, E; E', -D] structure requires np == nq (got np=%d, nq=%d).", np, nq);
    end
    n = np;


    % ----- Optimised multiplier: D = D^T >= 0, E = -E^T, eta fixed -----
                                 % (avoids sdpvar(...,'skew'), whose assign()
                                 % support is flaky in some YALMIP versions)

    % Reuse the last PI (from a previous multiplier solve, or from
    % whatever fixed PI was last used in the controller stage) as a warm
    % start, by splitting it back into D, E. Re-symmetrise/-skew defensively
    % in case of numerical noise from the previous solve.
    D_ws = [];
    E_ws = [];
    if ~isempty(opts.warmstartPI)
        PI_ws = opts.warmstartPI;
        if size(PI_ws, 1) ~= 2*n || size(PI_ws, 2) ~= 2*n
            error("IQC_reduced_PI:warmstartPISize", ...
                "warmstartPI must be %dx%d to match [eta^2*D, E; E', -D].", 2*n, 2*n);
        end
        D_from_PI = PI_ws(1:n, 1:n) / eta^2;
        E_from_PI = PI_ws(1:n, n+1:end);
        D_ws = (D_from_PI + D_from_PI') / 2;
        E_ws = (E_from_PI - E_from_PI') / 2;
    end
    if ~isempty(opts.warmstartD)
        D_ws = opts.warmstartD;   % explicit D takes precedence
    end
    if ~isempty(opts.warmstartE)
        E_ws = opts.warmstartE;   % explicit E takes precedence
    end
    % IMPORTANT: an sdpvar that is never assign()-ed evaluates to NaN under
    % double(), not 0. Since M below is built from D, E and later passed
    % through double(M) for the SVD-based conditioning step, an
    % unassigned D/E silently turns the whole M into NaN and svd() fails.
    % So always assign a default, even with no warm start.

    if isempty(D_ws)
        D_ws = cell(nomegas, 1);
        for i=1:nomegas
            D_ws{i} = eye(n);
        end
    end
    if isempty(E_ws)
        E_ws = cell(nomegas, 1);
        for i=1:nomegas
            E_ws{i} = 1e-6*eye(n);
        end
    end

    D = cell(nomegas, 1);
    E = cell(nomegas, 1);
    for i = 1:nomegas
        D{i} = sdpvar(n, n, 'hermitian', 'complex');
        E{i} = sdpvar(n, n, 'full', 'complex');

        warmstart(D{i}, D_ws{i});
        warmstart(E{i}, E_ws{i});
    end

    if opts.verbose
        fprintf('D warm-start:\n'); disp(double(D{1}));
        fprintf('E warm-start:\n'); disp(double(E{1}));
    end
    
    multiplierConstraints = cell(length(omegas), 1);
    PI = cell(length(omegas), 1);
    PI_plus = cell(length(omegas), 1);
    PI_minus_22 = cell(length(omegas), 1);
    for i_W = 1:length(omegas)
        Dk=D{i_W};
        Ek=E{i_W};
        PI_11 = eta^2 * Dk;
        PI_12 = Ek;
        PI_21 = pagectranspose(Ek);
        PI_22 = -Dk;
        PI{i_W} = [PI_11, PI_12; PI_21, PI_22];   % constant (frequency-independent) multiplier
    
        PI_plus{i_W} = [ ...
            PI_11 + alpha * eye(np) , PI_12                       ;
            PI_21                   , (beta + eps(1.0)) * eye(nq) ;
            ];
        PI_minus_22{i_W} = PI_22 - beta * eye(nq);
    
        % Dk >= 0 (strictly, to avoid PI_plus singularity - same convention as
        % the alpha_PI/beta_PI >= ~1e-4 regularisation used elsewhere).
        % PI_plus >= eps*I is required for the congruence transform below to
        % be a valid equivalence (it plays the role that "PI_plus was known
        % and could be inverted" played in the original function).
        multiplierConstraints{i_W} = [ ...
            Dk >= opts.epsD * eye(n) ; ...
            PI_plus{i_W} >= opts.epsPIplus * eye(np + nq) ; ...
            Ek + Ek' == zeros(n) ; ...   % enforces Ek = -Ek'
            ];
    end

    % ----- SDPVAR definitions -----
    if opts.addSlack
        slack = sdpvar(1);
        warmstart(slack, 1);
    else
        slack = 0;
    end

    % ----- Frequency responses (plant only - PI is constant in frequency) -----
    Gf = freqresp(plant, omegas);

    % ===== Constraints =====
    constraints_tmp = cell(length(omegas), 1);

    for i_W = 1:length(omegas)

        % Prepare different values
        [X, Y] = freqresp(controller, omegas(i_W));

        if i_W == 1 && (isa(X, 'sdpvar') || isa(Y, 'sdpvar'))
            warning("IQC_reduced_PI:controllerNotFixed", ...
                ['Controller frequency response still contains free decision ' ...
                 'variables. Jointly optimising PI (D,E) with a free controller ' ...
                 'makes the PI_plus*L block bilinear, which is not a solvable ' ...
                 'SDP. Fix the controller (pass its current numeric values) ' ...
                 'before calling this function in "optimise the multiplier" mode.']);
        end

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
            Phi_c = double(Phi);   % numeric, since controller is fixed for this call

            Psi = eye(nq) - pinv(G21)*G21;

            L  = [...
                G11 * Phi + G12 * X , G11 * Psi ;
                Phi                 , Psi       ;
                ];
            L = double(L);   % controller fixed -> L is a numeric constant

            % frakL is affine in D via PI_minus_22 (tmp); everything else
            % here is numeric because the controller is fixed.
            tmp = PI_minus_22{i_W};
            frakL11 = Phi'*tmp*Phi_c + Phi_c'*tmp*Phi - Phi_c'*tmp*Phi_c;
            frakL12 = Phi'*tmp*Psi - alpha * (G11 * Phi + G12 * X)' * (G11*Psi);
            frakL22 = Psi'*tmp*Psi - alpha * (G11 * Psi)'           * (G11*Psi);
            frakL = [ frakL11, frakL12; frakL12', frakL22];

            % Congruence-transformed LMI (avoids inv(PI_plus), which would
            % be non-affine now that PI_plus depends on D, E):
            %   M' = [PI_plus, PI_plus*L; L'*PI_plus, -frakL] >= 0
            % equivalent to the original M >= 0, given PI_plus > 0.
            PI_plus_iw=PI_plus{i_W};
            M = [PI_plus_iw, PI_plus_iw*L; L'*PI_plus_iw, -frakL];

        elseif controller.factorisation == "Left"
            error("IQC not defined for Left factorisation! Yet!");
        end

        M = M + slack * eye(size(M));
        % Normalisation for numerical stability. double(M) uses the
        % current assignment of D, E (warm start / previous solve), same
        % pattern as the original function used for controller variables.
        if opts.epsNormalisation > 0
            M_num = double(M);
            if any(~isfinite(M_num(:)))
                % warning("IQC_reduced_PI:unassignedVariable")
                     % "double(M) contains NaN/Inf at omega index %d - likely an ", ...
                     % "sdpvar in M (D, E, or slack) has no assign()-ed value. ", ...
                     % "Skipping SVD conditioning for this frequency point.");
                M = (M + M') / 2;
            else
                [s,v,d] = svd(M_num);
                v = max(diag(v), opts.epsNormalisation);
                Lnorm = s * sqrt(diag(1./v)) * d';
                M = ((Lnorm' * M * Lnorm) + (Lnorm' * M * Lnorm)') / 2;
            end
        else
            M = (M + M') / 2;
        end

        if opts.addSlack
            % Stability condition for slacked variation (unaffected by PI
            % change - depends only on the fixed controller).
            Phi = (Y - G22*X);
            Phi_c = double(Phi);
            stability = Phi'*Phi_c + Phi_c'*Phi - Phi_c'*Phi_c;
            if opts.epsNormalisation > 0
                [s,v,d] = svd(double(stability));
                v = max(diag(v), opts.epsNormalisation);
                Lnorm = s * sqrt(diag(1./v)) * d';
                stability = ((Lnorm' * stability * Lnorm) + (Lnorm' * stability * Lnorm)') / 2;
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
        constraints = [multiplierConstraints{:}, constraints_tmp{:}, slack >= 0];
    else
        constraints = [multiplierConstraints{:}, constraints_tmp{:}];
    end

end