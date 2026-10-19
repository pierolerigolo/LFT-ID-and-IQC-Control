function [constraints, slack, D, E, PI] = IQC_reduced_PI_par_fix(plant, controller, eta, alpha, beta, omegas, opts)
    %IQC_REDUCED_PI Slacked, reduced IQC constraint with an OPTIMISED
    %(full-block, real-uncertainty) multiplier
    %
    %       PI(w) = [ eta^2*D(w) ,  E(w)  ;
    %                   E(w)'    , -D(w)  ]
    %
    %   with D = D' >= 0, E = -E', and 'eta' a FIXED scalar bound on the
    %   uncertainty. D and E are YALMIP decision variables, one (n x n)
    %   matrix per frequency, returned as CELL ARRAYS (nomegas x 1).
    %
    %   The controller MUST be fixed (numeric) for this call: the congruence
    %   transformed LMI
    %       M' = [PI_plus, PI_plus*L; L'*PI_plus, -frakL] >= 0
    %   is affine in (D,E) only if L and frakL's controller parts are numeric.
    %
    %   ------------------------------------------------------------------------
    %   HOW TO USE THE SLACK TO GET THE "BEST" MULTIPLIER
    %   ------------------------------------------------------------------------
    %   With "addSlack" true the LMI is  M + slack*I >= 0.
    %     - "slackLower" = 0 (default): slack >= 0, so minimising slack only
    %       tells you HOW MUCH the constraint is violated (0 = feasible).
    %     - "slackLower" < 0 (e.g. -1): slack may become NEGATIVE; minimising
    %       it returns the multiplier that satisfies the constraint with the
    %       most room to spare (slack* < 0 = room, slack* > 0 = infeasible).
    %   Because the feasible set is a cone in (D,E), "DMax" bounds D <= DMax*I
    %   so the minimisation cannot go to -infinity by just scaling D, E up.
    %
    %   Name-value options:
    %     "epsNormalisation", "epsSVD", "addSlack" - as before.
    %     "epsD"       - D >= epsD*I. Default 1e-4.
    %     "DMax"       - D <= DMax*I. Default Inf (no bound). Use 1 together
    %                    with slackLower < 0 (D = 1, E = 0 is then the fixed
    %                    scalar multiplier, which is inside the feasible set).
    %     "slackLower" - slack >= slackLower. Default 0.
    %     "epsPIplus"  - PI_plus >= epsPIplus*I. Default 1e-6.
    %     "warmstartD", "warmstartE" - numeric n x n, or cell {nomegas x 1}.
    %
    %   Note: the "stability" LMI of IQC_reduced is NOT added here - with a
    %   fixed controller it does not depend on D, E (it is a numeric constant).
    %
    %   ------------------------------------------------------------------------
    %   References
    %     [1] V. Gupta, E. Klauser, and A. Karimi, "Non-parametric IQC
    %     multipliers in data-driven robust controller synthesis," 2024.
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
        opts.epsD (1, 1) double = 1e-2;
        opts.DMax (1, 1) double = Inf;
        opts.slackLower (1, 1) double = 0;
        opts.epsPIplus (1, 1) double = 1e-4;
        opts.warmstartD = []
        opts.warmstartE = []
        opts.warmstartPI double = [];   % previous full PI (numeric, 2n x 2n)
        opts.verbose (1, 1) logical = false;
    end

    % Ensure that frequencies are ordered!
    omegas = sort(omegas);
    nomegas = length(omegas);
    nq = length(plant.InputGroup.U1);
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

    % ----- Warm starts (always assigned: an unassigned sdpvar is NaN under double()) -----
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
    if ~isempty(opts.warmstartD), D_ws = opts.warmstartD; end
    if ~isempty(opts.warmstartE), E_ws = opts.warmstartE; end

    % Make them cells {nomegas x 1}
    if isempty(D_ws)
        D_ws = repmat({eye(n)}, nomegas, 1);
    elseif ~iscell(D_ws)
        D_ws = repmat({D_ws}, nomegas, 1);
    end
    if isempty(E_ws)
        E_ws = repmat({1e-6*eye(n)}, nomegas, 1);
    elseif ~iscell(E_ws)
        E_ws = repmat({E_ws}, nomegas, 1);
    end

    % ----- Optimised multiplier: one D_k, E_k per frequency -----
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

    multiplierConstraints = cell(nomegas, 1);
    PI = cell(nomegas, 1);
    PI_plus = cell(nomegas, 1);
    PI_minus_22 = cell(nomegas, 1);
    for i_W = 1:nomegas
        Dk = D{i_W};
        Ek = E{i_W};
        PI_11 = eta^2 * Dk;
        PI_12 = Ek;
        PI_21 = Ek';
        PI_22 = -Dk;
        PI{i_W} = [PI_11, PI_12; PI_21, PI_22];

        PI_plus{i_W} = [ ...
            PI_11 + alpha * eye(np) , PI_12                       ;
            PI_21                   , (beta + eps(1.0)) * eye(nq) ;
            ];
        PI_minus_22{i_W} = PI_22 - beta * eye(nq);

        multiplierConstraints{i_W} = [ ...
            Dk >= opts.epsD * eye(n) ; ...
            PI_plus{i_W} >= opts.epsPIplus * eye(np + nq) ; ...
            Ek + Ek' == zeros(n) ; ...   % enforces Ek = -Ek'
            ];
        if isfinite(opts.DMax)
            multiplierConstraints{i_W} = [multiplierConstraints{i_W}; ...
                Dk <= opts.DMax * eye(n)];
        end
    end

    % ----- SDPVAR definitions -----
    if opts.addSlack
        slack = sdpvar(1);
        warmstart(slack, 1);
    else
        slack = 0;
    end

    % ----- Frequency responses -----
    Gf = freqresp(plant, omegas);

    % ===== Constraints =====
    constraints_tmp = cell(nomegas, 1);

    for i_W = 1:nomegas

        % Prepare different values
        [X, Y] = freqresp(controller, omegas(i_W));

        if i_W == 1 && (isa(X, 'sdpvar') || isa(Y, 'sdpvar'))
            warning("IQC_reduced_PI:controllerNotFixed", ...
                ['Controller frequency response still contains free decision ' ...
                 'variables. Fix the controller before calling this function ' ...
                 'in "optimise the multiplier" mode.']);
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

            tmp = PI_minus_22{i_W};
            frakL11 = Phi'*tmp*Phi_c + Phi_c'*tmp*Phi - Phi_c'*tmp*Phi_c;
            frakL12 = Phi'*tmp*Psi - alpha * (G11 * Phi + G12 * X)' * (G11*Psi);
            frakL22 = Psi'*tmp*Psi - alpha * (G11 * Psi)'           * (G11*Psi);
            frakL = [ frakL11, frakL12; frakL12', frakL22];

            % Congruence-transformed LMI (equivalent to M >= 0 given PI_plus > 0)
            PI_plus_iw = PI_plus{i_W};
            M = [PI_plus_iw, PI_plus_iw*L; L'*PI_plus_iw, -frakL];

        elseif controller.factorisation == "Left"
            error("IQC not defined for Left factorisation! Yet!");
        end

        M = M + slack * eye(size(M));
        % Conditioning by a constant congruence (built from the warm-start
        % values). A congruence preserves "M + slack*I >= 0" exactly, so the
        % meaning of the slack is unchanged.
        if opts.epsNormalisation > 0
            M_num = double(M);
            if any(~isfinite(M_num(:)))
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

        constraints_tmp{i_W} = M >= 0;
    end

    if opts.addSlack
        constraints = [multiplierConstraints{:}, constraints_tmp{:}, slack >= opts.slackLower];
    else
        constraints = [multiplierConstraints{:}, constraints_tmp{:}];
    end

end