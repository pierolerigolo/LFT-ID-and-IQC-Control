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
    %   With "multiplier","constant" a single D, E is shared by all
    %   frequencies (cells then hold the same sdpvar) and the multiplier
    %   constraints (D >= epsD, PI_plus >= epsPIplus, E = -E', D <= DMax)
    %   are added ONCE, exactly like the original constant IQC_reduced_PI.
    %
    %   The controller MUST be fixed (numeric) for this call: the congruence
    %   transformed LMI
    %       M' = [PI_plus, PI_plus*L; L'*PI_plus, -frakL] >= 0
    %   is affine in (D,E) only if L and frakL's controller parts are numeric.
    %
    %   ------------------------------------------------------------------------
    %   PURE FEASIBILITY (addSlack = false) vs SLACK
    %   ------------------------------------------------------------------------
    %   For a square/SISO delta-channel Psi = 0 and the last rows/cols of M
    %   are identically zero. With "dropNullRows" true (default) they are
    %   removed so that M >= 0 has an interior. NOTE: without an objective the
    %   multiplier returned by the solver is NOT unique, and which point you
    %   get depends on the formulation (number of constraints, size of M,
    %   epsD, ...). Set "dropNullRows" false to reproduce the original
    %   4x4 formulation exactly.
    %
    %   With "addSlack" true the LMI is  M + slack*I >= 0.
    %     - "slackLower" = 0 (default): slack >= 0 (violation measure).
    %     - "slackLower" < 0: slack may be negative -> margin.
    %
    %   Name-value options:
    %     "epsNormalisation", "epsSVD", "addSlack" - as before.
    %     "multiplier"   - "frequency" (default) or "constant".
    %     "dropNullRows" - drop structurally zero rows/cols of M. Default true.
    %     "epsD"         - D >= epsD*I. Default 1e-4.
    %     "DMax"         - D <= DMax*I. Default Inf (no bound).
    %     "slackLower"   - slack >= slackLower. Default 0.
    %     "epsPIplus"    - PI_plus >= epsPIplus*I. Default 1e-6.
    %     "warmstartD", "warmstartE" - numeric n x n, or cell {nomegas x 1}.
    %
    %   Note: the "stability" LMI of IQC_reduced is NOT added here - with a
    %   fixed controller it does not depend on D, E (it is a numeric constant).
    %
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
        opts.epsD (1, 1) double = 1e-4;
        opts.DMax (1, 1) double = Inf;
        opts.slackLower (1, 1) double = 0;
        opts.epsPIplus (1, 1) double = 1e-6;
        opts.multiplier (1, 1) string {mustBeMember(opts.multiplier, ["frequency", "constant"])} = "frequency";
        opts.dropNullRows (1, 1) logical = true;
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
        E_ws = repmat({zeros(n)}, nomegas, 1);
    elseif ~iscell(E_ws)
        E_ws = repmat({E_ws}, nomegas, 1);
    end

    % ----- Optimised multiplier: one D_k, E_k per frequency -----
    D = cell(nomegas, 1);
    E = cell(nomegas, 1);
    for i = 1:nomegas
        if opts.multiplier == "constant" && i > 1
            % one D, E shared by all frequencies (the original fixed version)
            D{i} = D{1};
            E{i} = E{1};
            continue
        end
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

        if opts.multiplier == "constant" && i_W > 1
            % Same variables as frequency 1 -> same constraints. Adding them
            % nomegas times would only re-weight the barrier of the solver
            % (and change which multiplier it returns), so add them once.
            multiplierConstraints{i_W} = [];
        else
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
    end
    multiplierConstraints = multiplierConstraints(~cellfun(@isempty, multiplierConstraints));

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

            % Square/SISO case: Psi = I - pinv(G21)*G21 = 0, so the last
            % size(Psi,2) rows/cols of M are identically zero. M then has no
            % strictly feasible point. Drop those rows/cols (M >= 0 is
            % mathematically unchanged, but the solver sees a different problem).
            if opts.dropNullRows && norm(Psi, 'fro') < 1e-9
                keep = 1:(size(M,1) - size(Psi,2));
                M = M(keep, keep);
            end

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
