function [constraints, slack, D, E, PI, Tau] = IQC_performance_combinedPI_par_fix(plant, PILTI, controller, eta, gamma, alpha, beta, omegas, opts)
    %IQC_PERFORMANCE_COMBINEDPI_PAR_FIX Robust-performance version of
    %IQC_reduced_combinedPI_par_fix: for a FIXED controller and a FIXED performance
    %level gamma, LMI constraints in the multiplier of the parametric block
    %    PI_m(w) = [eta^2 D(w), E(w); E(w)', -D(w)],  D = D' >= 0, E = -E'
    %combined with a GIVEN multiplier PILTI for the LTI block (same ordering and
    %options as IQC_reduced_combinedPI_par_fix), such that the robust performance
    %condition of datadriven.Objective.IQC_performance holds:
    %    || T_zw ||^2 <= gamma  for every uncertainty satisfying the IQC.
    %
    %   plant: three-port, InputGroup U1 = q, U2 = w, U3 = u; OutputGroup Y1 = p,
    %   Y2 = z, Y3 = y (as for IQC_performance).
    %   With K and gamma fixed, the condition of IQC_performance
    %       [PIrp_plus^-1, L; L', -frakL] >= 0
    %   is written in its congruence-transformed form (affine in D, E)
    %       [PIrp_plus, PIrp_plus*L; L'*PIrp_plus, -frakL] >= 0,
    %   PIrp_plus = blkdiag over (p, z | q, w) of PI_plus, gamma^-1 I, beta_perf I.
    %   No PI normalisation (use pinorm = false in the matching IQC_performance call).
    %   "addSlack": M + slack I >= 0, slack >= slackLower (negative = margin).
    %   "freeLTIScale" (true): the LTI block multiplier is tau(w)*PILTI with tau(w) >= epsTau
    %   a decision variable (valid for an LTI block), so the scale of the uncertainty
    %   multipliers against the performance channel is optimised instead of fixed (tau = 1).
    %   Tau (cell, one per frequency) is returned; PI contains tau.
    %
    %   Adapted 2026 from IQC_reduced_combinedPI_par_fix and IQC_performance.
    %

    arguments
        plant lti
        PILTI
        controller datadriven.Controller.Controller
        eta (1, 1) double
        gamma (1, 1) double
        alpha (1, 1) double
        beta (1, 1) double
        omegas (:, 1) double
        opts.epsNormalisation (1, 1) double = 1e-6;
        opts.epsSVD (1, 1) double = 1e-6;
        opts.addSlack (1, 1) logical = true;
        opts.ltiIdx double = [];
        opts.epsD (1, 1) double = 1e-4;
        opts.DMax (1, 1) double = Inf;
        opts.slackLower (1, 1) double = 0;
        opts.epsPIplus (1, 1) double = 1e-6;
        opts.multiplier (1, 1) string {mustBeMember(opts.multiplier, ["frequency", "constant"])} = "frequency";
        opts.dropNullRows (1, 1) logical = true;
        opts.warmstartD = []
        opts.warmstartE = []
        opts.warmstartPI double = [];   % previous COMBINED PI (numeric, 2n x 2n)
        opts.verbose (1, 1) logical = false;
        opts.betaPerformance (1, 1) double = 1e-2;
        opts.freeLTIScale (1, 1) logical = false;
        opts.epsTau (1, 1) double = 1e-4;
    end

    % Ensure that frequencies are ordered!
    [omegas, idxSort] = sort(omegas);
    nomegas = length(omegas);
    nq = length(plant.InputGroup.U1);
    np = length(plant.OutputGroup.Y1);
    nw = length(plant.InputGroup.U2);
    nz = length(plant.OutputGroup.Y2);
    ny = length(plant.OutputGroup.Y3);

    % Check if base assumption is satisfied
    if controller.factorisation == "Right"
        if ny > np
            error("Base assumption is not satisfied!!!");
        end
    elseif controller.factorisation == "Left"
        error("IQC not defined for Left factorisation! Yet!");
    end

    % The structure requires a square delta-channel
    if np ~= nq
        error("IQC_reduced_combinedPI:dimensionMismatch", ...
            "The PI structure requires np == nq (got np=%d, nq=%d).", np, nq);
    end
    n = np;

    % ----- Split the delta channels: LTI part (given) / optimised part -----
    nl = size(PILTI, 1) / 2;
    % nl = 0 (PILTI = zeros(0)): no LTI block, the optimised block takes all channels
    if nl ~= round(nl) || size(PILTI, 2) ~= size(PILTI, 1) || nl < 0 || nl >= n
        error("IQC_reduced_combinedPI:PILTISize", ...
            "PILTI must be (2*nl x 2*nl) with 0 <= nl < n = %d (got %d x %d).", ...
            n, size(PILTI, 1), size(PILTI, 2));
    end
    idxL = opts.ltiIdx;
    if isempty(idxL), idxL = 1:nl; end
    if numel(idxL) ~= nl || numel(unique(idxL)) ~= nl || any(idxL < 1) || ...
            any(idxL > n) || any(idxL ~= round(idxL))
        error("IQC_reduced_combinedPI:ltiIdx", ...
            "ltiIdx must hold %d distinct integers in 1..%d.", nl, n);
    end
    idxM = setdiff(1:n, idxL);
    nm = numel(idxM);                 % size of the optimised block
    I_n = eye(n);
    SL = I_n(:, idxL);                % n x nl   (selects LTI channels)
    SM = I_n(:, idxM);                % n x nm   (selects optimised channels)

    PL = local_PILTI_pages(PILTI, omegas, idxSort, nomegas, nl);

    % ----- Warm starts (always assigned: an unassigned sdpvar is NaN under double()) -----
    D_ws = [];
    E_ws = [];
    if ~isempty(opts.warmstartPI)
        PI_ws = opts.warmstartPI;
        if size(PI_ws, 1) ~= 2*n || size(PI_ws, 2) ~= 2*n
            error("IQC_reduced_combinedPI:warmstartPISize", ...
                "warmstartPI must be %dx%d (combined PI).", 2*n, 2*n);
        end
        D_from_PI = PI_ws(idxM, idxM) / eta^2;
        E_from_PI = PI_ws(idxM, n + idxM);
        D_ws = (D_from_PI + D_from_PI') / 2;
        E_ws = (E_from_PI - E_from_PI') / 2;
    end
    if ~isempty(opts.warmstartD), D_ws = opts.warmstartD; end
    if ~isempty(opts.warmstartE), E_ws = opts.warmstartE; end

    % Make them cells {nomegas x 1}
    if isempty(D_ws)
        D_ws = repmat({eye(nm)}, nomegas, 1);
    elseif ~iscell(D_ws)
        D_ws = repmat({D_ws}, nomegas, 1);
    end
    if isempty(E_ws)
        E_ws = repmat({zeros(nm)}, nomegas, 1);
    elseif ~iscell(E_ws)
        E_ws = repmat({E_ws}, nomegas, 1);
    end

    % ----- Optimised multiplier (parametric block): one D_k, E_k per frequency -----
    D = cell(nomegas, 1);
    E = cell(nomegas, 1);
    for i = 1:nomegas
        if opts.multiplier == "constant" && i > 1
            % one D, E shared by all frequencies
            D{i} = D{1};
            E{i} = E{1};
            continue
        end
        D{i} = sdpvar(nm, nm, 'hermitian', 'complex');
        E{i} = sdpvar(nm, nm, 'full', 'complex');

        warmstart(D{i}, D_ws{i});
        warmstart(E{i}, E_ws{i});
    end

    % ----- Scale of the LTI-block multiplier: tau(w) (fixed 1 unless freeLTIScale) -----
    Tau = repmat({1}, nomegas, 1);
    if opts.freeLTIScale
        for i = 1:nomegas
            if opts.multiplier == "constant" && i > 1
                Tau{i} = Tau{1};   continue
            end
            Tau{i} = sdpvar(1);
            warmstart(Tau{i}, 1);
        end
    end

    if opts.verbose
        fprintf('D warm-start:\n'); disp(double(D{1}));
        fprintf('E warm-start:\n'); disp(double(E{1}));
    end

    multiplierConstraints = cell(nomegas, 1);
    PI = cell(nomegas, 1);
    PI_plus = cell(nomegas, 1);
    PI_minus_22 = cell(nomegas, 1);
    nBadLTI = 0;
    for i_W = 1:nomegas
        Dk = D{i_W};
        Ek = E{i_W};

        tk   = Tau{i_W};
        PLk  = PL(:, :, i_W);
        PL11 = tk*PLk(1:nl,     1:nl);
        PL12 = tk*PLk(1:nl,     nl+1:end);
        PL21 = tk*PLk(nl+1:end, 1:nl);
        PL22 = tk*PLk(nl+1:end, nl+1:end);

        % Combined multiplier on (p ; q): LTI block (given) + parametric block (optimised)
        PI_11 = SL*PL11*SL' + SM*(eta^2 * Dk)*SM';
        PI_12 = SL*PL12*SL' + SM*Ek*SM';
        PI_21 = SL*PL21*SL' + SM*Ek'*SM';
        PI_22 = SL*PL22*SL' - SM*Dk*SM';
        PI{i_W} = [PI_11, PI_12; PI_21, PI_22];

        PI_plus{i_W} = [ ...
            PI_11 + alpha * eye(np) , PI_12                       ;
            PI_21                   , (beta + eps(1.0)) * eye(nq) ;
            ];
        PI_minus_22{i_W} = PI_22 - beta * eye(nq);

        % LTI part of PI_plus: checked when fixed, constrained when tau is free
        PLplus = [PL11 + alpha*eye(nl), PL12; PL21, (beta + eps(1.0))*eye(nl)];
        if ~opts.freeLTIScale && min(real(eig((PLplus + PLplus') / 2))) <= 0
            nBadLTI = nBadLTI + 1;
        end

        if opts.multiplier == "constant" && i_W > 1
            % Same variables as frequency 1 -> same constraints (add once).
            multiplierConstraints{i_W} = [];
        else
            % PI_plus is block diagonal in (LTI, parametric) after permutation, so
            % PI_plus >= 0 splits: constrain the optimised block only.
            PI_plus_m = [ ...
                eta^2 * Dk + alpha * eye(nm) , Ek                          ;
                Ek'                          , (beta + eps(1.0)) * eye(nm) ;
                ];
            multiplierConstraints{i_W} = [ ...
                Dk >= opts.epsD * eye(nm) ; ...
                PI_plus_m >= opts.epsPIplus * eye(2*nm) ; ...
                Ek + Ek' == zeros(nm) ; ...   % enforces Ek = -Ek'
                ];
            if isfinite(opts.DMax)
                multiplierConstraints{i_W} = [multiplierConstraints{i_W}; ...
                    Dk <= opts.DMax * eye(nm)];
            end
            if opts.freeLTIScale
                multiplierConstraints{i_W} = [multiplierConstraints{i_W}; ...
                    tk >= opts.epsTau; ...
                    (PLplus + PLplus')/2 >= opts.epsPIplus * eye(2*nl)];
            end
        end
    end
    multiplierConstraints = multiplierConstraints(~cellfun(@isempty, multiplierConstraints));

    if nBadLTI > 0
        warning("IQC_reduced_combinedPI:PILTIplusNotPD", ...
            ['[PILTI11+alpha*I, PILTI12; PILTI21, beta*I] is not positive definite at ' ...
             '%d of %d frequencies: the congruence transform needs PI_plus > 0. ' ...
             'Adjust alpha/beta or PILTI.'], nBadLTI, nomegas);
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
    p_idx = [plant.OutputGroup.Y1, plant.OutputGroup.Y2];   % [p; z]
    q_idx = [plant.InputGroup.U1,  plant.InputGroup.U2 ];   % [q; w]
    bp    = opts.betaPerformance;

    % ===== Constraints =====
    constraints_tmp = cell(nomegas, 1);
    for i_W = 1:nomegas
        [X, Y] = freqresp(controller, omegas(i_W));
        if i_W == 1 && (isa(X, 'sdpvar') || isa(Y, 'sdpvar'))
            warning("IQC_performance_combinedPI:controllerNotFixed", ...
                'Controller frequency response still contains free decision variables: fix the controller.');
        end
        G11 = Gf(p_idx, q_idx, i_W);
        G12 = Gf(p_idx, plant.InputGroup.U3, i_W);
        G21 = Gf(plant.OutputGroup.Y3, q_idx, i_W);
        G22 = Gf(plant.OutputGroup.Y3, plant.InputGroup.U3, i_W);

        if controller.factorisation ~= "Right"
            error("IQC not defined for Left factorisation! Yet!");
        end
        if (opts.epsSVD >= 0)
            s = svds(G21, 1, "smallest");
            G21 = G21 / s;
            G12 = G12 * s;
        end
        Phi   = double(pinv(G21) * (Y - G22*X));   % numeric (K fixed)
        Psi   = eye(nq + nw) - pinv(G21)*G21;
        L     = double([G11*Phi + G12*X, G11*Psi;
                        Phi,              Psi    ]);

        % PIrp_plus on (p, z | q, w): PI_plus, 1/gamma I, beta_perf I
        Pp  = PI_plus{i_W};
        P11 = Pp(1:np, 1:np);   P12 = Pp(1:np, np+1:end);
        P21 = Pp(np+1:end, 1:np);   P22 = Pp(np+1:end, np+1:end);
        PIrp_plus = [blkdiag(P11, eye(nz)/gamma), blkdiag(P12, zeros(nz, nw));
                     blkdiag(P21, zeros(nw, nz)), blkdiag(P22, bp*eye(nw))];
        tmp = blkdiag(PI_minus_22{i_W}, -(1 + bp)*eye(nw));

        Gp = G11*Phi + G12*X;
        frakL11 = Phi'*tmp*Phi;
        frakL12 = Phi'*tmp*Psi - alpha * Gp' * (G11*Psi);
        frakL22 = Psi'*tmp*Psi - alpha * (G11*Psi)' * (G11*Psi);
        frakL   = [frakL11, frakL12; frakL12', frakL22];

        M = [PIrp_plus, PIrp_plus*L; L'*PIrp_plus, -frakL];
        M = M + slack * eye(size(M));
        if opts.epsNormalisation > 0
            M_num = double(M);
            if any(~isfinite(M_num(:)))
                M = (M + M') / 2;
            else
                [s, v, d] = svd(M_num);
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

function PL = local_PILTI_pages(PILTI, omegas, idxSort, nW, nl)
    %LOCAL_PILTI_PAGES  PILTI as a (2nl x 2nl x nW) numeric array at omegas.
    if isa(PILTI, 'DynamicSystem')            % frd, ss, tf, ...
        PL = freqresp(PILTI, omegas);
    elseif isnumeric(PILTI)
        if ndims(PILTI) < 3 || size(PILTI, 3) == 1
            PL = repmat(PILTI, [1 1 nW]);
        elseif size(PILTI, 3) == nW
            PL = PILTI(:, :, idxSort);        % follow the sorted frequency order
        else
            error("IQC_reduced_combinedPI:PILTIPages", ...
                "Numeric PILTI has %d pages but omegas has %d entries.", size(PILTI, 3), nW);
        end
    else
        error("IQC_reduced_combinedPI:PILTIType", ...
            "PILTI must be an frd/DynamicSystem or a numeric array.");
    end
    PL = double(PL);
    if ~isequal(size(PL, [1 2]), [2*nl, 2*nl]) || size(PL, 3) ~= nW
        error("IQC_reduced_combinedPI:PILTIPages", ...
            "PILTI evaluated to size %s, expected %dx%dx%d.", mat2str(size(PL)), 2*nl, 2*nl, nW);
    end
    % A multiplier must be Hermitian; symmetrise numerical noise, warn on real asymmetry.
    asym = 0;
    for k = 1:nW
        Pk = PL(:, :, k);
        asym = max(asym, norm(Pk - Pk', 'fro') / max(norm(Pk, 'fro'), eps));
        PL(:, :, k) = (Pk + Pk') / 2;
    end
    if asym > 1e-8
        warning("IQC_reduced_combinedPI:PILTINotHermitian", ...
            "PILTI is not Hermitian (relative asymmetry %.2e); it was symmetrised.", asym);
    end
end
