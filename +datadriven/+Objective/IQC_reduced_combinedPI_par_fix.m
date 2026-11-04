function [constraints, slack, D, E, PI] = IQC_reduced_combinedPI_par_fix(plant, PILTI, controller, eta, alpha, beta, omegas, opts)
    %IQC_REDUCED_COMBINEDPI_PAR_FIX Slacked, reduced IQC constraint where the
    %multiplier is the COMBINATION of
    %   - a GIVEN multiplier PILTI(w) for the (kind of) LTI uncertainty, and
    %   - an OPTIMISED full-block multiplier for the parametric uncertainty
    %       PI_m(w) = [ eta^2*D(w) ,  E(w)  ;
    %                     E(w)'    , -D(w)  ]
    %     with D = D' >= 0, E = -E', 'eta' a FIXED bound on the parameter.
    %
    %   ------------------------------------------------------------------------
    %   ORDERING
    %   ------------------------------------------------------------------------
    %   The delta channel has n = np = nq signals: p = plant outputs Y1 (into
    %   Delta, size np) and q = plant inputs U1 (out of Delta, size nq). The
    %   multiplier acts on [p ; q] (p FIRST), e.g. for n = 2
    %       PI = [ PI_11  PI_12 ;      acting on (p1, p2, q1, q2)
    %              PI_21  PI_22 ]
    %   p1,q1 belong to PILTI, p2,q2 to the parametric block, and there are
    %   no cross terms between the two blocks (E only couples p2 and q2).
    %   PILTI is the multiplier of the LTI block ALONE, ordered (p_lti ; q_lti):
    %       PILTI = [ PL11 PL12 ; PL21 PL22 ]      (2*nl x 2*nl)
    %   and is embedded on the channels "ltiIdx" (default 1:nl). The remaining
    %   channels belong to the optimised block (size nm = n - nl):
    %       PI_11 = S_L*PL11*S_L' + S_M*(eta^2*D)*S_M'
    %       PI_12 = S_L*PL12*S_L' + S_M*E*S_M'
    %       PI_21 = S_L*PL21*S_L' + S_M*E'*S_M'
    %       PI_22 = S_L*PL22*S_L' - S_M*D*S_M'
    %   With PILTI = diag(1,-1), D = 1, E = 0 this is diag(1, eta^2, -1, -1).
    %   Check the block order of your LFT with the 3rd output of lftdata
    %   ([M,Delta,BlkStruct] = lftdata(...)); if the parametric block comes
    %   first pass "ltiIdx", 2.
    %
    %   PILTI can be an frd / any DynamicSystem (evaluated at omegas), a numeric
    %   2nl x 2nl (same at every frequency) or a numeric 2nl x 2nl x nomegas
    %   array (page k belongs to omegas(k) as given by the caller).
    %
    %   D and E are returned as CELL ARRAYS (nomegas x 1) of the optimised
    %   block only (nm x nm). PI is a cell array of the COMBINED multiplier
    %   (2n x 2n each).
    %   With "multiplier","constant" a single D, E is shared by all
    %   frequencies and the multiplier constraints are added ONCE.
    %
    %   The controller MUST be fixed (numeric) for this call: the congruence
    %   transformed LMI
    %       M' = [PI_plus, PI_plus*L; L'*PI_plus, -frakL] >= 0
    %   is affine in (D,E) only if L and frakL's controller parts are numeric.
    %   It also needs PI_plus > 0, which for the LTI block means
    %   [PL11+alpha*I, PL12; PL21, beta*I] > 0 (a warning is raised if not).
    %
    %   ------------------------------------------------------------------------
    %   PURE FEASIBILITY (addSlack = false) vs SLACK
    %   ------------------------------------------------------------------------
    %   For a square/SISO delta-channel Psi = 0 and the last rows/cols of M
    %   are identically zero. With "dropNullRows" true (default) they are
    %   removed (only happens when Psi = 0). NOTE: without an objective the
    %   multiplier returned by the solver is NOT unique, and which point you
    %   get depends on the formulation.
    %
    %   With "addSlack" true the LMI is  M + slack*I >= 0.
    %     - "slackLower" = 0 (default): slack >= 0 (violation measure).
    %     - "slackLower" < 0: slack may be negative -> margin.
    %
    %   Name-value options:
    %     "epsNormalisation", "epsSVD", "addSlack" - as before.
    %     "ltiIdx"       - channels (1..n) of the delta channel owned by PILTI.
    %                      Default 1:nl.
    %     "multiplier"   - "frequency" (default) or "constant".
    %     "dropNullRows" - drop structurally zero rows/cols of M. Default true.
    %     "epsD"         - D >= epsD*I. Default 1e-4.
    %     "DMax"         - D <= DMax*I. Default Inf (no bound).
    %     "slackLower"   - slack >= slackLower. Default 0.
    %     "epsPIplus"    - PI_plus_m >= epsPIplus*I (optimised block). Default 1e-6.
    %     "warmstartD", "warmstartE" - numeric nm x nm, or cell {nomegas x 1}.
    %     "warmstartPI"  - previous COMBINED PI (numeric 2n x 2n), split back
    %                      into D, E.
    %
    %   Note: the "stability" LMI of IQC_reduced is NOT added here - with a
    %   fixed controller it does not depend on D, E (it is a numeric constant).
    %
    %   References
    %     [1] V. Gupta, E. Klauser, and A. Karimi, "Non-parametric IQC
    %     multipliers in data-driven robust controller synthesis," 2024.
    %
    %   Copyright 2024 Vaibhav Gupta, DDMAC, EPFL (MIT License)
    %   Adapted 2026 to optimise the IQC multiplier itself and to combine it
    %   with a given PILTI.
    %

    arguments
        plant lti
        PILTI
        controller datadriven.Controller.Controller
        eta (1, 1) double
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
    end

    % Ensure that frequencies are ordered!
    [omegas, idxSort] = sort(omegas);
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

    % The structure requires a square delta-channel
    if np ~= nq
        error("IQC_reduced_combinedPI:dimensionMismatch", ...
            "The PI structure requires np == nq (got np=%d, nq=%d).", np, nq);
    end
    n = np;

    % ----- Split the delta channels: LTI part (given) / optimised part -----
    nl = size(PILTI, 1) / 2;
    if nl ~= round(nl) || size(PILTI, 2) ~= size(PILTI, 1) || nl < 1 || nl >= n
        error("IQC_reduced_combinedPI:PILTISize", ...
            "PILTI must be (2*nl x 2*nl) with 1 <= nl < n = %d (got %d x %d).", ...
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

        PLk  = PL(:, :, i_W);
        PL11 = PLk(1:nl,     1:nl);
        PL12 = PLk(1:nl,     nl+1:end);
        PL21 = PLk(nl+1:end, 1:nl);
        PL22 = PLk(nl+1:end, nl+1:end);

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

        % The LTI part of PI_plus is a constant: it is only CHECKED, never constrained.
        PLplus = [PL11 + alpha*eye(nl), PL12; PL21, (beta + eps(1.0))*eye(nl)];
        if min(real(eig((PLplus + PLplus') / 2))) <= 0
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

    % ===== Constraints =====
    constraints_tmp = cell(nomegas, 1);

    for i_W = 1:nomegas

        % Prepare different values
        [X, Y] = freqresp(controller, omegas(i_W));

        if i_W == 1 && (isa(X, 'sdpvar') || isa(Y, 'sdpvar'))
            warning("IQC_reduced_combinedPI:controllerNotFixed", ...
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
