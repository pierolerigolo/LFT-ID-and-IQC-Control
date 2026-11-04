function [eta, D_val, E_val, PI_matrix, diagnostics] = IQC_reduced_combinedPI_par_fix_eta_bisection(plant, PILTI, controller, alpha, beta, omegas, opts)
    %IQC_REDUCED_COMBINEDPI_PAR_FIX_ETA_BISECTION Find the largest eta (bound of
    %the parametric uncertainty) for which an optimised multiplier D(w), E(w)
    %exists, TOGETHER with a given PILTI for the LTI part, for a FIXED
    %controller, by bisecting eta.
    %
    %   PILTI is the multiplier of the (kind of) LTI block alone, ordered
    %   (p_lti ; q_lti), as frd / numeric (2nl x 2nl [x nomegas]). It is combined
    %   automatically with the optimised parametric block into the full PI in
    %   (p1,p2,...,q1,q2,...) ordering - see IQC_reduced_combinedPI_par_fix.
    %
    %   At every candidate eta the multiplier problem is solved either as a
    %   pure feasibility problem (default) or, with "useSlack", with a SLACK
    %   that is allowed to be negative and is MINIMISED (D <= DMax*I keeps it
    %   bounded):
    %       min slack   s.t.   M(w) + slack*I >= 0   for all w
    %   - slack <= slackTol  -> a multiplier exists at this eta (feasible)
    %   - slack <  0         -> it exists with room to spare
    %   - slack >  slackTol  -> infeasible
    %
    %   Name-value options:
    %     "ltiIdx"        - channels of the delta channel owned by PILTI. Default 1:nl.
    %     "etaLower"      - initial lower (assumed feasible) bound. Default 1e-5.
    %     "etaUpperInit"  - initial upper-bound guess to test/expand from. Default 1.
    %     "etaHardCap"    - ceiling on eta. Default Inf.
    %     "relTol"        - stop when (etaUpper-etaLower)/etaLower < relTol. Default 1e-2.
    %     "slackTol"      - feasible iff minimised slack <= slackTol. Default 1e-6.
    %     "useSlack"      - solve with the slack. Default false.
    %     "slackLower"    - lower bound on the slack. Default -1.
    %     "DMax"          - D <= DMax*I. Default Inf.
    %     "multiplier"    - "frequency" or "constant" (one D,E for all w).
    %     "dropNullRows"  - drop structurally zero rows/cols of M. Default true.
    %     "maxExpansions", "maxBisections", "solver", "verbose" - as before.
    %     "epsD", "epsPIplus", "epsNormalisation", "epsSVD" - passed through.
    %     "D0", "E0"      - warm starts for the first solve (numeric or cell,
    %                       size of the optimised block).
    %
    %   Outputs D_val, E_val are cell arrays (nomegas x 1) of the optimised
    %   block; PI_matrix is the COMBINED multiplier (2n x 2n x nomegas).
    %

    arguments
        plant lti
        PILTI
        controller datadriven.Controller.Controller
        alpha (1, 1) double
        beta (1, 1) double
        omegas (:, 1) double
        % NB: PI_plus_m >= epsPIplus*I needs eta^2*D >= epsPIplus, and D <= DMax,
        % so etaLower must satisfy etaLower^2*DMax >= epsPIplus.
        opts.ltiIdx double = [];
        opts.etaLower (1, 1) double = 1e-5;
        opts.etaUpperInit (1, 1) double = 1;
        opts.etaHardCap (1, 1) double = Inf;
        opts.relTol (1, 1) double = 1e-2;
        opts.slackTol (1, 1) double = 1e-6;
        opts.useSlack (1, 1) logical = false;   % false: pure feasibility
        opts.multiplier (1, 1) string = "frequency";   % or "constant" (one D,E for all w)
        opts.dropNullRows (1, 1) logical = true;
        opts.slackLower (1, 1) double = -1;
        opts.DMax (1, 1) double = Inf;
        opts.maxExpansions (1, 1) double = 20;
        opts.maxBisections (1, 1) double = 60;
        opts.solver (1, 1) string = "mosek";
        opts.verbose (1, 1) logical = true;
        opts.epsD (1, 1) double = 1e-2;
        opts.epsPIplus (1, 1) double = 1e-6;
        opts.epsNormalisation (1, 1) double = 1e-6;
        opts.epsSVD (1, 1) double = 1e-6;
        opts.D0 = [];
        opts.E0 = [];
    end

    solverOpts = sdpsettings('verbose', 0, 'debug', 0, 'solver', char(opts.solver));

    n  = length(plant.OutputGroup.Y1);
    nl = size(PILTI, 1) / 2;            % channels owned by PILTI
    nm = n - nl;                        % size of the optimised block
    nomegas = length(omegas);

    D0 = opts.D0;
    E0 = opts.E0;
    if isempty(D0), D0 = repmat({eye(nm)}, nomegas, 1); end
    if isempty(E0), E0 = repmat({zeros(nm)}, nomegas, 1); end

    % ----- Inner solve for D, E at a FIXED eta -----
    function [feasible, D_val, E_val, slackVal, PI_val] = check_feasibility(eta_try, D_ws, E_ws)
        if nargin < 2, D_ws = []; end
        if nargin < 3, E_ws = []; end

        warmArgs = {};
        if ~isempty(D_ws), warmArgs = [warmArgs, {'warmstartD', D_ws}]; end %#ok<AGROW>
        if ~isempty(E_ws), warmArgs = [warmArgs, {'warmstartE', E_ws}]; end %#ok<AGROW>

        [LMIs, slack, D_var, E_var, PI_var] = datadriven.Objective.IQC_reduced_combinedPI_par_fix( ...
            plant, PILTI, controller, eta_try, alpha, beta, omegas, ...
            "addSlack", opts.useSlack, "slackLower", opts.slackLower, "DMax", opts.DMax, ...
            "ltiIdx", opts.ltiIdx, ...
            "multiplier", opts.multiplier, "dropNullRows", opts.dropNullRows, ...
            "verbose", opts.verbose, ...
            "epsD", opts.epsD, "epsPIplus", opts.epsPIplus, ...
            "epsNormalisation", opts.epsNormalisation, "epsSVD", opts.epsSVD, ...
            warmArgs{:});

        if opts.useSlack
            sol = optimize(LMIs, slack, solverOpts);
            slackVal = value(slack);
            feasible = (sol.problem == 0) && isfinite(slackVal) && (slackVal <= opts.slackTol);
        else
            % Pure feasibility.
            sol = optimize(LMIs, [], solverOpts);
            slackVal = 0;
            feasible = (sol.problem == 0);
        end
        D_val  = cellfun(@double, D_var,  'UniformOutput', false);
        E_val  = cellfun(@double, E_var,  'UniformOutput', false);
        PI_val = cellfun(@double, PI_var, 'UniformOutput', false);   % combined PI per frequency
    end

    % ----- Sanity check: etaLower must be feasible -----
    [feasible0, D_best, E_best, slack0, PI_best] = check_feasibility(opts.etaLower, D0, E0);
    if ~feasible0
        error("IQC_reduced_combinedPI_eta_bisection:infeasibleAtEtaLower", ...
            "No valid multiplier at eta = etaLower = %.3g (minimised slack = %.4g). " + ...
            "Check that etaLower^2*DMax >= epsPIplus, then PILTI / the controller / plant.", ...
            opts.etaLower, slack0);
    end
    if opts.verbose
        fprintf('  [etaLower] eta = %.4g -> true (slack = %.4g)\n', opts.etaLower, slack0);
        fprintf('    D = %s\n', mat2str(D_best{1}, 4));
        fprintf('    E = %s\n', mat2str(E_best{1}, 4));
    end
    etaLower = opts.etaLower;
    eta_best = etaLower;

    % ----- Bracket expansion: grow eta_upper until infeasible (or capped) -----
    etaUpper = opts.etaUpperInit;
    expansionHistory = etaUpper;
    nExpansions = 0;
    [feasibleUpper, D_at_upper, E_at_upper, slackUpper, PI_at_upper] = check_feasibility(etaUpper, D_best, E_best);
    if opts.verbose
        fprintf('  [expand] eta = %.4g -> %s (slack = %.4g)\n', etaUpper, string(feasibleUpper), slackUpper);
        if feasibleUpper
            fprintf('    D = %s\n', mat2str(D_at_upper{1}, 4));
            fprintf('    E = %s\n', mat2str(E_at_upper{1}, 4));
        end
    end
    while feasibleUpper && nExpansions < opts.maxExpansions
        eta_best = etaUpper;
        D_best = D_at_upper;
        E_best = E_at_upper;
        PI_best = PI_at_upper;
        etaLower = etaUpper;
        if etaUpper >= opts.etaHardCap
            if opts.verbose
                fprintf('IQC_reduced_combinedPI_eta_bisection: hit etaHardCap = %.4g while still feasible; stopping expansion.\n', opts.etaHardCap);
            end
            break
        end
        etaUpper = min(etaUpper * 2, opts.etaHardCap);
        expansionHistory(end+1) = etaUpper; %#ok<AGROW>
        nExpansions = nExpansions + 1;
        [feasibleUpper, D_at_upper, E_at_upper, slackUpper, PI_at_upper] = check_feasibility(etaUpper, D_best, E_best);
        if opts.verbose
            fprintf('  [expand] eta = %.4g -> %s (slack = %.4g)\n', etaUpper, string(feasibleUpper), slackUpper);
            if feasibleUpper
                fprintf('    D = %s\n', mat2str(D_at_upper{1}, 4));
                fprintf('    E = %s\n', mat2str(E_at_upper{1}, 4));
            end
        end
    end
    if feasibleUpper
        eta = etaUpper;
        D_val = D_at_upper;
        E_val = E_at_upper;
        PI_matrix = cat(3, PI_at_upper{:});
        diagnostics = struct('expansionHistory', expansionHistory, ...
            'bisectionHistory', [], 'hitCeiling', true);
        warning("IQC_reduced_combinedPI_eta_bisection:noInfeasibleUpperBoundFound", ...
            "Still feasible at eta = %.4g (etaHardCap/maxExpansions reached): " + ...
            "returned eta is a lower bound on the true maximum.", eta);
        return
    end

    % ----- Bisection within the verified [etaLower(feasible), etaUpper(infeasible)] bracket -----
    bisectionHistory = zeros(0, 3);   % [eta_mid, feasible, slack]
    iter = 0;
    while (etaUpper - etaLower)/etaLower > opts.relTol && iter < opts.maxBisections
        iter = iter + 1;
        etaMid = (etaLower + etaUpper) / 2;
        [feasibleMid, D_mid, E_mid, slackMid, PI_mid] = check_feasibility(etaMid, D_best, E_best);
        bisectionHistory(end+1, :) = [etaMid, feasibleMid, slackMid]; %#ok<AGROW>
        if opts.verbose
            fprintf('  [bisect %2d] eta = %.6g -> %s (slack = %.4g)\n', iter, etaMid, string(feasibleMid), slackMid);
            if feasibleMid
                fprintf('    D = %s\n', mat2str(D_mid{1}, 4));
                fprintf('    E = %s\n', mat2str(E_mid{1}, 4));
            end
        end
        if feasibleMid
            etaLower = etaMid;
            eta_best = etaMid;
            D_best = D_mid;
            E_best = E_mid;
            PI_best = PI_mid;
        else
            etaUpper = etaMid;
        end
    end

    eta = eta_best;
    D_val = D_best;
    E_val = E_best;
    PI_matrix = cat(3, PI_best{:});
    diagnostics = struct('expansionHistory', expansionHistory, ...
        'bisectionHistory', bisectionHistory, 'hitCeiling', false, ...
        'iterations', iter);

    if opts.verbose
        fprintf('IQC_reduced_combinedPI_eta_bisection: converged, eta = %.6g\n', eta);
    end

end
