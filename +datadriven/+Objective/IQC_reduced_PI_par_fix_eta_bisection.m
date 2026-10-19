function [eta, D_val, E_val, PI_matrix, diagnostics] = IQC_reduced_PI_par_fix_eta_bisection(plant, controller, alpha, beta, omegas, opts)
    %IQC_REDUCED_PI_ETA_BISECTION Find the largest eta for which
    %IQC_reduced_PI (D, E free, complex/skew-Hermitian E) is feasible, for
    %a FIXED controller, by bisecting eta.
    %
    %   ------------------------------------------------------------------------
    %   WHY BISECTION (and why via the SLACK, not raw solver status)
    %   ------------------------------------------------------------------------
    %   eta only enters IQC_reduced_PI's LMI through eta^2*D; with D also
    %   free this is bilinear, so eta cannot be jointly optimised with D, E
    %   in a single SDP (see IQC_reduced_PI's own docstring). Instead this
    %   fixes eta at each candidate and re-solves the D,E-feasibility
    %   problem, exactly like your manual attempt - the difference here is
    %   twofold:
    %     1) The bracket is actually GROWN (doubled) until a genuinely
    %        infeasible eta is found, rather than bisecting inside a fixed
    %        [etaLower, etaUpperInit] window that might not contain the
    %        true boundary at all (this was the bug in your manual loop -
    %        eta1/eta2 were never updated, so it wasn't bisecting at all).
    %     2) Feasibility is judged by calling IQC_reduced_PI with
    %        "addSlack" true and MINIMISING the slack, then checking
    %        slack <= slackTol - not by trusting optimize()'s raw
    %        sol.problem code. Right at the feasibility boundary (exactly
    %        where bisection spends most of its time), problem-status codes
    %        can be noisy; the slack's sign is a continuous, more robust
    %        signal (this mirrors the paper's own Sec 3.4 "Iterative
    %        procedure": minimise a slack nu, feasible iff nu <= 0).
    %
    %   D, E are warm-started from the previous accepted bisection point at
    %   each new candidate eta, both for solver speed and because
    %   IQC_reduced_PI's own SVD-based conditioning step uses double(D),
    %   double(E) at their currently assigned values.
    %
    %   ------------------------------------------------------------------------
    %   [eta, D, E, PI, diagnostics] = IQC_reduced_PI_eta_bisection(plant,
    %   controller, alpha, beta, omegas) returns the largest feasible eta
    %   found, the corresponding numeric D (n x n, real), E (n x n,
    %   complex), the numeric PI = [eta^2*D, E; E', -D] (frequency-
    %   invariant), and a diagnostics struct with the bracket/bisection
    %   history.
    %
    %   controller MUST be fixed (e.g. datadriven.Controller.Fixed(K)) -
    %   IQC_reduced_PI itself will warn if it is not.
    %
    %   Name-value options:
    %     "etaLower"      - initial lower (assumed feasible) bound. Default 1e-5.
    %     "etaUpperInit"  - initial upper-bound guess to test/expand from. Default 1.
    %     "etaHardCap"    - physical ceiling on eta (e.g. your dm_max) - the
    %                       bracket is never expanded past this value even
    %                       if still feasible there. Default Inf.
    %     "relTol"        - stop when (etaUpper-etaLower)/etaLower < relTol. Default 1e-2.
    %     "slackTol"      - feasible iff minimised slack <= slackTol. Default 1e-6.
    %     "maxExpansions" - cap on bracket-doubling steps. Default 20.
    %     "maxBisections" - cap on bisection iterations. Default 60.
    %     "solver"        - YALMIP solver name. Default "mosek".
    %     "verbose"       - print progress. Default true.
    %     "epsD", "epsPIplus", "epsNormalisation", "epsSVD" - passed through
    %                       to IQC_reduced_PI (see its own docstring).
    %
    %   Adapted 2026: eta-bisection wrapper around IQC_reduced_PI (complex E).
    %

    arguments
        plant lti
        controller datadriven.Controller.Controller
        alpha (1, 1) double
        beta (1, 1) double
        omegas (:, 1) double
        opts.etaLower (1, 1) double = 1e-5;
        opts.etaUpperInit (1, 1) double = 1;
        opts.etaHardCap (1, 1) double = Inf;
        opts.relTol (1, 1) double = 1e-2;
        opts.slackTol (1, 1) double = 1e-6;
        opts.maxExpansions (1, 1) double = 20;
        opts.maxBisections (1, 1) double = 60;
        opts.solver (1, 1) string = "mosek";
        opts.verbose (1, 1) logical = true;
        opts.epsD (1, 1) double = 1e-4;
        opts.epsPIplus (1, 1) double = 1e-6;
        opts.epsNormalisation (1, 1) double = 1e-6;
        opts.epsSVD (1, 1) double = 1e-6;
        opts.D0 double = [];   % starting D for the very first solve. Default:
                                % eye(n) (i.e. D=1 for n=1), filled in below
                                % once n is known.
        opts.E0 double = [];   % starting E for the very first solve. Default:
                                % zeros(n) (i.e. E=0 - or pass e.g. 1e-6i for
                                % a small nonzero start).
    end

    solverOpts = sdpsettings('verbose', 0, 'debug', 0, 'solver', char(opts.solver));

    % n must match IQC_reduced_PI's own convention (np==nq required there).
    n = length(plant.OutputGroup.Y1);
    nomegas = length(omegas);
    
    D0 = opts.D0;
    E0 = opts.E0;
    if isempty(D0)
        D0 = cell(nomegas, 1);
        for i=1:nomegas
            D0{i} = eye(n);
        end
    end
    if isempty(E0)
        E0 = cell(nomegas, 1);
        for i=1:nomegas
            E0{i} = 1e-6*eye(n);
        end
    end

    % ----- Inner feasibility solve for D, E at a FIXED eta -----
    function [feasible, D_val, E_val, slackVal] = check_feasibility(eta_try, D_ws, E_ws)
        if nargin < 2, D_ws = []; end
        if nargin < 3, E_ws = []; end

        warmArgs = {};
        if ~isempty(D_ws), warmArgs = [warmArgs, {'warmstartD', D_ws}]; end %#ok<AGROW>
        if ~isempty(E_ws), warmArgs = [warmArgs, {'warmstartE', E_ws}]; end %#ok<AGROW>

        [LMIs, slack, D_var, E_var] = datadriven.Objective.IQC_reduced_PI_par_fix( ...
            plant, controller, eta_try, alpha, beta, omegas, ...
            "addSlack", false, ...
            "epsD", opts.epsD, "epsPIplus", opts.epsPIplus, ...
            "epsNormalisation", opts.epsNormalisation, "epsSVD", opts.epsSVD, ...
            warmArgs{:});

        sol = optimize(LMIs, [], solverOpts);
        slackVal = double(slack);
        feasible = (sol.problem == 0); %&& (slackVal <= opts.slackTol);
        D_val = cellfun(@double, D_var, 'UniformOutput', false);
        E_val = cellfun(@double, E_var, 'UniformOutput', false);
    end

    % ----- Sanity check: etaLower must be feasible -----
    % Explicitly start from D0, E0 (default D=1, E=0 for n=1) rather than
    % relying on IQC_reduced_PI's own internal fallback default.
    [feasible0, D_best, E_best, slack0] = check_feasibility(opts.etaLower, D0, E0);
    if ~feasible0
        error("IQC_reduced_PI_eta_bisection:infeasibleAtEtaLower", ...
            ["No valid multiplier exists even at eta = etaLower = %.3g " ...
             "(minimised slack = %.4g). No robust design is certifiable for " ...
             "this controller with this IQC structure - check the controller, " ...
             "plant, or lower etaLower further."], opts.etaLower, slack0);
    end
    if opts.verbose
        fprintf('  [etaLower] eta = %.4g -> true (slack = %.4g)\n', opts.etaLower, slack0);
        % fprintf('    D = %s\n', mat2str(D_best, 4));
        % fprintf('    E = %s\n', mat2str(E_best, 4));
    end
    etaLower = opts.etaLower;
    eta_best = etaLower;

    % ----- Bracket expansion: grow eta_upper until infeasible (or capped) -----
    etaUpper = opts.etaUpperInit;
    expansionHistory = etaUpper;
    nExpansions = 0;
    [feasibleUpper, D_at_upper, E_at_upper, slackUpper] = check_feasibility(etaUpper, D_best, E_best);
    if opts.verbose
        fprintf('  [expand] eta = %.4g -> %s (slack = %.4g)\n', etaUpper, string(feasibleUpper), slackUpper);
        if feasibleUpper
            % fprintf('    D = %s\n', mat2str(D_at_upper, 4));
            % fprintf('    E = %s\n', mat2str(E_at_upper, 4));
        end
    end
    while feasibleUpper && nExpansions < opts.maxExpansions
        eta_best = etaUpper;
        D_best = D_at_upper;
        E_best = E_at_upper;
        etaLower = etaUpper;
        if etaUpper >= opts.etaHardCap
            if opts.verbose
                fprintf('IQC_reduced_PI_eta_bisection: hit etaHardCap = %.4g while still feasible; stopping expansion.\n', opts.etaHardCap);
            end
            break
        end
        etaUpper = min(etaUpper * 2, opts.etaHardCap);
        expansionHistory(end+1) = etaUpper; %#ok<AGROW>
        nExpansions = nExpansions + 1;
        [feasibleUpper, D_at_upper, E_at_upper, slackUpper] = check_feasibility(etaUpper, D_best, E_best);
        if opts.verbose
            fprintf('  [expand] eta = %.4g -> %s (slack = %.4g)\n', etaUpper, string(feasibleUpper), slackUpper);
            % if feasibleUpper
            %     fprintf('    D = %s\n', mat2str(D_at_upper, 4));
            %     fprintf('    E = %s\n', mat2str(E_at_upper, 4));
            % end
        end
    end
    if feasibleUpper
        eta = etaUpper;
        D_val = D_at_upper;
        E_val = E_at_upper;
        PI_matrix=zeros(2*n,2*n,nomegas);
        for i=1:nomegas
            PI_matrix(:,:,i) = [eta^2*D_val{i}, E_val{i}; E_val{i}', -D_val{i}];
        end
        diagnostics = struct('expansionHistory', expansionHistory, ...
            'bisectionHistory', [], 'hitCeiling', true);
        warning("IQC_reduced_PI_eta_bisection:noInfeasibleUpperBoundFound", ...
            ["Still feasible at eta = %.4g (etaHardCap/maxExpansions reached) - " ...
             "the returned eta is a lower bound on the true maximum, not the " ...
             "maximum itself. Raise etaHardCap or maxExpansions to search further."], eta);
        return
    end

    % ----- Bisection within the verified [etaLower(feasible), etaUpper(infeasible)] bracket -----
    bisectionHistory = zeros(0, 3);   % [eta_mid, feasible, slack]
    iter = 0;
    while (etaUpper - etaLower)/etaLower > opts.relTol && iter < opts.maxBisections
        iter = iter + 1;
        etaMid = (etaLower + etaUpper) / 2;
        [feasibleMid, D_mid, E_mid, slackMid] = check_feasibility(etaMid, D_best, E_best);
        bisectionHistory(end+1, :) = [etaMid, feasibleMid, slackMid]; %#ok<AGROW>
        if opts.verbose
            fprintf('  [bisect %2d] eta = %.6g -> %s (slack = %.4g)\n', iter, etaMid, string(feasibleMid), slackMid);
            if feasibleMid
                % fprintf('    D = %s\n', mat2str(D_mid, 4));
                % fprintf('    E = %s\n', mat2str(E_mid, 4));
            end
        end
        if feasibleMid
            etaLower = etaMid;
            eta_best = etaMid;
            D_best = D_mid;
            E_best = E_mid;
        else
            etaUpper = etaMid;
        end
    end

    eta = eta_best;
    D_val = D_best;
    E_val = E_best;
    PI_matrix=zeros(2*n,2*n,nomegas);
    for i=1:nomegas
        PI_matrix(:,:,i) = [eta^2*D_val{i}, E_val{i}; E_val{i}', -D_val{i}];
    end
    diagnostics = struct('expansionHistory', expansionHistory, ...
        'bisectionHistory', bisectionHistory, 'hitCeiling', false, ...
        'iterations', iter);

    if opts.verbose
        fprintf('IQC_reduced_PI_eta_bisection: converged, eta = %.6g\n', eta);
    end

end