function model = find_elliptical_uncertainty(frd_models, w, options)
%find_elliptical_uncertainty Find the elliptical uncertainty set
% 
%   model = find_elliptical_uncertainty(frd_models) returns a structure 
%       containing the nominal model and the uncertainty set. Here,`frd_models` 
%       is the array of the FRF models under consideration.
% 
%   model = find_elliptical_uncertainty(frd_models, w) returns a structure
%       containing the nominal model and the uncertainty set computed at `w` 
%       (in [rad/s]).
% 
%   model = find_elliptical_uncertainty(..., "G", G) returns a structure with
%       the nominal model set as `G`.
% 
%   model = find_elliptical_uncertainty(..., "solver", solver) sets the solver
%       to use for convex optimisation. Tested solvers are: Mosek (default) and
%       SeDuMi.
% 
%   model = find_elliptical_uncertainty(..., "verbosity", level) sets the
%       verbosity level.
% 
%   ----------------------------------------------------------------------------
%   Output:
% 
%     model
%       ├ NominalModel              % Estimated/fixed nominal model
%       ├ UncertaintySet
%       │  ├ Type = "elliptical"    % Type of uncertainty set
%       │  ├ A                      % Cell array of matrix A for the set
%       │  ├ Area                   % Area of uncertainty
%       │  └ Frequency              % Freq. at which set is computed [rad/s]
%       └ EllapsedTime              % Time taken
% 
%   ----------------------------------------------------------------------------
%   Verbosity Levels:
%       0: No output
%       1: Basic logs
%       2: Complete logs with optimisation logs
% 
%   ----------------------------------------------------------------------------
%   Description:
%       Implements the method to find the elliptical uncertainty set for
%       LTI-SISO systems. 
%       
%           G(jω) = ⎰ Ĝ(jω) + ∆ ⏐ ║A(jω) ⎡Re{∆}⎤║  ≤ 1  ⎱
%                   ⎱           ⏐ ║      ⎣Im{∆}⎦║₂      ⎰
% 
%       Ĝ(jω) = Nominal Model
%       ∆     = Uncertainty
% 
%   ----------------------------------------------------------------------------
%   Copyright 2024 Vaibhav Gupta, DDMAC, EPFL (MIT License)
% 

    arguments
        frd_models (:, :, :) lti
        w (:, 1) double = []
        options.G lti = tf();
        options.solver (1, 1) string = "mosek"
        options.verbosity (1, 1) int32 {mustBeMember(options.verbosity, [0, 1, 2])} = 1
    end
    
    % In case of error, return empty struct
    model = struct();
    
    % Logger handle
    logger = @(type, varargin) utils.logger_verbose(options.verbosity, type, varargin{:});
    logger("START", mfilename);
    cleanUp = onCleanup(@() logger("STOP", mfilename));
    
    [ny, nu, n_models] = size(frd_models);

    if isempty(w)
        % If freq. points are not specified
        w = logspace(...
                log10(frd_models.Frequency(1)), ...
                log10(frd_models.Frequency(end)), ...
                500).';
        logger("DEBUG", "Using 500 logarithmically spaced frequency points for set estimation!")
    end
    nw = length(w);
    
    Ts = frd_models.Ts;
    
    %% Elliptical uncertainty set & Nominal model
    tic;
       
    if ~isempty(options.G)
        % Fixed nominal model --------------------------------------------------
        logger("INFO", "Finding elliptical uncertainty set with fixed nominal model...");
        
        % MVEE Optimisation
        localMVEE = MVEE_optimizer_centered(2, n_models, options.solver, options.verbosity);
        
        A_all = cell(ny, nu);
        % Hamdard uncertianty modelling
        for i_output = 1:ny
            for i_input = 1:nu
                % Freq. points
                p = freqresp(frd_models(i_output, i_input), w) - freqresp(options.G(i_output, i_input), w);
                P = [real(p); imag(p)];
                P = permute(P, [1, 4, 3, 2]);
                
                A = cell(nw, 1);
                for i_w = 1:nw
                    A{i_w} = localMVEE(P(:,:,i_w));
                end
                
                A_all{i_output, i_input} = A;
            end
        end
        
        % Generate a FRD model
        G = frd(freqresp(options.G, w), w, Ts);
    else
        % Find nominal model ---------------------------------------------------
        logger("INFO", "Finding elliptical uncertainty set and the nominal model...");
        
        % MVEE Optimisation
        localMVEE = MVEE_optimizer(2, n_models, options.solver, options.verbosity);
        
        A_all = cell(ny, nu);
        G_all = zeros(ny, nu, nw);
        % Hamdard uncertianty modelling
        for i_output = 1:ny
            for i_input = 1:nu
                % Freq. points
                p = freqresp(frd_models(i_output, i_input), w);
                P = [real(p); imag(p)];
                P = permute(P, [1, 4, 3, 2]);
                
                output = cell(nw, 1);
                for i_w = 1:nw
                    output{i_w} = localMVEE(P(:,:,i_w));
                end
                A = cellfun(@(x) x{1}, output, 'UniformOutput', 0);
                c = cellfun(@(x) - (x{1} \ x{2}), output, 'UniformOutput', 0);

                % Nominal model
                G_freqresp = cellfun(@(x) x(1) + 1i*x(2), c);
                
                A_all{i_output, i_input} = A;
                G_all(i_output, i_input, :) = G_freqresp;
            end
        end
        % Generate a FRD model
        G = frd(G_all, w, Ts);
    end
    ellapsed_time = toc;
    
    % Area of Uncertainty
    area = cellfun( @(A) pi * cellfun(@(x) 1 / det(x), A) , A_all, "UniformOutput", 0);
    
    logger("INFO", "Elliptical uncertainty set found in %.4f s", ellapsed_time);
    
    %% Output
    model = struct(...
        "NominalModel", G, ...
        "UncertaintySet", struct(...
            "Type"      , "elliptical", ...
            "A"         , {A_all}, ...
            "Area"      , {area}, ...
            "Frequency" , w ...
        ), ...
        "EllapsedTime", ellapsed_time ...
        );

end

%% Local Helper Functions ------------------------------------------------------
function MVEE = MVEE_optimizer(d, N, solver, verbosity)
%MVEE_OPTIMIZER Generate an optimizer object for ellipsoid fit
% 
%   MVEE = MVEE_optimizer(d, N) returns the optimizer object which takes in a
%   (d, N) matrix and returns the pair {A, b}. Here, `d` is the dimension of the
%   points to fit and `N` is the number of points.
% 
%   ----------------------------------------------------------------------------
%   Description:
%       The minimum volume enclosing ellipsoid (MVEE) of the data points P is to
%       be found. The following dual optimization problem is solved,
%               min   {- log(det(A))}
%           s.t.	||A * xi + b||_2 <= 1
%       in variables A and b with xi as the the i-th column of the matrix P. 
% 
%   ----------------------------------------------------------------------------
%   Copyright:
% 

    arguments
        d (1, 1) double
        N (1, 1) double
        solver (1, 1) string
        verbosity (1, 1) int32
    end

    P = sdpvar(d, N, 'full');
    A = sdpvar(d, d, 'symmetric');
    b = sdpvar(d, 1);

    constraints = [...
        A >= 0;
        splitapply(@(x) cone(A*x + b, 1.1), P, 1:N);
        ];

    objective =  - log(det(A));

    opt = sdpsettings(...
        'solver', char(solver), ...
        'verbose', verbosity==2 ...
        ); 
    MVEE = optimizer(constraints, objective, opt, P, {A, b});

end

function MVEE = MVEE_optimizer_centered(d, N, solver, verbosity)
%MVEE_OPTIMIZER Generate an optimizer object for centered ellipsoid fit
% 
%   MVEE = MVEE_optimizer(d, N) returns the optimizer object which takes in a
%   (d, N) matrix and returns the matrix A. Here, `d` is the dimension of the
%   points to fit and `N` is the number of points.
% 
%   ----------------------------------------------------------------------------
%   Description:
%       The minimum volume enclosing ellipsoid (MVEE), which is centered at 
%       origin, of the data points P is to be found. The following dual
%       optimization problem is solved,
%               min   {- log(det(A))}
%           s.t.	||A * xi||_2 <= 1
%       in variables A with xi as the the i-th column of the matrix P. 
% 
%   ----------------------------------------------------------------------------
%   Copyright:
% 

    arguments
        d (1, 1) double
        N (1, 1) double
        solver (1, 1) string
        verbosity (1, 1) int32
    end

    P = sdpvar(d, N, 'full');
    A = sdpvar(d, d, 'symmetric');
   
    constraints = [...
        A >= 0;
        splitapply(@(x) cone(A*x, 1), P, 1:N);
        ];

    objective =  - log(det(A));

    opt = sdpsettings(...
        'solver', char(solver), ...
        'verbose', verbosity==2 ...
        ); 
    MVEE = optimizer(constraints, objective, opt, P, A);

end
