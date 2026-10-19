classdef Fixed < datadriven.Controller.SS
    %FIXED A frozen (numeric) copy of a datadriven.Controller.SS
    %
    %   K_fixed = datadriven.Controller.Fixed(K) builds a controller object
    %   that is still of type datadriven.Controller.Controller (so it
    %   passes the "controller datadriven.Controller.Controller" argument
    %   check used throughout the Objective functions), but whose internal
    %   A, B, C, D matrices are frozen to double() of K's CURRENT values
    %   instead of being live sdpvars.
    %
    %   Use this instead of ss(K) when you need to hold a controller fixed
    %   for a stage that optimises something else (e.g. the IQC multiplier
    %   in IQC_reduced_PI) - ss(K) returns a plain MATLAB "ss" object,
    %   which does not satisfy the datadriven.Controller.Controller type
    %   and will raise "Invalid argument ... Value must be of type
    %   datadriven.Controller.Controller".
    %
    %   freqresp(K_fixed, omegas) / evalfr(K_fixed, omega) then return
    %   purely numeric X, Y (no sdpvar), which is exactly what you want
    %   when the controller is meant to be held fixed while D, E (the IQC
    %   multiplier) are optimised.
    %
    %   Example (alternating scheme):
    %       [LMIs_c, slack_c] = datadriven.Objective.IQC_reduced_PI(G, K, PI_frd, alpha, beta, omegas);
    %       optimize(LMIs_c, slack_c, opts_solver);
    %       K_fixed = datadriven.Controller.Fixed(K);   % <-- freeze K here
    %       [LMIs_pi, slack_pi, D, E] = datadriven.Objective.IQC_reduced_PI(G, K_fixed, eta, alpha, beta, omegas, ...);
    %
    methods
        function obj = Fixed(K)
            arguments
                K datadriven.Controller.SS
            end
            % Run the normal SS constructor first (this wastes a random
            % initialise() call), then immediately overwrite A, B, C, D
            % with double() of K's current values, breaking the sdpvar
            % link. Being a subclass of SS, we have protected access to
            % K's A, B, C, D even though they are not publicly readable.
            obj@datadriven.Controller.SS(K.nx, K.ny, K.nu, K.Ts, K.factorisation);
            obj.A = double(K.A);
            obj.B = double(K.B);
            obj.C = double(K.C);
            obj.D = double(K.D);
        end
    end
end
