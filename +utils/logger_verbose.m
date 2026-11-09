function logger_verbose(verbosity, type, varargin)
% logger_verbose Helper function for logging with verbosity level checks
% 
%   logger_verbose(verbosity, type, vargin) logs based on verbosity level
%       requirements and type.
% 
%   ----------------------------------------------------------------------------
%   Copyright 2024 Vaibhav Gupta, DDMAC, EPFL (MIT License)
% 

    if verbosity > 0
        % START, STOP, and IMP Print (Prints always)
        switch type
            case "START"
                fprintf("%s\n", repmat('-', 80, 1));
                if length(varargin) == 1
                    fprintf("[%5s] Running `%s` ... \n", "IMP", varargin{:});
                elseif length(varargin) > 1
                    fprintf("[%5s] %s\n", "IMP", sprintf(varargin{:}));
                end
            case "STOP"
                if length(varargin) == 1
                    fprintf("[%5s] `%s` finished. \n", "IMP", varargin{:});
                elseif length(varargin) > 1
                    fprintf("[%5s] %s\n", "IMP", sprintf(varargin{:}));
                end
                fprintf("%s\n", repmat('-', 80, 1));
            case "IMP"
                fprintf("[%5s] %s\n", "IMP", sprintf(varargin{:}));
        end
    end
    
    % Print nothing at verbosity level 0
    switch verbosity
        case 0
            % Only print errors
            if strcmpi(type, "ERROR")
                fprintf("[%5s] %s\n", type, sprintf(varargin{:}));
            end
        case 1
            % Print errors & infos
            if any(strcmpi(type, ["ERROR", "INFO", "WARN"]))
                fprintf("[%5s] %s\n", type, sprintf(varargin{:}));
            end
        case 2
            % Print everything
            fprintf("[%5s] %s\n", type, sprintf(varargin{:}));
    end
end