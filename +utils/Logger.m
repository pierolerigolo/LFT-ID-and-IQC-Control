classdef Logger < handle
    %LOGGER Class to log data-driven approach 
    
    properties
        fileID = 1;
        n_props
        syntax
        err_syntax
        stage_syntax
        divider
        title
        
        iter_syntax = "| % 4d |";
        iter_printed = 0;
        
        printf = @(fid, varargin) arrayfun(@(f) fprintf(f,  varargin{:}), fid, "UniformOutput", 0);
    end
    
    methods
        function obj = Logger(fileID, prop_name, prop_format)
            arguments
                fileID
            end
            arguments (Repeating)
                prop_name string
                prop_format string
            end
            
            obj.fileID = fileID;
            obj.n_props = length(prop_name);
            obj.syntax = sprintf(" %s |", prop_format{:}) + "\n";
            
            tmp = @(x) regexp(x, "(?:%(?:\d\$)?[-+ 0#]?)(\d+)[^a-z$.]?", "tokens");
            field_widths = cellfun(@(y) cellfun(@double, tmp(y)), prop_format);
            
            obj.divider = "|------|" + sprintf("%s|", arrayfun( ...
                @(y) string(repelem('-', y+2)), ...
                field_widths)) + "\n";
            
            obj.title = "| Iter |" + ...
                sprintf(sprintf(" %%-%ds |", field_widths), prop_name{:}) + ...
                "\n";
            
            val = sum(field_widths) + length(field_widths)*3 - 3;
            obj.err_syntax = " %-" + string(val) + "s |\n";
            obj.stage_syntax = "|======| %-" + string(val) + "s |\n";
        end
        
        function init(obj)
            obj.printf(obj.fileID, obj.divider);
            obj.printf(obj.fileID, obj.title);
            obj.printf(obj.fileID, obj.divider);
        end
        
        function cleanup(obj)
            obj.printf(obj.fileID, obj.divider);
        end
        
        function log_iter(obj, iter)
            if obj.iter_printed ~= 0
                obj.printf(obj.fileID, repmat('\b', 1, obj.iter_printed));
            end
            obj.printf(obj.fileID, obj.iter_syntax, iter);
            obj.iter_printed = 1;
        end
        
        function log(obj, varargin)
            if obj.iter_printed == 0
                obj.log_iter(0);
            end
            if length(varargin) ~= obj.n_props
                warning("Incorrect number of data");
                return;
            end
            obj.printf(obj.fileID, obj.syntax, varargin{:});  
            obj.iter_printed = 0;
        end
        
        function log_err(obj, varargin)
            if obj.iter_printed == 0
                obj.log_iter(0);
            end
            if length(varargin) ~= 1
                warning("Only single error message is allowed.");
                return;
            end
            obj.printf(obj.fileID, obj.err_syntax, varargin{:});
            obj.iter_printed = 0;
        end
        
        function log_stage(obj, varargin)
            if length(varargin) ~= 1
                warning("Only single stage message is allowed.");
                return;
            end
            
            if obj.iter_printed ~= 0
                obj.printf(obj.fileID, repmat('\b', 1, obj.iter_printed));
                obj.iter_printed = 0;
            end
            obj.printf(obj.fileID, obj.stage_syntax, varargin{:});
        end
    end
end

