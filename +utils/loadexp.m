function [y, d, u] = loadexp(path)
%LOADEXP Signals of one closed-loop experiment (recorder .mat file)
%
%   [y, d, u] = utils.loadexp(path)
%     y : payload acceleration     (Y(10), m/s^2)
%     d : disturbance command      (Y(13), V)
%     u : controller output        (Y(1),  zero while the controller is off)
%
    loaddata = load(path);
    varnames = fieldnames(loaddata);
    Y = loaddata.(varnames{1}).Y;
    y = Y(10).Data(:);
    d = Y(13).Data(:);
    u = Y(1).Data(:);
end
