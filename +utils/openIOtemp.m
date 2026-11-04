function [u, y] = openIOtemp(path, itmin, itmax)
%OPENIOTEMP Load the input (command) and output (measurement) of one experiment
%
%   [u, y] = utils.openIOtemp(path)              whole record
%   [u, y] = utils.openIOtemp(path, itmin, itmax) samples itmin..itmax
%                                                 (itmax = [] -> until the end)
%
%   The .mat file holds one recorder struct whose field Y(2) is the command
%   (V) and Y(4) the measurement (m/s^2). Both signals are plotted.
%
arguments (Input)
    path
    itmin = 1
    itmax = []
end

    loaddata = load(path);
    varnames = fieldnames(loaddata);
    data     = loaddata.(varnames{1}).Y;
    u = data(2).Data;
    y = data(4).Data;

    if isempty(itmax)
        itmax = numel(u);
    end
    u = u(itmin:itmax);
    y = y(itmin:itmax);

    figure(Name=path)
    plot(y); hold on; plot(u)
    xlabel('i')
    ylabel('command (V) - measurements (m/s^2)')
end
