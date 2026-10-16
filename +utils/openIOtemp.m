function [u,y] = openIOtemp(path,itmin,itmax,opt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    path
    itmin=1
    itmax=[]
    opt.oldrec=false
end

arguments (Output)
    u
    y
end

    %loaddata=load("earlytest data/"+path+".mat");
    loaddata=load(path)
    varnames=fieldnames(loaddata);
    rawdata=loaddata.(varnames{1});
    data=rawdata.Y;
    if opt.oldrec
        u=data(1).Data;
        y=data(2).Data;
    else
        u=data(2).Data;
        y=data(4).Data;
    end
    if isempty(itmax)
        u=u(itmin:end);
        y=y(itmin:end);
    else
        u=u(itmin:itmax);
        y=y(itmin:itmax);
    end
    figure(Name=path)
    plot(y)
    hold on
    plot(u)
    xlabel('i')
    ylabel('command (V) - measurements (m/s^2)')
end