function [h, h_ell] = plot_ellipse(c, A, opts)
%PLOT_ELLIPSE Plot ellipse and its center
% 
%   PLOT_ELLIPSE(c, A) plot ellipse with center `c` and shape matrix `A`.
% 
%   ----------------------------------------------------------------------------
%   Description:
%       Plots ellipse of form
%           ║A ⎡x - Re{c}⎤║  ≤ 1
%           ║  ⎣y - Im{c}⎦║₂    
% 
%   ----------------------------------------------------------------------------
%   Copyright 2024 Vaibhav Gupta, DDMAC, EPFL (MIT License)
% 

    arguments
        c (1, 1) double
        A (2, 2) double
        opts.Color = "#050505"
        opts.Marker = "."
        opts.MarkerSize = 10
        opts.LineStyle = "--"
        opts.LineWidth = 2
    end


    [V, D] = eig(A);
    phi = atan2(V(2, 1), V(1, 1));
    a = 1/D(1, 1);
    b = 1/D(2, 2);
    
    t = linspace(0,2*pi)';
    X = a*cos(t) + 1i * (b*sin(t));
    C = X * exp(1i*phi) + c;

    h = plot(real(c), imag(c), ...
        "LineStyle", "none", ...
        "Marker", opts.Marker, ...
        "MarkerSize",  opts.MarkerSize, ...
        "Color",  opts.Color, ...
        "MarkerFaceColor",  opts.Color);
    h_ell = plot(real(C), imag(C), ...
        "LineStyle", opts.LineStyle, ...
        "LineWidth", opts.LineWidth, ...
        "Color", opts.Color);
end

