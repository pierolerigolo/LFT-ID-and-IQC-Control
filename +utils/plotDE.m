function plotDE(D_val, E_val, W, eta, name, file)
%PLOTDE Mass multiplier [eta^2 D, E; E', -D] of the stage-3 bisection across frequency
%
%   utils.plotDE(D_val, E_val, W, eta)
%   utils.plotDE(D_val, E_val, W, eta, name, file)   title and pdf file (optional)
%
%   D_val, E_val : cells (one entry per frequency of W, rad/s) as returned by the
%                  bisection (saved in the controller files); scalar mass block here.
%   Plots D(w) and eta^2 D(w) (D > 0), and E(w) (E = -E': purely imaginary for a scalar
%   block, so Im E is shown, with |Re E| as a check). A constant D (= its warm start)
%   and E = 0 mean stage 3 did not run or did not move the multiplier.
%
    if nargin < 5, name = ""; end
    if nargin < 6, file = ""; end
    if ~iscell(D_val) || isempty(D_val)
        warning('plotDE:empty', 'No D(w), E(w) to plot.');   return
    end
    f  = W(:)/2/pi;
    Dw = cellfun(@(x) real(x(1)), D_val(:));
    Ew = cellfun(@(x) x(1), E_val(:));
    figure('Name', 'Mass multiplier D, E')
    tl = tiledlayout(2, 1, 'TileSpacing', 'compact');
    nexttile
    semilogx(f, Dw, 'b', f, eta^2*Dw, 'b--'); grid on; xlim([f(1) f(end)])
    ylabel('D(\omega)'); legend('D', sprintf('\\eta^2 D  (\\eta = %.3g)', eta), 'Location', 'best')
    nexttile
    semilogx(f, imag(Ew), 'r', f, abs(real(Ew)), 'k:'); grid on; xlim([f(1) f(end)])
    ylabel('E(\omega)'); xlabel('Frequency (Hz)'); legend('Im E', '|Re E| (should be 0)', 'Location', 'best')
    title(tl, "Mass multiplier [\eta^2 D, E; E^*, -D]" + ifelse(strlength(name) > 0, " - " + name, ""))
    fprintf('D: min %.3g, max %.3g, median %.3g | Im E: min %.3g, max %.3g | max |Re E| %.2g\n', ...
        min(Dw), max(Dw), median(Dw), min(imag(Ew)), max(imag(Ew)), max(abs(real(Ew))));
    if strlength(file) > 0
        exportgraphics(gcf, file, 'ContentType', 'vector');
    end
end

function s = ifelse(c, a, b)
    if c, s = a; else, s = b; end
end
