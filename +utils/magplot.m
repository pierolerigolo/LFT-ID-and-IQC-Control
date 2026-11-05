function magplot(model, plotoptions)
%MAGPLOT Magnitude (dB) vs frequency (Hz, log scale) on the current axes
%
%   utils.magplot(sys)                     sys: frd, or tf/ss (evaluated
%                                          from 10 Hz to Nyquist)
%   utils.magplot(sys, 'Color', c, 'DisplayName', name, ...)
%                                          line properties passed to plot
%
arguments (Input)
    model
    plotoptions.DisplayName = ''
    plotoptions.Color
    plotoptions.LineStyle
    plotoptions.LineWidth = 1
    plotoptions.HandleVisibility = 'on'
end
    hold on
    if isa(model, 'tf') || isa(model, 'ss')
        model = frd(model, linspace(2*pi*1, pi/model.Ts, 700));
    end
    plot(model.Frequency/2/pi, mag2db(squeeze(abs(model.ResponseData))), plotoptions)
    xlim([10 1/2/model.Ts])
    ylabel('Magnitude (dB)')
    xlabel('Frequency (Hz)')
    xscale('log')
    grid on
    if ~isempty(plotoptions.DisplayName)
        legend('Location', 'best')
    end
end
