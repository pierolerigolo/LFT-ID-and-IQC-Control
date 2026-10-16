function magplot(model,plotoptions,options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    model=''
    plotoptions.DisplayName = ''
    plotoptions.Color
    plotoptions.Linestyle
    plotoptions.Linewidth = 1
    plotoptions.HandleVisibility = 'on'
    options.Title =''
    options.iuiy=[]
    options.ZoomFreq = false
    options.Freq = []
    options.phase=false
    options.freqphacorr=150 %rad/s
    options.degphacorrobj=0
    options.ylabel=true
    options.NRMSE_var=[]
    options.xlinleg = 'on'
    options.xlabel=true
end
    hold on
    if options.phase==false
        if ~isempty(model)
            if isa(model,'tf')|isa(model,'ss')
                w_grid=linspace(2*pi*10,pi/model.Ts,700);
                model=frd(model,w_grid);
            end
            plot(model.Frequency/2/pi,mag2db(squeeze(abs(model.ResponseData))), ...
                plotoptions)
            xlim([10 1/2/model.Ts])
        end
        if options.ylabel
            ylabel('Magnitude (dB)')
        end
    else
        if ~isempty(model)
            [~,phase]=bode(model);
            phase=squeeze(phase);
            freq=squeeze(model.Frequency);
            degphacorr=phase(freq>options.freqphacorr & ...
                freq<(options.freqphacorr+5));
            shiftphase=360*floor(degphacorr(1)/360)-options.degphacorrobj;
            plot(freq/2/pi,phase-shiftphase,plotoptions)
            xlim([10 model.Frequency(end)/2/pi])
        end 
        if options.ylabel
            ylabel('Phase (deg)')
        end
    end
    xscale('log')
    if options.xlabel ; xlabel('Frequency (Hz)');end
    grid on
    if ~isempty(plotoptions.DisplayName)
        legend('Location','best')
    end
    if ~isempty(options.iuiy)
        subtitle(sprintf('From: u_%d to: y_%d',options.iuiy))
    end
    if ~isempty(options.Title)
        title(options.Title)
    end
    if ~isempty(options.Freq) & ~options.ZoomFreq
        xline(options.Freq(1)/2/pi,'DisplayName','Varying part', ...
            'HandleVisibility',options.xlinleg,'Label',sprintf('%1.f hz',options.Freq(1)/2/pi))
        xline(options.Freq(2)/2/pi,'DisplayName','freq. window', ...
            'HandleVisibility',options.xlinleg,'Label',sprintf('%1.f hz',options.Freq(2)/2/pi))
    end
    if (options.ZoomFreq)
        subtitle('Zoomed')
        xlim(options.Freq)
        legend off
    end
    if ~isempty(options.NRMSE_var)
        line=findobj(gca,'type','line');
        for i=1:5
            set(line(2*(6-i)-1),'DisplayName', ...
                sprintf('NRMSE_{var} %.1f%%',options.NRMSE_var(i)))
        end
        legend('Location','best','Visible','on','NumColumns',1)
    end
end