function fitpercent = NRMSE(ymodel,ymeas,options)
arguments
    ymodel 
    ymeas 
    options.inFofFreq = false
    options.absDB = false
    options.freqrange=[]
    options.magnitudeonly=false
end
% Normalized Root Mean Squared Error (NRMSE)
%compare(ymeas,ymodel)
if isa(ymeas,'frd')&isa(ymodel,'frd')
    if ~isempty(options.freqrange)
        fgrid=ymeas.Frequency;
        freqmin=options.freqrange(1);
        freqmax=options.freqrange(2);
        ymeas=interp(ymeas,fgrid((fgrid>freqmin)&(fgrid<freqmax)));
        ymodel=interp(ymodel,fgrid((fgrid>freqmin)&(fgrid<freqmax)));
    end
    ymeas=squeeze(ymeas.ResponseData);
    ymodel=squeeze(ymodel.ResponseData);

end

if isa(ymodel,'tf')
    ymodel=frd(ymodel,ymeas.frequency);
    ymeas=squeeze(ymeas.ResponseData);
    ymodel=squeeze(ymodel.ResponseData);    
end

nonnanid=(~isnan(ymeas))&(~isnan(ymodel));
ymeas=ymeas(nonnanid);
ymodel=ymodel(nonnanid);
      temp=ymodel;
    ymodel=ymeas;
    ymeas=temp;
if options.inFofFreq
    if options.absDB
        fitpercent=mag2db(abs(ymeas)) - mag2db(abs(ymodel));
    else
        fitpercent=100*(1-abs(ymeas-ymodel))./(abs(ymeas-mean(ymeas)));
    end
elseif options.magnitudeonly
    fitpercent=100*(1-norm(abs(ymodel)-abs(ymeas))/norm(abs(ymodel)-mean(abs(ymeas))));
else

        fitpercent=100*(1-norm(ymodel-ymeas)/norm(ymodel-mean(ymeas)));
end
end