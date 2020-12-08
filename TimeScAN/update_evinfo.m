function update_evinfo(traceidx, keepev)

%% Declare globally shared variables
global PLOTPEAKS FTIME evINFO PREPOINTS POSTPOINTS BASEPOINTS

%% Update information
onsetidx = evINFO(traceidx).onsetidx;
peakidx = evINFO(traceidx).peakidx;

delonset = onsetidx(~keepev);
delpeaks = peakidx(~keepev);

evINFO(traceidx).binaryonset(delonset) = 0;
evINFO(traceidx).binarypeaks(delpeaks) = 0;
evINFO(traceidx).onsetidx = evINFO(traceidx).onsetidx(keepev);
evINFO(traceidx).peakidx = evINFO(traceidx).peakidx(keepev);
evINFO(traceidx).amps = evINFO(traceidx).amps(keepev);
evINFO(traceidx).amps_dFoF = evINFO(traceidx).amps_dFoF(keepev);
evINFO(traceidx).baselinevalues = evINFO(traceidx).baselinevalues(keepev);
evINFO(traceidx).risetimes = evINFO(traceidx).risetimes(keepev);
evINFO(traceidx).decaytimes = evINFO(traceidx).decaytimes(keepev);
evINFO(traceidx).preframes = PREPOINTS;
evINFO(traceidx).postframes = POSTPOINTS;
evINFO(traceidx).baseframes = BASEPOINTS;

%% Recalculate
if numel(evINFO(traceidx).onsetidx) > 1
    if PLOTPEAKS, evINFO(traceidx).ieis = diff(evINFO(traceidx).peakidx).*FTIME;
    else, evINFO(traceidx).ieis = diff(evINFO(traceidx).onsetidx).*FTIME;
    end
else
    evINFO(traceidx).ieis = NaN;
end

evINFO(traceidx).cviei = std(evINFO(traceidx).ieis)/mean(evINFO(traceidx).ieis);
%     eventrate = numel(onsetidx)/(numel(tmptrace)*FTIME); % Eventrate calculated from number of events [Hz]
evINFO(traceidx).eventrate = 1/median(evINFO(traceidx).ieis); % Eventrate calculated from intereventintervals [Hz]

evINFO(traceidx).avamp = mean(evINFO(traceidx).amps);
evINFO(traceidx).avamp_dFoF = mean(evINFO(traceidx).amps_dFoF);
evINFO(traceidx).aviei = mean(evINFO(traceidx).ieis);
evINFO(traceidx).accepted = 'accepted';


end