function update_evinfo(traceidx, keepev)

%% Declare globally shared variables
global PLOTPEAKS FTIME evINFO

%% Update information
crossidx = evINFO(traceidx).crossidx;
peakidx = evINFO(traceidx).peakidx;

delcross = crossidx(~keepev);
delpeaks = peakidx(~keepev);

evINFO(traceidx).binarycross(delcross) = 0;
evINFO(traceidx).binarypeaks(delpeaks) = 0;
evINFO(traceidx).crossidx = evINFO(traceidx).crossidx(keepev);
evINFO(traceidx).peakidx = evINFO(traceidx).peakidx(keepev);
evINFO(traceidx).amps = evINFO(traceidx).amps(keepev);
evINFO(traceidx).amps_dFoF = evINFO(traceidx).amps_dFoF(keepev);
evINFO(traceidx).baselinevalues = evINFO(traceidx).baselinevalues(keepev);
evINFO(traceidx).risetimes = evINFO(traceidx).risetimes(keepev);
evINFO(traceidx).decaytimes = evINFO(traceidx).decaytimes(keepev);

%% Recalculate
if PLOTPEAKS, evINFO(traceidx).ieis = diff(evINFO(traceidx).peakidx).*FTIME;
else, evINFO(traceidx).ieis = diff(evINFO(traceidx).crossidx).*FTIME;
end
evINFO(traceidx).cviei = std(evINFO(traceidx).ieis)/mean(evINFO(traceidx).ieis);
%     eventrate = numel(crossidx)/(numel(tmptrace)*FTIME); % Eventrate calculated from number of events [Hz]
evINFO(traceidx).eventrate = 1/median(evINFO(traceidx).ieis); % Eventrate calculated from intereventintervals [Hz]

evINFO(traceidx).avamp = mean(evINFO(traceidx).amps);
evINFO(traceidx).avamp_dFoF = mean(evINFO(traceidx).amps_dFoF);
evINFO(traceidx).aviei = mean(evINFO(traceidx).ieis);
evINFO(traceidx).accepted = 'accepted';
end