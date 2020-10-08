function traces_overview(showalltr, showevents)

% Declare globally shared variables
global WHTRCFIG POSITIONTRCFIG FONTSIZE CLIMRAW TRCYLABEL TRCXLABEL TRCPLOTCOL1...
    TRCALPHA THRESHLW THRESHALPHA EVPLOTCOL PLOTPEAKS THRESHOLD PEAKTHRESHOLD SMTHWIN...
    figINFO roiINFO traceINFO FTIME SAVEPARENT FNAME...
    TRACEID TRACEDATA ROILIST RANGE TOTP evINFO LW1 LW2


% Check input
if ~isa(showevents, 'logical'), showevents = false; end

% Traces to be displayed
if showalltr, seltrace = 1:size(TRACEDATA,2);
else, seltrace = TRACEID;
end
disptr = TRACEDATA(RANGE(1):RANGE(2), seltrace);
dispid = ROILIST(seltrace);
numtr = size(disptr,2);
timeptstot = (1:TOTP).*FTIME;
timepts = timeptstot(RANGE(1):RANGE(2));
if SMTHWIN ~= 0
    showsmth = true;
    disptrsmth = zeros(size(disptr));
    for iTr = 1:size(disptr,2), disptrsmth(:,iTr) = smooth_data(disptr(:,iTr),SMTHWIN); end
else
    showsmth = false;
end

%% Display parameters
if showevents, spcol = 7; else, spcol = 1; end
if numtr > 6, sprow = 6; else, sprow = numtr; end
hspace = 40;
vspace = 40;
if showevents, whsp = [floor((WHTRCFIG(1)-3*hspace)/spcol) (WHTRCFIG(2)-(sprow+2)*vspace)/sprow];
else, whsp = [WHTRCFIG(1)-2*hspace (WHTRCFIG(2)-(sprow+1)*vspace)/sprow];
end
whan = [120 50];
whsp = whsp./WHTRCFIG; hspace = hspace/WHTRCFIG(1); vspace = vspace/WHTRCFIG(2); whan=whan./WHTRCFIG;
sppos1 = zeros(sprow,4);
sppos2 = zeros(sprow,4);
ypos = 0;
for iR = 1:sprow
    if showevents, sppos1(iR,:) = [hspace ypos+vspace whsp(1)*4 whsp(2)];
    else, sppos1(iR,:) = [hspace ypos+vspace whsp(1) whsp(2)];
    end
    if showevents, sppos2(iR,:) = [sppos1(iR,1)+sppos1(iR,3)+hspace ypos+vspace whsp(1) whsp(2)]; end
    ypos = ypos+vspace+whsp(2);
end
sppos1 = flipud(sppos1); sppos2 = flipud(sppos2);

figure('Position', POSITIONTRCFIG, 'Name', 'Traces 1');
iFig = 1; iSP = 1;
for iTr = 1:numtr
    sp1 = subplot('Position', sppos1(iSP,:),'Parent', gcf);
    t1 = plot(timepts, disptr(:,iTr),'linewidth', LW1, 'color', TRCPLOTCOL1); t1.Color(4) = TRCALPHA;
    if showsmth, hold on; t2 = plot(timepts, disptrsmth(:,iTr),'linewidth', LW2, 'color', 'black'); end
    title(dispid{iTr});
    if mod(iTr,sprow) == 0, xlabel(TRCXLABEL); end
    ylabel(TRCYLABEL);
    set(gca, 'FONTSIZE', FONTSIZE);
    xlim([timepts(1), timepts(end)]);
    if showevents
        hold on, ty = yline(eventINFO(iTr).threshold, 'color', [0 0 0], 'linewidth', THRESHLW); ty.Color(4) = THRESHALPHA;
        ty = yline(eventINFO(iTr).peakthreshold, 'color', [0 0 0], 'linewidth', THRESHLW/2); ty.Color(4) = THRESHALPHA;
        if PLOTPEAKS, evidx = timeptstot(eventINFO(iTr).peakidx);
        else, evidx = timeptstot(eventINFO(iTr).crossidx);
        end
        
        for iE = 1:numel(evidx), hold on; xline(evidx(iE), 'Linewidth',lw3, 'Color', EVPLOTCOL); end
        annotxt = {sprintf('Average | SD: %s  |  %s', num2str(eventsINFO(iTr).average,3), num2str(eventsINFO(iTr).sd,3)),...
            sprintf('Threshold: %s', num2str(eventsINFO(iTr).threshold,3)),...
            sprintf('Threshold peak: %s', num2str(eventsINFO(iTr).peakthreshold,3)),...
            sprintf('Event type: %s', eventsINFO(iTr).eventtype),...
            sprintf('Av. Inter-event-interval: %s s', num2str(eventsINFO(iTr).aviei,3)), ...
            sprintf('CV IEI: %s', num2str(eventsINFO(iTr).cviei,3)),...
            sprintf('Av. Amplitude: %s', num2str(eventsINFO(iTr).avamp,3)),...
            sprintf('Av. Eventrate: %s Hz', num2str(eventsINFO(iTr).eventrate,3))};
        hold off;
        
        sp2 =  subplot('Position', sppos2(iSP,:),'Parent', gcf); axis off;
        evinfo = text(0.02,0.5, annotxt, 'FONTSIZE', FONTSIZE, 'backgroundcolor', [1 1 1], 'edgecolor', [0 0 0]);
    end
    if mod(iTr,sprow) == 0 && iTr < numtr
        iFig = iFig+1; iSP = 1;
        figure('Position', POSITIONTRCFIG, 'Name', sprintf('Traces %i',iFig));
    else
        iSP = iSP+1;
    end
end

end
