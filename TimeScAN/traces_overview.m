function traces_overview(showalltr, showevents)

% Declare globally shared variables
global WHTRCFIG POSITIONTRCFIG FONTSIZE TRCYLABEL TRCXLABEL TRCPLOTCOL1...
    TRCALPHA THRESHLW THRESHALPHA EVPLOTCOL PLOTPEAKS SMTHWIN...
    TRACEID TRACEDATA ROILIST RANGE evINFO LW1 LW2 FTIMEVEC XTYPE


% Check input
if ~isa(showevents, 'logical'), showevents = false; end

% Traces to be displayed
if showalltr, seltrace = 1:size(TRACEDATA,2);
else, seltrace = TRACEID;
end
disptr = TRACEDATA(RANGE(1):RANGE(2), seltrace);
dispid = ROILIST(seltrace);
numtr = size(disptr,2);
if strcmp(XTYPE,'s'), timepts = FTIMEVEC(RANGE(1):RANGE(2));
else, timepts = RANGE(1):RANGE(2);
end
if SMTHWIN ~= 0
    showsmth = true;
    disptrsmth = zeros(size(disptr));
    for iTr = 1:size(disptr,2), disptrsmth(:,iTr) = smooth_data(disptr(:,iTr),SMTHWIN); end
else
    showsmth = false;
end

%% Display parameters
if showevents
    spcol = 7;
    if numtr > 4, sprow = 4; else, sprow = numtr; end
else
    spcol = 1;
    if numtr > 6, sprow = 6; else, sprow = numtr; end
end
hspace = 40;
vspace = 50;
if showevents, whsp = [floor((WHTRCFIG(1)-3*hspace)/spcol) (WHTRCFIG(2)-(sprow+2)*vspace)/sprow];
else, whsp = [WHTRCFIG(1)-2*hspace (WHTRCFIG(2)-(sprow+1)*vspace)/sprow];
end
whan = [120 50];
whsp = whsp./WHTRCFIG; hspace = hspace/WHTRCFIG(1); vspace = vspace/WHTRCFIG(2); whan=whan./WHTRCFIG;
sppos1 = zeros(sprow,4);
sppos2 = zeros(sprow,4);
ypos = 0;
for iR = 1:sprow
    if showevents, sppos1(iR,:) = [hspace ypos+vspace whsp(1)*6 whsp(2)];
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
    t1 = plot(timepts, disptr(:,iTr),'linewidth', LW1, 'color', TRCPLOTCOL1);
    if showsmth
        t1.Color(4) = TRCALPHA;
        hold on; t2 = plot(timepts, disptrsmth(:,iTr),'linewidth', LW2, 'color', 'black');
    end
    title(dispid{iTr});
    if mod(iTr,sprow) == 0 || iTr == numtr, xlabel(TRCXLABEL); end
    ylabel(TRCYLABEL);
    set(gca, 'FONTSIZE', FONTSIZE);
    xlim([timepts(1), timepts(end)]);
    if showevents && strcmp(evINFO(seltrace(iTr)).accepted,'accepted')
        if strcmp(evINFO(seltrace(iTr)).thresholdtype, 'sd')
            hold on, ty = yline(evINFO(seltrace(iTr)).threshold, 'color', [0 0 0], 'linewidth', THRESHLW);
            ty.Color(4) = THRESHALPHA;
            ty = yline(evINFO(seltrace(iTr)).peakthreshold, 'color', [0 0 0], 'linewidth', THRESHLW/2);
            ty.Color(4) = THRESHALPHA;
        end
        if PLOTPEAKS, evidx = timepts(evINFO(seltrace(iTr)).peakidx);
        else, evidx = timepts(evINFO(seltrace(iTr)).onsetidx);
        end
        
        for iE = 1:numel(evidx), hold on; xline(evidx(iE), 'Linewidth',LW1, 'Color', EVPLOTCOL); end
        annotxt = {sprintf('Average | SD: %s  |  %s', num2str(evINFO(seltrace(iTr)).average,3), num2str(evINFO(seltrace(iTr)).sd,3)),...
            sprintf('Threshold: %s', num2str(evINFO(seltrace(iTr)).threshold,3)),...
            sprintf('Threshold peak: %s', num2str(evINFO(seltrace(iTr)).peakthreshold,3)),...
            sprintf('Threshold type: %s', evINFO(seltrace(iTr)).thresholdtype),...
            sprintf('Av. Inter-event-interval: %s s', num2str(evINFO(seltrace(iTr)).aviei,3)), ...
            sprintf('CV IEI: %s', num2str(evINFO(seltrace(iTr)).cviei,3)),...
            sprintf('Av. Amp. (dFoF): %s (%s)', num2str(evINFO(seltrace(iTr)).avamp,3), num2str(evINFO(seltrace(iTr)).avamp_dFoF,3)),...
            sprintf('Av. Eventrate: %s Hz', num2str(evINFO(seltrace(iTr)).eventrate,3))};
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
