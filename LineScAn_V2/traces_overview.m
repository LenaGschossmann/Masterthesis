function traces_overview()

% Declare globally shared variables
global WHTRCFIG POSITIONTRCFIG FONTSIZE PLOTRANGE TRCYLABEL TRCXLABEL TRCPLOTCOL1...
    TRCALPHA SCMAP EVDATA roiINFO traceINFO FTIME IMHW

[traceidx, roiidx] = access_trc();
numrois = numel(roiidx);
trcxrange = 1*FTIME :FTIME: IMHW(1)*FTIME;
range = PLOTRANGE(1):PLOTRANGE(2);
% cscm = figINFO(figidx).cscmap;
% cscl = figINFO(figidx).csclimits;

if ~isempty(traceidx)
    % Display parameters
    spcol = 12; sprow = 3;
    hspace = 20;
    vspace = 30;
    whbut = [100 25];
    whsp = [floor((WHTRCFIG(1)-7*hspace)/spcol) floor((WHTRCFIG(2)-(sprow-1+4)*vspace-whbut(2))/sprow)];
    whrbg = [whsp(1)*2 round(whsp(2)/2)];
    whrb = [round(0.8*whrbg(1)) round(whrbg(2)/3)];
    whan = [120 50];
    rbbotpos = [0.15*whrbg(1) round(whrb(2)*0.5) whrb];
    rbtoppos = [0.15*whrbg(1) round(whrb(2)*1.5) whrb];
    whsp = whsp./WHTRCFIG; whrbg = whrbg./WHTRCFIG; hspace = hspace/WHTRCFIG(1); vspace = vspace/WHTRCFIG(2); whbut=whbut./WHTRCFIG; whrb=whrb./WHTRCFIG; whan=whan./WHTRCFIG;
    spposleft = zeros(sprow,4);
    spposmid = zeros(sprow,4);
    spposright = zeros(sprow,4);
    bgpos = zeros(sprow,4);
    rbsavetrcpos = zeros(sprow,4);
    rbshowallpos = zeros(sprow,4);
    ypos = vspace*2+whbut(2);
    for iR = 1:sprow
        spposleft(iR,:) = [hspace ypos+vspace whsp(1)*2 whsp(2)];
        spposmid(iR,:) = [spposleft(iR,1)+spposleft(iR,3)+2.5*hspace ypos+vspace whsp(1)*8 whsp(2)];
        %         spposright(iR,:) = [spposmid(iR,1)+spposmid(iR,3)+hspace ypos+vspace whsp(1)*2 whsp(2)];
        bgpos(iR,:) = [spposmid(iR,1)+spposmid(iR,3)+2*hspace ypos+vspace+whsp(2)/2 whrbg];
        rbsavetrcpos(iR,:) = [spposmid(iR,1)+spposmid(iR,3)+2*hspace ypos+vspace whrb(1)*1.25 whrb(2)];
        rbshowallpos(iR,:) = [rbsavetrcpos(iR,1) rbsavetrcpos(iR,2)+whrb(2) whrb(1)*1.25 whrb(2)];
        ypos = ypos+vspace+whsp(2);
    end
    spposleft = flipud(spposleft); spposmid = flipud(spposmid); spposright = flipud(spposright); bgpos = flipud(bgpos); rbsavetrcpos = flipud(rbsavetrcpos); rbshowallpos = flipud(rbshowallpos);
    
    figure('Position', POSITIONTRCFIG, 'Name', 'Traces 1');

    iFig = 1; iSP = 1;
    for iRoi = 1:numrois
        tmproi = roiidx(iRoi);
        sp1 = subplot('Position', spposleft(iSP,:),'Parent', gcf);
        colormap(SCMAP);
        snap = traceINFO(traceidx(iRoi)).plotmarked{1};
        imagesc(snap, [0 max(snap,[],'all')]);
        axis off;
        
        sp2 = subplot('Position', spposmid(iSP,:),'Parent', gcf);
        t1 = plot(trcxrange(range), traceINFO(traceidx(iRoi)).roi_av{1}(range),'linewidth', 1, 'color', TRCPLOTCOL1); t1.Color(4) = TRCALPHA;
        if size(EVDATA,1) >= tmproi && ~(isempty(EVDATA{tmproi,6}) && isempty(EVDATA{tmproi,7}))
            save_idx = EVDATA{tmproi,7}; save_idx = save_idx(save_idx > PLOTRANGE(1) & save_idx < PLOTRANGE(2));
            save_pts = (save_idx-PLOTRANGE(1)+1) .* FTIME;
            rev_idx = EVDATA{tmproi,6}; rev_idx = rev_idx(rev_idx > PLOTRANGE(1) & rev_idx < PLOTRANGE(2));
            rev_pts = (rev_idx-PLOTRANGE(1)+1) .* FTIME;
            hold on;
            for iEv = 1:numel(save_pts), xline(save_pts(iEv), 'Color',[0 0 0], 'LineWidth', 1.2, 'Alpha',.4); end
            for iEv = 1:numel(rev_pts), xline(rev_pts(iEv), 'r', 'LineWidth', 1.2, 'Alpha', .4); end
            hold off;
        end
        
        %         for iE = 1:numel(evidx), hold on; xline(evidx(iE), 'Linewidth',2, 'Color', EVPLOTCOL); end
        %         annotxt = {sprintf('Average | SD: %s  |  %s', num2str(traceINFO(traceidx(iRoi)).events.average,3), num2str(traceINFO(traceidx(iRoi)).events.sd,3)),...
        %             sprintf('Threshold: %s', num2str(traceINFO(traceidx(iRoi)).events.threshold,3)),...
        %             sprintf('Threshold peak: %s', num2str(traceINFO(traceidx(iRoi)).events.peakthreshold,3)),...
        %             sprintf('Event type: %s', traceINFO(traceidx(iRoi)).events.eventtype),...
        %             sprintf('Av. Inter-event-interval: %s s', num2str(traceINFO(traceidx(iRoi)).events.aviei,3)), ...
        %             sprintf('CV IEI: %s', num2str(traceINFO(traceidx(iRoi)).events.cviei,3)),...
        %             sprintf('Av. Amplitude: %s', num2str(traceINFO(traceidx(iRoi)).events.avamp,3)),...
        %             sprintf('Av. Eventrate: %s Hz', num2str(traceINFO(traceidx(iRoi)).events.eventrate,3))};
        hold off;
        ylabel(TRCYLABEL); set(gca, 'FONTSIZE', FONTSIZE);
        xlabel(TRCXLABEL); xlim([trcxrange(PLOTRANGE(1)) trcxrange(PLOTRANGE(2))]);
        
        %         sp3 =  subplot('Position', spposright(iSP,:),'Parent', gcf); axis off;
        %         evinfo = text(0.02,0.5, annotxt, 'FONTSIZE', FONTSIZE, 'backgroundcolor', [1 1 1], 'edgecolor', [0 0 0]);
        bg = uibuttongroup('parent', gcf, 'visible', 'off', 'position', bgpos(iSP,:), 'SelectionChangedFcn', {@cb_traceswitchrb, sp2, traceidx(iRoi)});
        showalltrcrb =  uicontrol('parent', gcf,'style', 'radiobutton',  'unit', 'normal','position', rbshowallpos(iSP,:),'string', '  Show full trace', 'Value', 0 ,'FONTSIZE', FONTSIZE, 'callback', {@cb_wholetrcrb, sp2, traceidx(iRoi)});
        botrb = uicontrol('parent', bg,'style', 'radiobutton', 'position', rbbotpos,'string', '  DeltaF/F', 'Value', 0 ,'FONTSIZE', FONTSIZE, 'handlevisibility', 'off');
        toprb = uicontrol('parent', bg, 'style', 'radiobutton', 'position', rbtoppos,'string', '  Average', 'Value', 1, 'FONTSIZE', FONTSIZE, 'handlevisibility', 'off');
        bg.Visible = 'on';
        if mod(iRoi,sprow) == 0 && iRoi < numrois
            iFig = iFig+1; iSP = 1;
            figure('Position', POSITIONTRCFIG, 'Name', sprintf('Traces %i',iFig));
        else
            iSP = iSP+1;
        end
    end
end

%% Local callbacks
    function cb_traceswitchrb(~,evdat, sp2, trace)
        traceINFO(trace).currmode = get(evdat.NewValue, 'String');
        update_trc(sp2, trace, traceINFO(trace).currmode, traceINFO(trace).showtot);
    end

    function cb_wholetrcrb(hObj, ~, sp2, trace)
        traceINFO(trace).showtot = logical(hObj.Value);
        update_trc(sp2, trace, traceINFO(trace).currmode, traceINFO(trace).showtot);
    end

end
