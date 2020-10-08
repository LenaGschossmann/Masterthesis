function ev_detection_overview(detectall)

% Declare globally shared variables
global SCRSZ WHFIG WHTRCFIG POSITIONFIG POSITIONTRCFIG POSITIONROISELECT FONTSIZE SCMAPDDITEMS...
    SCMAP WINSZDDITEMS WINSZ PLOTRANGE CLIMRAW CLIMUI TRCYLABEL TRCYLABELDFOF TRCXLABEL TRCPLOTCOL1...
    TRCPLOTCOL2 TRCALPHA THRESHLW THRESHALPHA EVPLOTCOL PLOTPEAKS THRESHOLDTYPE THRESHOLD PEAKTHRESHOLD SMTHWIN...
    figINFO roiINFO traceINFO ROICNTID FIGCOUNTER CURRFILE...
    SAVEPARENT IMHW FTIME COMPOSITE2D SAVEPATH FULLFILENAMES NUMFILES FNAME IMMETA...
    TRACEDATA ROILIST RANGE TOTP TRACEID ISBGSUB WHREVFIG evINFO LW1 LW2 BASEPOINTS PREFRAMES POSTFRAMES

% Initialize variables
currdisp = 1;

% Traces to be displayed
if detectall, seltrace = 1:size(TRACEDATA,2);
else, seltrace = TRACEID;
end
selectedtr = TRACEDATA(RANGE(1):RANGE(2), seltrace);
detectid = ROILIST(seltrace);
numtr = size(selectedtr,2);
timeptstot = (1:TOTP).*FTIME;
timepts = timeptstot(RANGE(1):RANGE(2));
if SMTHWIN ~= 0
    detecttr = zeros(size(selectedtr));
    for iTr = 1:size(selectedtr,2), detecttr(:,iTr) = smooth_data(selectedtr(:,iTr),SMTHWIN); end
else
    detecttr = selectedtr;
end

% List of traces processed with current settings
idxprocessed = find([evINFO(:).rangein] >= RANGE(1) & [evINFO(:).rangeout] <= RANGE(2) &...
    [evINFO(:).smthwin] == SMTHWIN & [evINFO(:).basepoints] == BASEPOINTS &...
    [evINFO(:).preframes] == PREFRAMES & [evINFO(:).postframes] == POSTFRAMES);
listprocessed = cell(numel(idxprocessed),1);
for ii = 1:numel(idxprocessed), listprocessed{1,ii} = evINFO(idxprocessed(ii)).id; end
% processedtr = zeros(numel(listprocessed),1);
% for ii = 1: numel(listprocessed), processedtr(ii) = strcmp(ROILIST, listprocessed{ii}); end
listcurrent = [ROILIST(seltrace)];

%% Event detection
detect_events(detecttr, seltrace);

%% Initialize display parameters
hspace = 40;
vspace = 40;
spcol = 6;
sprow = 16;

whtxtbox = [(WHTRCFIG(1)-4)/spcol ((WHTRCFIG(2)-2*vspace)/sprow)*8];
whtracewin = [((WHTRCFIG(1)-4)/spcol)*4 ((WHTRCFIG(2)-3*vspace)/sprow)*9];
whhistwin = [((WHTRCFIG(1)-4)/spcol)*2 ((WHTRCFIG(2)-3*vspace)/sprow)*5];
whbutsm = [40 20];
whbutbg = [60 20];

whtxtbox = whtxtbox./WHTRCFIG; whtracewin = whtracewin./WHTRCFIG; whhistwin = whhistwin./WHTRCFIG;
whbutsm = whbutsm./WHTRCFIG; whbutbg = whbutbg./WHTRCFIG;
hspace = hspace/WHTRCFIG(1); vspace = vspace/WHTRCFIG(2);

txtbox1pos = [hspace 1-vspace-whtxtbox(2) whtxtbox];
txtbox2pos = [hspace txtbox1pos(2)-whtxtbox(2) whtxtbox];
histwinpos = [txtbox1pos(1)+txtbox1pos(3)+hspace vspace whhistwin];
tracewinpos = [txtbox1pos(1)+txtbox1pos(3)+hspace 1-vspace-whtracewin(2) whtracewin];
fwbutpos = [tracewinpos(1)+tracewinpos(3)/2-hspace/3-whbutsm(1) tracewinpos(2)-vspace*1.5 whbutsm];
backbutpos = [tracewinpos(1)+tracewinpos(3)/2+hspace/3 tracewinpos(2)-vspace*1.5 whbutsm];

evdetfig = figure('Position',POSITIONTRCFIG, 'Name', 'Event Detection Overview');

txtbox1 = uicontrol('parent', evdetfig, 'style', 'text', 'string', ['Processed with current settings', '***', listprocessed],...
    'units', 'normalized','position', txtbox1pos, 'fontsize', FONTSIZE+1,'backgroundcolor',[0.8 0.8 0.8]);
txtbox2 = uicontrol('parent', evdetfig, 'style', 'text', 'string', ['Currently processing', '***', listcurrent],...
    'units', 'normalized', 'position', txtbox2pos, 'fontsize', FONTSIZE+1,'backgroundcolor',[0.8 0.8 0.8]);
% Trace
tracewin = subplot('Position', tracewinpos, 'parent', evdetfig);
t1 = plot(timepts, detecttr(:,currdisp),'linewidth', LW1, 'color', TRCPLOTCOL1); t1.Color(4) = TRCALPHA;
title(detectid{currdisp}); set(gca, 'FONTSIZE', FONTSIZE);
xlabel(TRCXLABEL); ylabel(TRCYLABEL);
xlim([timepts(1), timepts(end)]);
if ~isnan(evINFO(seltrace(currdisp)).crossidx)
    if PLOTPEAKS, evidx = timeptstot(evINFO(seltrace(currdisp)).peakidx);
    else, evidx = timeptstot(evINFO(seltrace(currdisp)).crossidx);
    end
    for iE = 1:numel(evidx), hold on; xline(evidx(iE), 'Linewidth',LW1, 'Color', [0 0 0]); end
end

% Histogram
deltas = diff(detecttr,1,1); deltas = reshape(deltas,[numel(deltas) 1]);
deltas = sort(deltas, 'descend'); deltas = deltas(1:round(0.25*numel(deltas)));
histwin = subplot('Position', histwinpos, 'parent', evdetfig);
h = histogram(deltas, 20);
% Trace buttons
fwbut = uicontrol('parent', evdetfig, 'style', 'pushbutton',  'units', 'normalized',...
    'position', fwbutpos,'string', '<-','backgroundcolor',[0.8 0.8 0.8],...
    'FONTSIZE', FONTSIZE+2, 'callback', {@cb_fwbacktrace});
backbut = uicontrol('parent', evdetfig, 'style', 'pushbutton',...
    'units', 'normalized','position', backbutpos,...
    'string', '->','FONTSIZE', FONTSIZE+2, 'backgroundcolor',[0.8 0.8 0.8],...
    'callback', {@cb_fwbacktrace});


%% Local Callbacks
    function cb_fwbacktrace(hObj,~)
        if strcmp(hObj.String, '->') && currdisp < numtr
            currdisp = currdisp+1;
        elseif strcmp(hObj.String, '<-') && currdisp > 1
            currdisp = currdisp-1;
        end
        t1.YData = detecttr(:,currdisp);
        axes(tracewin);
        childs = get(gca,'children');
        delidx = [];
        for ii = 1:numel(childs)
            if strcmp(childs(ii).Type, 'constantline'), delidx = [delidx ii]; end
        end
        delete(childs(delidx));
        if ~isnan(evINFO(seltrace(currdisp)).crossidx)
            if PLOTPEAKS, evidx = timeptstot(evINFO(seltrace(currdisp)).peakidx);
            else, evidx = timeptstot(evINFO(seltrace(currdisp)).crossidx);
            end
            for iE = 1:numel(evidx), hold on; xline(evidx(iE), 'Linewidth',LW1, 'Color', [0 0 0]); end
        end
        title(detectid{currdisp});
    end



end