function ev_detection_overview(detectall, singlesave)

% Declare globally shared variables
global WHTRCFIG POSITIONTRCFIG FONTSIZE TRCYLABEL TRCXLABEL TRCPLOTCOL1...
   PLOTPEAKS SMTHWIN FTIME TRACEDATA ROILIST RANGE TRACEID evINFO LW1 LW2...
   FTIMEVEC XTYPE THRESHOLD PEAKTHRESHOLD THRESHOLDTYPE XLINEALPHA

%% Traces to be displayed
if detectall, seltrace = 1:size(TRACEDATA,2);
else, seltrace = TRACEID;
end
selectedtr = TRACEDATA(RANGE(1):RANGE(2), seltrace);
detectid = ROILIST(seltrace);
numtr = size(selectedtr,2);
if strcmp(XTYPE,'s'),timepts = FTIMEVEC(RANGE(1):RANGE(2));
else, timepts = RANGE(1):RANGE(2);
end
if SMTHWIN ~= 0
    detecttr = zeros(size(selectedtr));
    for iTr = 1:size(selectedtr,2), detecttr(:,iTr) = smooth_data(selectedtr(:,iTr),SMTHWIN); end
else
    detecttr = selectedtr;
end

%% Initialize variables
currdisp = 1;
keepev = cell(numtr,1);
revised = false(numtr,1);

%% Extract settings of current trace if already analyzed
if strcmp(evINFO(seltrace(currdisp)).accepted, 'accepted')
    info = {sprintf('Range: %i - %i', evINFO(seltrace(currdisp)).rangein, evINFO(seltrace(currdisp)).rangeout),...
    sprintf('Smoothing window: %i', evINFO(seltrace(currdisp)).smthwin),...
    sprintf('Baseline frames: %i', evINFO(seltrace(currdisp)).baseframes),...
    sprintf('Threshold: %i | type: %s', round(evINFO(seltrace(currdisp)).threshold), evINFO(seltrace(currdisp)).thresholdtype)};
else
    info = {'not yet processed'};
end

% processedtr = zeros(numel(listprocessed),1);
% for ii = 1: numel(listprocessed), processedtr(ii) = strcmp(ROILIST, listprocessed{ii}); end
listcurrent = [ROILIST(seltrace)];

%% Event detection
detectidx = find(strcmp({evINFO(seltrace).accepted}, 'pending'));
detect_events(detecttr(:,detectidx), seltrace(detectidx));
for iTr = 1:numtr
    if any(iTr == detectidx)
        keepev{iTr,1} = false(numel(evINFO(iTr).crossidx),1);
    else
        keepev{iTr,1} = true(numel(evINFO(seltrace(iTr)).crossidx),1);
        revised(iTr,1) = true;
    end
end

%% Initialize display parameters
fs = FONTSIZE+1;
hspace = 40;
vspace = 40;
spcol = 6;
sprow = 16;

whtxtbox = [(WHTRCFIG(1))/spcol ((WHTRCFIG(2)-2*vspace)/sprow)*8];
whtracewin = [((WHTRCFIG(1))/spcol)*4 ((WHTRCFIG(2)-3*vspace)/sprow)*9];
whhistwin = [((WHTRCFIG(1))/spcol)*2 ((WHTRCFIG(2)-3*vspace)/sprow)*5];
whbutsm = [40 20];
whbutbg = [85 30];
whin = [50 20];
whtxt = [200 20];
whrb = [50 20];

whtxtbox = whtxtbox./WHTRCFIG; whtracewin = whtracewin./WHTRCFIG; whhistwin = whhistwin./WHTRCFIG;
whbutsm = whbutsm./WHTRCFIG; whbutbg = whbutbg./WHTRCFIG;
hspace = hspace/WHTRCFIG(1); vspace = vspace/WHTRCFIG(2);
whin = whin./WHTRCFIG; whtxt = whtxt./WHTRCFIG; whrb = whrb./WHTRCFIG;

txtbox1pos = [hspace 1-vspace-whtxtbox(2) whtxtbox];
txtbox2pos = [hspace txtbox1pos(2)-whtxtbox(2) whtxtbox];
histwinpos = [txtbox1pos(1)+txtbox1pos(3)+hspace vspace whhistwin];
tracewinpos = [txtbox1pos(1)+txtbox1pos(3)+hspace 1-vspace-whtracewin(2) whtracewin];
fwbutpos = [tracewinpos(1)+tracewinpos(3)/2-hspace/3-whbutsm(1) tracewinpos(2)-vspace*1.5 whbutsm];
backbutpos = [tracewinpos(1)+tracewinpos(3)/2+hspace/3 tracewinpos(2)-vspace*1.5 whbutsm];
threshtypetxtpos = [histwinpos(1)+histwinpos(3)+hspace*0.7 histwinpos(2)+histwinpos(4)-whtxt(2)-vspace/2 whtxt(1)*0.5 whtxt(2)];
threshtyperbpos = [threshtypetxtpos(1) threshtypetxtpos(2)-whtxt(2) whrb;...
    threshtypetxtpos(1)+whrb(1)*1.5 threshtypetxtpos(2)-whtxt(2) whrb];
threshtxtpos = [threshtypetxtpos(1) threshtyperbpos(2,2)-whrb(2)-vspace/3 whtxt(1)*0.8 whtxt(2)];
threshinpos = [threshtxtpos(1) threshtxtpos(2)-whin(2)-vspace/8 whin];
peakinpos = [threshinpos(1)+threshinpos(3)+hspace/5 threshtxtpos(2)-whin(2)-vspace/8 whin];

statetxtpos = [threshtxtpos(1)+threshtxtpos(3)+0.6*hspace histwinpos(2)+histwinpos(4) whtxt(1)*0.7 whtxt(2)*2];
acceptbutpos = [statetxtpos(1)+hspace*0.8 statetxtpos(2)-whbutbg(2)-vspace/3 whbutbg];
revisebutpos = [acceptbutpos(1) acceptbutpos(2)-whbutbg(2)-vspace/3 whbutbg];
closebutpos = [acceptbutpos(1) revisebutpos(2)-whbutbg(2)-vspace/3 whbutbg];

evdetfig = figure('Position',POSITIONTRCFIG, 'Name', 'Event Detection Overview');
% Overview
txtbox1 = uicontrol('parent', evdetfig, 'style', 'text', 'string', ['Processed with following settings:', '***', info],...
    'units', 'normalized','position', txtbox1pos, 'fontsize', fs,'backgroundcolor',[0.8 0.8 0.8]);
txtbox2 = uicontrol('parent', evdetfig, 'style', 'text', 'string', ['Currently processing', '***', listcurrent],...
    'units', 'normalized', 'position', txtbox2pos, 'fontsize', fs,'backgroundcolor',[0.8 0.8 0.8]);

% Buttons and Input fields
threshtxt = uicontrol('parent', evdetfig, 'style', 'text','units', 'normalized','position', threshtxtpos,'string', 'Threshold [deltaF]:','FONTSIZE', FONTSIZE);
threshin = uicontrol('parent', evdetfig, 'style', 'edit','units', 'normalized','position', threshinpos,'FONTSIZE', FONTSIZE, 'String', num2str(THRESHOLD),'Callback', {@cb_threshin});
peakin = uicontrol('parent', evdetfig, 'style', 'edit','units', 'normalized','position', peakinpos, 'Enable', 'off', 'FONTSIZE', FONTSIZE, 'String', num2str(PEAKTHRESHOLD),'Callback', {@cb_peakthreshin});

threshtypetxt = uicontrol('parent', evdetfig, 'style', 'text','units', 'normalized','position', threshtypetxtpos,'string', 'Threshold type:','FONTSIZE', FONTSIZE);
threshtyperb1 = uicontrol('parent', evdetfig, 'style', 'radiobutton', 'units', 'normalized','position', threshtyperbpos(1,:),'string', 'delta F', 'Value', 1,'FONTSIZE', FONTSIZE, 'callback', {@cb_threshtyperb});
threshtyperb2 = uicontrol('parent', evdetfig, 'style', 'radiobutton', 'units', 'normalized','position', threshtyperbpos(2,:),'string', 'sd', 'Value', 0,'FONTSIZE', FONTSIZE, 'callback', {@cb_threshtyperb});

statetxt = uicontrol('parent', evdetfig, 'style', 'text','units', 'normalized','position', statetxtpos,...
    'string', {'Event detection state:', evINFO(seltrace(currdisp)).accepted},'FONTSIZE', fs+1,...
    'BackgroundColor', [0.2 0.2 0.2], 'ForegroundColor', [1 1 1]);
acceptbut = uicontrol('parent', evdetfig, 'style', 'pushbutton',  'units', 'normalized','position', acceptbutpos,'string', 'Accept events',...
    'FONTSIZE', fs, 'callback', {@cb_acceptbut});
revisebut = uicontrol('parent', evdetfig, 'style', 'pushbutton',  'units', 'normalized','position', revisebutpos,'string', 'Revision',...
    'FONTSIZE', fs, 'callback', {@cb_revisebut});
closebut = uicontrol('parent', evdetfig, 'style', 'pushbutton',  'units', 'normalized','position', closebutpos,'string', 'Close',...
    'FONTSIZE', fs, 'callback', {@cb_closebut});
if isempty(evINFO(seltrace(currdisp)).crossidx)
    acceptbut.Enable = 'off'; revisebut.Enable = 'off';
elseif strcmp(evINFO(seltrace(currdisp)).accepted,'accepted')
    acceptbut.String = 'Re-run detection';
else
    acceptbut.Enable = 'on'; revisebut.Enable = 'on';
end
        
% Trace
tracewin = subplot('Position', tracewinpos, 'parent', evdetfig);
plot_trace();

% Histogram
deltas = diff(detecttr(:,currdisp));
deltas = sort(deltas, 'descend'); deltas = deltas(1:round(0.5*numel(deltas)));
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
        plot_trace();
        statetxt.String = {'Event detection state:', evINFO(seltrace(currdisp)).accepted};
        if isempty(evINFO(seltrace(currdisp)).crossidx)
            acceptbut.Enable = 'off'; revisebut.Enable = 'off';
            acceptbut.String = 'Accept events';
            info = {'no events detected'};
            txtbox1.String = ['Processed with following settings:', '***', info];
        elseif strcmp(evINFO(seltrace(currdisp)).accepted,'accepted')
            info = {sprintf('Range: %i - %i', evINFO(seltrace(currdisp)).rangein, evINFO(seltrace(currdisp)).rangeout),...
                sprintf('Smoothing window: %i', evINFO(seltrace(currdisp)).smthwin),...
                sprintf('Baseline frames: %i', evINFO(seltrace(currdisp)).baseframes),...
                sprintf('Threshold: %i | type: %s', round(evINFO(seltrace(currdisp)).threshold), evINFO(seltrace(currdisp)).thresholdtype)};
            txtbox1.String = ['Processed with following settings:', '***', info];
            acceptbut.String = 'Re-run detection';
        else
            acceptbut.Enable = 'on'; revisebut.Enable = 'on';
            info = {'not yet accepted'};
            txtbox1.String = ['Processed with following settings:', '***', info];
            acceptbut.String = 'Accept events';
        end
        % Histogram update
        deltas = diff(detecttr(:,currdisp));
        deltas = sort(deltas, 'descend'); deltas = deltas(1:round(0.5*numel(deltas)));
        axes(histwin); h = histogram(deltas, 20);
    end

    function cb_threshtyperb(hObj,~)
        if strcmp(hObj.String, 'delta F')
            THRESHOLDTYPE = 'dF';
            threshtyperb2.Value = 0;
            peakin.Enable = 'off';
            threshtxt.String = 'Threshold [deltaF]:';
        else
            THRESHOLDTYPE = 'sd';
            threshtyperb1.Value = 0;
            peakin.Enable = 'on';
            threshtxt.String = 'Threshold [s.d.]: start | peak';
        end
    end

    function cb_threshin(hObj,~)
        THRESHOLD = str2double(get(hObj,'String'));
        % Apply new threshold to detection
        detectidx = find(strcmp({evINFO(seltrace).accepted}, 'pending'));
        detect_events(detecttr(:,detectidx), seltrace(detectidx));
        for iTr = 1:numtr
            if any(iTr == detectidx)
                keepev{iTr,1} = false(numel(evINFO(iTr).crossidx),1);
            else
                keepev{iTr,1} = true(numel(evINFO(seltrace(iTr)).crossidx),1);
                revised(iTr,1) = true;
            end
        end
        plot_trace();
    end

    function cb_peakthreshin(hObj,~)
        PEAKTHRESHOLD = str2double(get(hObj,'String'));
        % Apply new threshold to detection
        detectidx = find(strcmp({evINFO(seltrace).accepted}, 'pending'));
        detect_events(detecttr(:,detectidx), seltrace(detectidx));
        for iTr = 1:numtr
            if any(iTr == detectidx)
                keepev{iTr,1} = false(numel(evINFO(iTr).crossidx),1);
            else
                keepev{iTr,1} = true(numel(evINFO(seltrace(iTr)).crossidx),1);
                revised(iTr,1) = true;
            end
        end
        plot_trace();
    end


    function cb_acceptbut(hObj,~)
        if strncmp(hObj.String, 'Re',2)
            detect_events(detecttr(:,currdisp), seltrace(currdisp));
            keepev{currdisp,1} = false(numel(evINFO(seltrace(currdisp)).crossidx),1);
            revised(currdisp,1) = false;
            plot_trace();
            acceptbut.String = 'Accept events';
            revisebut.Enable = 'on';
            evINFO(seltrace(currdisp)).accepted = 'pending';
        else
            if ~revised(currdisp,1)
                keepev{currdisp,1} = true(numel(evINFO(seltrace(currdisp)).crossidx),1);
            end
            update_evinfo(seltrace(currdisp), keepev{currdisp,1});
            keepev{currdisp,1} = true(sum(keepev{currdisp,1}),1);
            if singlesave, save_files(seltrace(currdisp)); end
            statetxt.String = {'Event detection state:', evINFO(seltrace(currdisp)).accepted};
            revisebut.Enable = 'off';
            acceptbut.String = 'Re-run detection';
        end
        info = {sprintf('Range: %i - %i', evINFO(seltrace(currdisp)).rangein, evINFO(seltrace(currdisp)).rangeout),...
                sprintf('Smoothing window: %i', evINFO(seltrace(currdisp)).smthwin),...
                sprintf('Baseline frames: %i', evINFO(seltrace(currdisp)).baseframes),...
                sprintf('Threshold: %i | type: %s', round(evINFO(seltrace(currdisp)).threshold), evINFO(seltrace(currdisp)).thresholdtype)};
        txtbox1.String = ['Processed with following settings:', '***', info];
    end

    function cb_revisebut(~,~)
        goOn = false;
        while ~goOn
            [goOn, currkeepev] = revise_events(detecttr(:,currdisp), seltrace(currdisp));
        end
        keepev{currdisp,1} = currkeepev;
        revised(currdisp,1) = true;
        plot_trace();
    end

    function cb_closebut(~,~)
        close gcf;
    end

    function plot_trace()
        axes(tracewin);
        cla;
        t1 = plot(timepts, detecttr(:,currdisp),'linewidth', LW1, 'color', TRCPLOTCOL1);
        title(detectid{currdisp}); set(gca, 'FONTSIZE', FONTSIZE);
        xlabel(TRCXLABEL); ylabel(TRCYLABEL);
        xlim([timepts(1), timepts(end)]);
        if ~revised(currdisp,1) && any(~isnan(evINFO(seltrace(currdisp)).crossidx))
            if PLOTPEAKS, evidx = timepts(evINFO(seltrace(currdisp)).peakidx);
            else, evidx = timepts(evINFO(seltrace(currdisp)).crossidx);
            end
            for iE = 1:numel(evidx), hold on; xl = xline(evidx(iE), 'Linewidth',LW2, 'Color', [0 0 0]); xl.Color(4) = XLINEALPHA; end
        elseif revised(currdisp,1) && sum(keepev{currdisp,1}) > 0
            if PLOTPEAKS, evidx = timepts(evINFO(seltrace(currdisp)).peakidx(keepev{currdisp,1}));
            else, evidx = timepts(evINFO(seltrace(currdisp)).crossidx(keepev{currdisp,1}));
            end
            for iE = 1:numel(evidx), hold on; xl = xline(evidx(iE), 'Linewidth',LW2, 'Color', [0 0 0]); xl.Color(4) = XLINEALPHA; end
        end
    end
end


%         t1.YData = detecttr(:,currdisp);
%         axes(tracewin);
%         childs = get(gca,'children');
%         delidx = [];
%         for ii = 1:numel(childs)
%             if strcmp(childs(ii).Type, 'constantline'), delidx = [delidx ii]; end
%         end
%         delete(childs(delidx));
%         if ~isnan(evINFO(seltrace(currdisp)).crossidx)
%             if PLOTPEAKS, evidx = FTIMEVEC(evINFO(seltrace(currdisp)).peakidx);
%             else, evidx = FTIMEVEC(evINFO(seltrace(currdisp)).crossidx);
%             end
%             for iE = 1:numel(evidx), hold on; xline(evidx(iE), 'Linewidth',LW2, 'Color', [0 0 0]); end
%         end
%         title(detectid{currdisp});
