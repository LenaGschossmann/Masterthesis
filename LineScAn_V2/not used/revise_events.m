function [goOn, keepev] = revise_events(worktrace, traceidx)
%% Declare globally shared variables
global SCRSZ FTIME TRCPLOTCOL1 TRCALPHA evINFO LW1 LW2 PREPOINTS POSTPOINTS FONTSIZE...
    XLINEALPHA FTIMEVEC

%% Initialize display parameters
fs = FONTSIZE+1;
whrevfig = [400 400];
hspace = 40;
vspace = 30;

whbut = [85 30];
whrb = [70 20];
whtxt = [200 20];
whtracewin = [whrevfig(1)-round(hspace*2.5) whrevfig(2)-vspace*3-whbut(2)*2];

whtracewin = whtracewin./whrevfig;
whbut = whbut./whrevfig;
hspace = hspace/whrevfig(1); vspace = vspace/whrevfig(2);
whrb = whrb./whrevfig; whtxt = whtxt./whrevfig;

revfigpos = [round(SCRSZ(1)/2-whrevfig(1)) round(SCRSZ(2)/2-whrevfig(2)/2) whrevfig];
trcplotpos = [hspace*1.5 vspace*2+2*whbut(2) whtracewin];
acceptbutpos = [0.5-hspace*0.2-whbut(1) vspace whbut];
rejectbutpos = [0.5+hspace*0.2 acceptbutpos(2) whbut];
showrbpos = [rejectbutpos(1)+rejectbutpos(3)+hspace/2 vspace whrb];

%% Initialize und unpack varaibles
numev = numel(evINFO(traceidx).crossidx);
currev = 1;
showall = false;
keepev = false(numev,1);
crossidx = evINFO(traceidx).crossidx;
baseline = evINFO(traceidx).baselinevalues;
croptrace = NaN(numev,PREPOINTS+POSTPOINTS+1);
xpoints = (0:PREPOINTS+POSTPOINTS).*FTIME.*1000; % in [ms]

%% Extrace cropped out region around events
for iE = 1:numev
    % Event preceding part
    if crossidx(iE)-PREPOINTS < 1
        croptrace(iE,PREPOINTS-crossidx(iE)+2:PREPOINTS+1) = worktrace(1:crossidx(iE))';
    else
        croptrace(iE,1:PREPOINTS+1) = worktrace(crossidx(iE)-PREPOINTS:crossidx(iE))';
    end
    % Event following part
    if crossidx(iE)+POSTPOINTS > numel(worktrace)
        croptrace(iE,PREPOINTS+2:PREPOINTS+1+numel(worktrace)-crossidx(iE)) = worktrace(crossidx(iE)+1:end);
    else
        croptrace(iE,PREPOINTS+2:end) = worktrace(crossidx(iE)+1:crossidx(iE)+POSTPOINTS);
    end
    % Delta F over F
    croptrace(iE,:) = (croptrace(iE,:)-baseline(iE))./baseline(iE);
end

%% Display
revfig = figure('Position',revfigpos, 'Name', 'Event Revision');
trcplot = subplot('position', trcplotpos);
acceptbut = uicontrol('parent', revfig, 'style', 'pushbutton','units', 'normalized', 'position', acceptbutpos,'string', 'ACCEPT','BackgroundColor', '#77AC30','FONTSIZE', fs, 'callback', {@cb_acceptbut});
rejectbut = uicontrol('parent', revfig, 'style', 'pushbutton','units', 'normalized', 'position', rejectbutpos,'string', 'REJECT','BackgroundColor', '#A2142F','FONTSIZE', fs, 'callback', {@cb_rejectbut});
showrb = uicontrol('parent', revfig, 'style', 'radiobutton','units', 'normalized', 'position', showrbpos,'string', 'Show all', 'Value', 0,'FONTSIZE', FONTSIZE, 'callback', {@cb_showrb});

plot_traces();

%% Apply revision
uiwait();

%% Local Callbacks
    function cb_acceptbut(~,~)
        keepev(currev) = true;
        if currev < numev
            currev = currev+1; plot_traces();
        else
            pause(0.1);
            goOn = true;
            close gcf;
        end
    end

    function cb_rejectbut(~,~)
        keepev(currev) = false;
        if currev < numev
            currev = currev+1;plot_traces();
        else
            pause(0.1);
            goOn = true;
            close gcf;
        end
    end

    function cb_showrb(hObj,~)
        showall = logical(hObj.Value);
        plot_traces();
    end

    function plot_traces()
        axes(trcplot);
        cla;
        if showall
            for iE = 1:numev
                p=plot(xpoints,croptrace(iE,:),'linewidth', LW1, 'color', TRCPLOTCOL1); p.Color(4) = TRCALPHA;
                hold on;
            end
        else
            if currev > 1
                plotidx = find(keepev(1:currev-1) == true);
                for iE = 1:numel(plotidx)
                    p=plot(xpoints,croptrace(plotidx(iE),:),'linewidth', LW1, 'color', TRCPLOTCOL1); p.Color(4) = TRCALPHA;
                    hold on;
                end
            end
        end
        % Current trace
        p=plot(xpoints,croptrace(currev,:),'linewidth', LW2, 'color', [0 0 0]);
        xl=xline(xpoints(PREPOINTS+1), 'Linewidth',LW2+1, 'Color', [0 0 0]); alpha(XLINEALPHA); %xl.Color(4) = XLINEALPHA;
        xlim([xpoints(1) xpoints(end)]);
        xtckmarks = xpoints(1:5:end);
        xlabs = cell(numel(xtckmarks),1);
        for ii = 1:numel(xtckmarks), xlabs{ii} = num2str(xtckmarks(ii)); end
        xticks(xtckmarks); set(gca, 'XTickLabel', xlabs, 'FontSize', FONTSIZE-1); xtickangle(45);
        xlabel('time [ms]'); ylabel('dF/F');
        yrange = [min(croptrace, [], 'all') max(croptrace, [], 'all')];
        ylim([yrange(1)-diff(yrange)/10 yrange(2)+diff(yrange)/10]);
    end


end