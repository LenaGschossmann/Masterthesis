function [] = update_trc(sp2, trace, mode, showall)

% Declare globally shared variables
global FONTSIZE TRCYLABEL TRCXLABEL TRCYLABELDFOF TRCPLOTCOL1 TRCPLOTCOL2 TRCALPHA...
    traceINFO IMHW FTIME EVDATA PLOTRANGE

subplot(sp2);cla;
tmproi = traceINFO(trace).roiID;

if strcmp(mode, '  Average') && ~showall
    timestamps = PLOTRANGE(1) *FTIME : FTIME : PLOTRANGE(2)*FTIME;
    t1 = plot(timestamps, traceINFO(trace).roi_av{1}(PLOTRANGE(1):PLOTRANGE(2)),'linewidth', 1, 'color', TRCPLOTCOL1); t1.Color(4) = TRCALPHA;
    plot_events(tmproi, PLOTRANGE, false);
    hold off;
    ylabel(TRCYLABEL); set(gca, 'FONTSIZE', FONTSIZE);
    xlim(PLOTRANGE.*FTIME);
    
elseif strcmp(mode, '  Average') && showall
    timestamps = 1*FTIME : FTIME : IMHW(1)*FTIME;
    ylimit = [min(traceINFO(trace).roi_av{1},[],'all') max(traceINFO(trace).roi_av{1},[],'all')]; ylimit = [ylimit(1)-0.1*ylimit(1) ylimit(2)+0.1*ylimit(2)];
    bgx = [PLOTRANGE(1)*FTIME PLOTRANGE(1)*FTIME, PLOTRANGE(2)*FTIME PLOTRANGE(2)*FTIME];
    bgy = ylimit; bgy = [bgy(1) bgy(2) bgy(2) bgy(1)];
    patch(bgx, bgy, [0.3 0.3 0.3], 'FaceAlpha', 0.2, 'LineStyle','none');
    hold on;
    t1 = plot(timestamps, traceINFO(trace).roi_av{1},'linewidth', 1, 'color', TRCPLOTCOL1); t1.Color(4) = TRCALPHA;
    plot_events(tmproi, [1 IMHW(1)],false);
    hold off;
    ylabel(TRCYLABEL); set(gca, 'FONTSIZE', FONTSIZE);
    xlim([1*FTIME IMHW(1)*FTIME]); ylim(ylimit);
    
elseif strcmp(mode, '  DeltaF/F') && ~showall
    timestamps = PLOTRANGE(1) *FTIME : FTIME : PLOTRANGE(2)*FTIME;
    t1 = plot(timestamps, traceINFO(trace).dFoF_roi_av{1}(PLOTRANGE(1):PLOTRANGE(2)),'linewidth', 1, 'color', TRCPLOTCOL2);
    ylabel(TRCYLABELDFOF); set(gca, 'FONTSIZE', FONTSIZE);
    xlim(PLOTRANGE.*FTIME);
    plot_events(tmproi, PLOTRANGE,true);
    
elseif strcmp(mode, '  DeltaF/F') && showall
    timestamps = 1*FTIME : FTIME : IMHW(1)*FTIME;
    ylimit = [min(traceINFO(trace).dFoF_roi_av{1},[],'all') max(traceINFO(trace).dFoF_roi_av{1},[],'all')]; ylimit = [ylimit(1)-0.1*ylimit(1) ylimit(2)+0.1*ylimit(2)];
    bgx = [PLOTRANGE(1)*FTIME PLOTRANGE(1)*FTIME, PLOTRANGE(2)*FTIME PLOTRANGE(2)*FTIME];
    bgy = ylimit; bgy = [bgy(1) bgy(2) bgy(2) bgy(1)];
    patch(bgx, bgy, [0.3 0.3 0.3], 'FaceAlpha', 0.2, 'LineStyle','none');
    hold on;
    plot(timestamps, traceINFO(trace).dFoF_roi_av{1},'linewidth', 1, 'color', TRCPLOTCOL2);
    ylabel(TRCYLABELDFOF); set(gca, 'FONTSIZE', FONTSIZE);
    xlim([1*FTIME IMHW(1)*FTIME]); ylim(ylimit);
    plot_events(tmproi, [1 IMHW(1)],true);
    hold off;
end

xlabel(TRCXLABEL);
% annotxt = {sprintf('Average | SD: %s  |  %s', num2str(events.average,3), num2str(events.sd,3)),...
%     sprintf('Threshold: %s', num2str(events.threshold,3)),...
%     sprintf('Threshold peak: %s', num2str(events.peakthreshold,3)),...
%     sprintf('Event type: %s', events.eventtype),...
%     sprintf('Av. Inter-event-interval: %s s', num2str(events.aviei,3)), ...
%     sprintf('CV IEI: %s', num2str(events.cviei,3)),...
%     sprintf('Av. Amplitude: %s', num2str(events.avamp,3)),...
%     sprintf('Av. Eventrate: %s Hz', num2str(events.eventrate,3))};
% set(evinfo, 'string', annotxt);

%% Local Callback
    function plot_events(tmproi, range, plotlines)
        if plotlines
            if size(EVDATA,1) >= tmproi && ~isempty(EVDATA{tmproi,10})
                hold on, yline(EVDATA{tmproi,10}(2)*EVDATA{tmproi,10}(3), 'cyan');
                text(1,EVDATA{tmproi,10}(2)*EVDATA{tmproi,10}(3), strcat('det.thresh. x ',num2str(EVDATA{tmproi,10}(3))), 'Fontsize', 6);
            end
            if size(EVDATA,1) >= tmproi && ~isempty(EVDATA{tmproi,10})
                yline(EVDATA{tmproi,10}(1)*EVDATA{tmproi,10}(4), 'red');
                text(1,EVDATA{tmproi,10}(1)*EVDATA{tmproi,10}(4), strcat('%-tile val. x ',num2str(EVDATA{tmproi,10}(4))), 'Fontsize', 6);
            end
            if size(EVDATA,1) >= tmproi && ~isempty(EVDATA{tmproi,10})
                hold on, yline(-EVDATA{tmproi,10}(1), 'Color',[.5 .5 .5]);
                yline(EVDATA{tmproi,10}(1), 'Color',[.5 .5 .5]);
            end
        end
        if size(EVDATA,1) >= tmproi && ~(isempty(EVDATA{tmproi,6}) && isempty(EVDATA{tmproi,7}))
            save_idx = EVDATA{tmproi,7}; save_idx = save_idx(save_idx > range(1) & save_idx < range(2));
            save_pts = save_idx*FTIME;
            rev_idx = EVDATA{tmproi,6}; rev_idx = rev_idx(rev_idx > range(1) & rev_idx < range(2));
            rev_pts = rev_idx.*FTIME;
            hold on;
            for iEv = 1:numel(save_pts), xline(save_pts(iEv), 'Color',[0 0 0], 'LineWidth', 1.2, 'Alpha',.4); end
            for iEv = 1:numel(rev_pts), xline(rev_pts(iEv), 'r', 'LineWidth', 1.2, 'Alpha', .4); end
            hold off;
        end
    end
end