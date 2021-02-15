function [] = update_trc(sp2, trace, xrange, mode, showall, evinfo)

% Declare globally shared variables
global FONTSIZE TRCYLABEL TRCXLABEL TRCYLABELDFOF TRCPLOTCOL1 TRCPLOTCOL2 TRCALPHA...
    THRESHLW THRESHALPHA EVPLOTCOL traceINFO IMHW FTIME PLOTPEAKS 

subplot(sp2);
if strcmp(mode, '  Average') && ~showall
    t1 = plot(traceINFO(trace).timestamp{1}, traceINFO(trace).binned_roi_av{1},'linewidth', 1, 'color', TRCPLOTCOL1); t1.Color(4) = TRCALPHA;
    hold on, t2 = plot(traceINFO(trace).timestamp{1}, traceINFO(trace).smoothed{1},'linewidth', 2, 'color', 'black');
%     ty = yline(traceINFO(trace).events.threshold, 'color', [0 0 0],'linewidth', THRESHLW); ty.Color(4) = THRESHALPHA;
%     ty = yline(traceINFO(trace).events.peakthreshold, 'color', [0 0 0], 'linewidth', THRESHLW/2); ty.Color(4) = THRESHALPHA;
%     events = traceINFO(trace).events;
%     if PLOTPEAKS, evidx = traceINFO(trace).timestamp{1}(events.peaks);
%     else, evidx = traceINFO(trace).timestamp{1}(events.crossings);
%     end
%     for iE = 1:numel(evidx), hold on; xline(evidx(iE),'Linewidth',2, 'Color', EVPLOTCOL); end
    hold off;
    ylabel(TRCYLABEL); set(gca, 'FONTSIZE', FONTSIZE);
    xlim(xrange);
    
elseif strcmp(mode, '  Average') && showall
    ylimit = [min(traceINFO(trace).binned_roi_av{1},[],'all') max(traceINFO(trace).binned_roi_av{1},[],'all')]; ylimit = [ylimit(1)-0.1*ylimit(1) ylimit(2)+0.1*ylimit(2)];
    bgx = [traceINFO(trace).timestamp{1}(1) traceINFO(trace).timestamp{1}(1),...
        traceINFO(trace).timestamp{1}(end) traceINFO(trace).timestamp{1}(end)];
    bgy = ylimit; bgy = [bgy(1) bgy(2) bgy(2) bgy(1)];
    patch(bgx, bgy, [0.3 0.3 0.3], 'FaceAlpha', 0.2, 'LineStyle','none');
    hold on;
    t1 = plot(traceINFO(trace).tot_timestamp{1}, traceINFO(trace).tot_binned_roi_av{1},'linewidth', 1, 'color', TRCPLOTCOL1); t1.Color(4) = TRCALPHA;
    t2 = plot(traceINFO(trace).tot_timestamp{1}, traceINFO(trace).tot_smoothed{1},'linewidth', 2, 'color', 'black');
%     ty = yline(traceINFO(trace).tot_events.threshold, 'color', [0 0 0],'linewidth', THRESHLW); ty.Color(4) = THRESHALPHA;
%     ty = yline(traceINFO(trace).tot_events.peakthreshold, 'color', [0 0 0], 'linewidth', THRESHLW/2); ty.Color(4) = THRESHALPHA;
%     events = traceINFO(trace).tot_events;
%     if PLOTPEAKS, evidx = traceINFO(trace).tot_timestamp{1}(events.peaks);
%     else, evidx = traceINFO(trace).tot_timestamp{1}(events.crossings);
%     end
%     for iE = 1:numel(evidx), hold on; xline(evidx(iE),'Linewidth',2, 'Color', EVPLOTCOL); end
    hold off;
    ylabel(TRCYLABEL); set(gca, 'FONTSIZE', FONTSIZE);
    xlim([1*FTIME IMHW(1)*FTIME]); ylim(ylimit);
    
elseif strcmp(mode, '  DeltaF/F') && ~showall
    t1 = plot(traceINFO(trace).timestamp{1}, traceINFO(trace).dFoF_roi_av{1},'linewidth', 1, 'color', TRCPLOTCOL2);
    ylabel(TRCYLABELDFOF); set(gca, 'FONTSIZE', FONTSIZE);
    xlim(xrange);
    events = traceINFO(trace).events;
elseif strcmp(mode, '  DeltaF/F') && showall
    ylimit = [min(traceINFO(trace).tot_dFoF_roi_av{1},[],'all') max(traceINFO(trace).tot_dFoF_roi_av{1},[],'all')]; ylimit = [ylimit(1)-0.1*ylimit(1) ylimit(2)+0.1*ylimit(2)];
    bgx = [traceINFO(trace).timestamp{1}(1) traceINFO(trace).timestamp{1}(1),...
        traceINFO(trace).timestamp{1}(end) traceINFO(trace).timestamp{1}(end)];
    bgy = ylimit; bgy = [bgy(1) bgy(2) bgy(2) bgy(1)];
    patch(bgx, bgy, [0.3 0.3 0.3], 'FaceAlpha', 0.2, 'LineStyle','none');
    hold on;
    plot(traceINFO(trace).tot_timestamp{1}, traceINFO(trace).tot_dFoF_roi_av{1},'linewidth', 1, 'color', TRCPLOTCOL2);
    ylabel(TRCYLABELDFOF); set(gca, 'FONTSIZE', FONTSIZE);
    xlim([1*FTIME IMHW(1)*FTIME]); ylim(ylimit);
%     events = traceINFO(trace).tot_events;
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
end