function [] = update_trc(sp2, trace, xrange, mode, showall)

% Declare globally shared variables
global FONTSIZE TRCYLABEL TRCXLABEL TRCYLABELDFOF TRCPLOTCOL1 TRCPLOTCOL2 TRCALPHA...
    traceINFO IMHW FTIME EVDATA COMPOSITE2D

subplot(sp2);
tmproi = traceINFO(trace).roiID;
pr = traceINFO(trace).fig_params{1,1};
prtot = [1 size(COMPOSITE2D,1)];

if strcmp(mode, '  Average') && ~showall
    t1 = plot(traceINFO(trace).timestamp{1}, traceINFO(trace).binned_roi_av{1},'linewidth', 1, 'color', TRCPLOTCOL1); t1.Color(4) = TRCALPHA;
    plot_events(tmproi, traceINFO(trace).timestamp{1}, pr,false);
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
    plot_events(tmproi, traceINFO(trace).tot_timestamp{1}, prtot,false);
    hold off;
    ylabel(TRCYLABEL); set(gca, 'FONTSIZE', FONTSIZE);
    xlim([1*FTIME IMHW(1)*FTIME]); ylim(ylimit);
    
elseif strcmp(mode, '  DeltaF/F') && ~showall
    t1 = plot(traceINFO(trace).timestamp{1}, traceINFO(trace).dFoF_roi_av{1},'linewidth', 1, 'color', TRCPLOTCOL2);
    ylabel(TRCYLABELDFOF); set(gca, 'FONTSIZE', FONTSIZE);
    xlim(xrange);
    plot_events(tmproi, traceINFO(trace).timestamp{1}, pr,true);
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
    plot_events(tmproi, traceINFO(trace).tot_timestamp{1}, prtot,true);
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
    function plot_events(tmproi, tpts, pr, plotlines)
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
            save_idx = EVDATA{tmproi,7}; save_idx = save_idx(save_idx > pr(1) & save_idx < pr(2));
            save_pts = tpts(save_idx-pr(1)+1);
            rev_idx = EVDATA{tmproi,6}; rev_idx = rev_idx(rev_idx > pr(1) & rev_idx < pr(2));
            rev_pts = tpts(rev_idx-pr(1)+1);
            hold on;
            for iEv = 1:numel(save_pts), xline(save_pts(iEv), 'Color',[0 0 0], 'LineWidth', 1.2, 'Alpha',.4); end
            for iEv = 1:numel(rev_pts), xline(rev_pts(iEv), 'r', 'LineWidth', 1.2, 'Alpha', .4); end
            hold off;
        end
    end
end