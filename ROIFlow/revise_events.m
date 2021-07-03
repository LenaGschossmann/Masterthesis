function [goOn, keepev] = revise_events(auto, fig3_pos, params, eventdata, roiidx, ft)

%% Parameters
ev_idx = eventdata(roiidx).onset_idx;
num_ev = numel(ev_idx);
trace = eventdata(roiidx).filtrd_FoFtrace;
rnd_peak_FoF =  eventdata(roiidx).rnd_peak_FoF;
ctrl_peak_FoF = eventdata(roiidx).ctrl_peak_FoF;
peak_FoF = trace(eventdata(roiidx).peak_idx);
binedges = [min([peak_FoF rnd_peak_FoF ctrl_peak_FoF]) max([peak_FoF rnd_peak_FoF ctrl_peak_FoF])];
binedges = linspace(binedges(1),binedges(2),20);

prepts = round(params.crp_trc_pre_s/ft);
postpts = round(params.crp_trc_post_s/ft);
totpts = prepts+postpts+1;
ev_croptrc = NaN(num_ev,totpts);

currev = 1;
if auto
    keepev = true(num_ev,1);
    autokeep = keepev;
else
    keepev = false(num_ev,1); % Will be decided by user
    autokeep = eventdata(roiidx).auto_keep_ev;
end
xpoints = (1:totpts).*ft.*1000;
goOn = [];

%% GUI
whfig3 = fig3_pos(3:4);
whunits = 1./[12 12];
whbut3 = [2 1].*whunits;
% whtxt = [];
whtrcsp = [10 7].*whunits;
whhistsp = [2 1.5].*whunits;
hvspace3 = [1 1].*whunits;

trccol1 = [0.6510 0.0667 0.2784];
trccol2 = [0 0 0];
fontsize = 10;
an=[]; an2 = [];

but3_pos_1 = [hvspace3(1)*3 hvspace3(2) whbut3];
but3_pos_2 = [but3_pos_1(1)+but3_pos_1(3)+0.5*hvspace3(1) hvspace3(2) whbut3];
but3_pos_3 = [but3_pos_2(1)+but3_pos_2(3)+0.5*hvspace3(1) hvspace3(2) whbut3];
trcsp_pos = [hvspace3(1)*1.25 but3_pos_1(2)+but3_pos_1(4)+2*hvspace3(2) whtrcsp];
hist_pos = [trcsp_pos(1)+trcsp_pos(3)*0.62 trcsp_pos(2)+0.72*trcsp_pos(4) whhistsp];
leg_pos = [.82 .75 0.1 .1];
an_pos = [.82 .8 .1 .1];
an_pos2 = [.12 .8 .1 .1];

revfig = figure('Name', 'Event revision', 'Position', fig3_pos, 'toolbar', 'none', 'menu', 'none'); axis off;
set(revfig,'WindowKeyPressFcn',@keyPressCallback);
if auto, str1 = '<-'; str2 = '->'; else, str1 = 'Accept(A)'; str2 = 'Reject (D)'; end
but3_1 = uicontrol('parent', revfig, 'style', 'pushbutton','units', 'normalized','position', but3_pos_1,'string', str1,'fontsize', fontsize-1, 'callback', {@cb_but3_accept});
but3_2 = uicontrol('parent', revfig, 'style', 'pushbutton','units', 'normalized','position', but3_pos_2,'string', str2,'fontsize', fontsize-1, 'callback', {@cb_but3_reject});
if ~auto
    but3_1.BackgroundColor = '#77AC30';
    but3_2.BackgroundColor = '#A2142F';
    but3_3 = uicontrol('parent', revfig, 'style', 'pushbutton','units', 'normalized','position', but3_pos_3,'string', 'Back (B)','fontsize', fontsize-1, 'callback', {@cb_but3_back});
end

%% Extract cropped out region around events
% Revision events
for iE = 1:numel(ev_idx)
    % Event preceding part
    if ev_idx(iE)-prepts < 1, ev_croptrc(iE,prepts-ev_idx(iE)+2:prepts+1) = trace(1:ev_idx(iE))';
    else, ev_croptrc(iE,1:prepts+1) = trace(ev_idx(iE)-prepts:ev_idx(iE))';
    end
    % Event following part
    if ev_idx(iE)+postpts > numel(trace), ev_croptrc(iE,prepts+2:prepts+1+numel(trace)-ev_idx(iE)) = trace(ev_idx(iE)+1:end);
    else, ev_croptrc(iE,prepts+2:end) = trace(ev_idx(iE)+1:ev_idx(iE)+postpts);
    end
end

av_crp_trc = mean(ev_croptrc,1, 'omitnan');
min_crp_trc = min(ev_croptrc,[],1);
max_crp_trc = max(ev_croptrc,[],1);

if ~isempty(av_crp_trc), plot_traces(); end

%% Apply revision
uiwait();
if isempty(goOn)
    goOn = true;
    keepev(currev:end) = [];
end

%% Local functions
    function keyPressCallback(~,ev)
        if ~auto && strcmp(ev.Key,'a')
            keepev(currev) = true;
            if currev < num_ev, currev = currev+1; plot_traces();
            else, pause(0.1); goOn = true; close gcf;
            end
        elseif ~auto && strcmp(ev.Key,'d')
            keepev(currev) = false;
            if currev < num_ev, currev = currev+1; plot_traces();
            else, pause(0.1); goOn = true; close gcf;
            end
        elseif ~auto && strcmp(ev.Key,'b')
            if currev > 1
                currev = currev-1;
                keepev(currev) = false;
                plot_traces();
            end
        end
    end

    function cb_but3_accept(~,~)
        if ~auto
            keepev(currev) = true;
            if currev < num_ev, currev = currev+1; plot_traces();
            else, pause(0.1); goOn = true; close gcf;
            end
        else
            if currev > 1
                currev = currev-1;
                keepev(currev) = false;
                plot_traces();
            end
        end
    end

    function cb_but3_reject(~,~)
        if ~auto, keepev(currev) = false; end
        if currev < num_ev, currev = currev+1; plot_traces();
        else, pause(0.1); goOn = true; close gcf;
        end
    end

    function cb_but3_back(~,~)
        if currev > 1
            currev = currev-1;
            keepev(currev) = false;
            plot_traces();
        end
    end


    function plot_traces()
        figure(revfig);
        subplot('Position', trcsp_pos);
        cla;
        if isempty(an), an = annotation('textbox', an_pos, 'FitBoxToText','on','EdgeColor', [1 1 1],'fontsize', fontsize); end
        if isempty(an2), an2 = annotation('textbox', an_pos2, 'FitBoxToText','on','EdgeColor', [1 1 1],'fontsize', fontsize+2); end

        %         for iE = 1:num_ev
%             p = plot(xpoints,save_croptrc(iE,:),'linewidth', 1, 'color', trccol1); p.Color(4) = .25;
%             hold on;
%         end

        % Plot mean/min/max
        hold on
        pm1 = plot(xpoints,min_crp_trc,'linewidth', .5, 'color', [.5 .5 .5]); pm1.Color(4) = .1;
        pm2 = plot(xpoints,max_crp_trc,'linewidth', .5, 'color', [.5 .5 .5]); pm2.Color(4) = .1;
        pa = patch([xpoints fliplr(xpoints)], [min_crp_trc fliplr(max_crp_trc)],'g');
        pa.FaceColor = [.5 .5 .5]; pa.FaceAlpha = .1; pa.EdgeAlpha = .05; 
        pm3 = plot(xpoints,av_crp_trc,'linewidth', 2.5, 'color', [.5 .5 .5]); pm3.Color(4) = .8;
        
        if currev > 1
            plotidx = find(keepev(1:currev-1) == true);
            for iE = 1:numel(plotidx)
                hold on; p = plot(xpoints,ev_croptrc(plotidx(iE),:),'linewidth', 1, 'color', trccol1); p.Color(4) = .25;
            end
        end
        
        % Current trace
        if ~all(isnan(ev_croptrc(currev,:)))
            hold on;
            plot(xpoints,ev_croptrc(currev,:),'linewidth', 1.2, 'color', trccol2);
            an.String = sprintf('%i/%i', currev,num_ev);
        end
        hold on;
        if ~auto, if autokeep(currev), an2.String = 'keep'; else, an2.String = 'discard'; end, end
        plot(xpoints,ev_croptrc(currev,:),'o', 'color', trccol2); hold on;
        xline(xpoints(prepts+1), 'Linewidth',1.2, 'Color', [0 0 0]); alpha(.2); %xl.Color(4) = XLINEALPHA;
        xlim([xpoints(1) xpoints(end)]);
        xtsep = floor(100/diff(xpoints(1:2))); xtckmarks = xpoints(1:xtsep:end);
        xlabs = cell(numel(xtckmarks),1);
        for ii = 1:numel(xtckmarks), xlabs{ii} = num2str(xtckmarks(ii)); end
        xticks(xtckmarks); set(gca, 'XTickLabel', xlabs, 'fontsize', fontsize-3); xtickangle(45);
        xlabel('time [ms]'); ylabel('F/F');set(gca,'FontSize',8); 
        yrange = [min([ev_croptrc], [], 'all') max([ev_croptrc], [], 'all')];
        ylim([yrange(1)-diff(yrange)/10 yrange(2)+diff(yrange)/10]);
        hold off;
        histax = axes('Position', hist_pos);
        histogram(histax,rnd_peak_FoF,'BinEdges', binedges,'Normalization', 'Probability', 'FACEALPHA',.5, 'FaceColor',[.4 .61 .86]);
        hold on, histogram(histax,ctrl_peak_FoF,'BinEdges', binedges,'Normalization', 'Probability', 'FACEALPHA',.5, 'FaceColor',[.4 .86 .58]);
        hold on,histogram(histax, max(ev_croptrc(currev,:)),'BinEdges', binedges,'Normalization', 'Probability', 'FACEALPHA',.5,'FaceColor',[1 .57 .28]);
        leg=legend('rnd', 'neg','curr'); set(leg,'units','normalized','Box', false,'Orientation','vertical','Position',leg_pos, 'FontSize', 5);
        set(histax,'ytick',[]); set(histax,'ycolor',[1 1 1]); set(histax,'box','off');
        hold off;
    end

    function cb_backbut(~,~)
        if currev > 1
            currev = currev-1;
            keepev(currev) = false;
            plot_traces();
        end
    end
end