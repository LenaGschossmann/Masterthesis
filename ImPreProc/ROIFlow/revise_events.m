function [goOn, keepev] = revise_events(fig3_pos, roi_evdata, trace, ft)

%% Parameters
save_idx = roi_evdata{1,7};
rev_idx = roi_evdata{1,6};
rand_amps = roi_evdata{1,8};
neg_amps = roi_evdata{1,9};
amps = roi_evdata{1,5}(:,1);
num_save = numel(save_idx);
num_rev = numel(rev_idx);
prepts = round(0.4/ft);
postpts = round(0.8/ft);
totpts = prepts+postpts+1;
save_croptrc = NaN(num_save,totpts);
rev_croptrc = NaN(num_save,totpts);
currev = 1;
keepev = false(num_rev,1);
xpoints = (1:totpts).*ft.*1000;
goOn = [];
binedges = [min([amps; rand_amps; neg_amps]) max([amps; rand_amps; neg_amps])];
binedges = linspace(binedges(1),binedges(2),20);

%% GUI
whfig3 = fig3_pos(3:4);
whunits = 1./[12 12];
whbut3 = [2 1].*whunits;
% whtxt = [];
whtrcsp = [10 7].*whunits;
whhistsp = [3 2].*whunits;
hvspace3 = [1 1].*whunits;

trccol1 = [0.6510 0.0667 0.2784];
trccol2 = [0 0 0];
fontsize = 10;
an=[];

but3_pos_1 = [hvspace3(1)*3 hvspace3(2) whbut3];
but3_pos_2 = [but3_pos_1(1)+but3_pos_1(3)+0.5*hvspace3(1) hvspace3(2) whbut3];
but3_pos_3 = [but3_pos_2(1)+but3_pos_2(3)+0.5*hvspace3(1) hvspace3(2) whbut3];
trcsp_pos = [hvspace3(1)*1.25 but3_pos_1(2)+but3_pos_1(4)+2*hvspace3(2) whtrcsp];
hist_pos = [trcsp_pos(1)+trcsp_pos(3)*0.65 trcsp_pos(2)+0.68*trcsp_pos(4) whhistsp];

revfig = figure('Name', 'Event revision', 'Position', fig3_pos, 'toolbar', 'none', 'menu', 'none'); axis off;
set(revfig,'WindowKeyPressFcn',@keyPressCallback);
but3_1 = uicontrol('parent', revfig, 'style', 'pushbutton','units', 'normalized','position', but3_pos_1,'string', 'Accept(A)','fontsize', fontsize-1, 'BackgroundColor', '#77AC30','callback', {@cb_but3_accept});
but3_2 = uicontrol('parent', revfig, 'style', 'pushbutton','units', 'normalized','position', but3_pos_2,'string', 'Reject (D)','fontsize', fontsize-1, 'BackgroundColor', '#A2142F','callback', {@cb_but3_reject});
but3_3 = uicontrol('parent', revfig, 'style', 'pushbutton','units', 'normalized','position', but3_pos_3,'string', 'Back (B)','fontsize', fontsize-1, 'callback', {@cb_but3_back});

%% Extract cropped out region around events
% Save events
for iE = 1:numel(save_idx)
    % Event preceding part
    if save_idx(iE)-prepts < 1, save_croptrc(iE,prepts-save_idx(iE)+2:prepts+1) = trace(1:save_idx(iE))';
    else, save_croptrc(iE,1:prepts+1) = trace(save_idx(iE)-prepts:save_idx(iE))';
    end
    % Event following part
    if save_idx(iE)+postpts > numel(trace), save_croptrc(iE,prepts+2:prepts+1+numel(trace)-save_idx(iE)) = trace(save_idx(iE)+1:end);
    else, save_croptrc(iE,prepts+2:end) = trace(save_idx(iE)+1:save_idx(iE)+postpts);
    end
end
% Revision events
for iE = 1:numel(rev_idx)
    % Event preceding part
    if rev_idx(iE)-prepts < 1, rev_croptrc(iE,prepts-rev_idx(iE)+2:prepts+1) = trace(1:rev_idx(iE))';
    else, rev_croptrc(iE,1:prepts+1) = trace(rev_idx(iE)-prepts:rev_idx(iE))';
    end
    % Event following part
    if rev_idx(iE)+postpts > numel(trace), rev_croptrc(iE,prepts+2:prepts+1+numel(trace)-rev_idx(iE)) = trace(rev_idx(iE)+1:end);
    else, rev_croptrc(iE,prepts+2:end) = trace(rev_idx(iE)+1:rev_idx(iE)+postpts);
    end
end

plot_traces();

%% Apply revision
uiwait();
if isempty(goOn)
    goOn = true;
    keepev(currev:end) = [];
end

%% Local functions
    function keyPressCallback(~,ev)
        if strcmp(ev.Key,'a')
            keepev(currev) = true;
            if currev < num_rev, currev = currev+1; plot_traces();
            else, pause(0.1); goOn = true; close gcf;
            end
        elseif strcmp(ev.Key,'d')
            keepev(currev) = false;
            if currev < num_rev, currev = currev+1; plot_traces();
            else, pause(0.1); goOn = true; close gcf;
            end
        elseif strcmp(ev.Key,'b')
            if currev > 1
                currev = currev-1;
                keepev(currev) = false;
                plot_traces();
            end
        end
    end

    function cb_but3_accept(~,~)
        keepev(currev) = true;
        if currev < num_rev, currev = currev+1; plot_traces();
        else, pause(0.1); goOn = true; close gcf;
        end
    end

    function cb_but3_reject(~,~)
        keepev(currev) = false;
        if currev < num_rev, currev = currev+1; plot_traces();
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
        if isempty(an), an = annotation('textbox', [.85 .8 .1 .1], 'FitBoxToText','on','EdgeColor', [1 1 1],'fontsize', fontsize); end
        for iE = 1:num_save
            p = plot(xpoints,save_croptrc(iE,:),'linewidth', 1, 'color', trccol1); p.Color(4) = .25;
            hold on;
        end
        if currev > 1
            plotidx = find(keepev(1:currev-1) == true);
            for iE = 1:numel(plotidx)
                p = plot(xpoints,rev_croptrc(plotidx(iE),:),'linewidth', 1, 'color', trccol1); p.Color(4) = .25;
                hold on;
            end
        end
        
        % Current trace
        if ~all(isnan(rev_croptrc(currev,:)))
            hold on;
            plot(xpoints,rev_croptrc(currev,:),'linewidth', 1.2, 'color', trccol2);
            an.String = sprintf('%i/%i', currev,num_rev);
        end
        hold on;
        plot(xpoints,rev_croptrc(currev,:),'o', 'color', trccol2); hold on;
        xline(xpoints(prepts+1), 'Linewidth',1.2, 'Color', [0 0 0]); alpha(.2); %xl.Color(4) = XLINEALPHA;
        xlim([xpoints(1) xpoints(end)]);
        xtsep = floor(100/diff(xpoints(1:2))); xtckmarks = xpoints(1:xtsep:end);
        xlabs = cell(numel(xtckmarks),1);
        for ii = 1:numel(xtckmarks), xlabs{ii} = num2str(xtckmarks(ii)); end
        xticks(xtckmarks); set(gca, 'XTickLabel', xlabs, 'fontsize', fontsize-3); xtickangle(45);
        xlabel('time [ms]'); ylabel('dF/F');set(gca,'FontSize',8); 
        yrange = [min([rev_croptrc;save_croptrc], [], 'all') max([rev_croptrc;save_croptrc], [], 'all')];
        ylim([yrange(1)-diff(yrange)/10 yrange(2)+diff(yrange)/10]);

        histax = axes('Position', hist_pos);
        histogram(histax,rand_amps,'BinEdges', binedges,'Normalization', 'Probability', 'FACEALPHA',.5, 'FaceColor',[.4 .61 .86]);
        hold on, histogram(histax,neg_amps,'BinEdges', binedges,'Normalization', 'Probability', 'FACEALPHA',.5, 'FaceColor',[.4 .86 .58]);
        hold on,histogram(histax, rev_croptrc(currev,prepts+1),'BinEdges', binedges,'Normalization', 'Probability', 'FACEALPHA',.5,'FaceColor',[1 .57 .28]);
        set(histax,'ytick',[]); set(histax,'ycolor',[1 1 1]); set(histax,'box','off');
    end

    function cb_backbut(~,~)
        if currev > 1
            currev = currev-1;
            keepev(currev) = false;
            plot_traces();
        end
    end
end