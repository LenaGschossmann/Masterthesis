function create_summary()

% Declare globally shared variables
global POSITIONTRCFIG FONTSIZE TRCPLOTCOL1...
   SMTHWIN FTIME TRACEDATA RANGE evINFO LW1 LW2...
   XLINEALPHA PREPOINTS POSTPOINTS...
   TRCALPHA

%% Initialize display parameters
fs = FONTSIZE+1;
nbins = 20;

overlaypos = [0.07 0.42 0.4 0.5];
risepos = [0.09 0.11 0.15 0.2];
decaypos = [0.27 0.11 0.15 0.2];
amphistpos = [0.51 0.7 0.25 0.25];
ampbarpos = [0.8 0.7 0.17 0.25];
ieihistpos = [0.51 0.35 0.25 0.25];
ieibarpos = [0.8 0.35 0.17 0.25];
cvieipos = [0.51 0.11 0.15 0.15];

summaryfig = figure('Position',POSITIONTRCFIG, 'Name', 'Event summary');

%% Summary data
% dFoF traces
traceidx = find(strcmp({evINFO(:).accepted},'accepted') == 1);
numtr = numel(traceidx); numallev = 0;
for iTr = 1:numtr, numallev = numallev + numel(evINFO(traceidx(iTr)).onsetidx); end
worktraces = TRACEDATA(RANGE(1):RANGE(2), traceidx);
if SMTHWIN ~= 0, for iTr = 1:numtr, worktraces(:,iTr) = smooth_data(worktraces(:,iTr),SMTHWIN); end, end
croptrace = NaN(numallev,PREPOINTS+POSTPOINTS+1);
xpoints = (0:PREPOINTS+POSTPOINTS).*FTIME.*1000; % in [ms]
rise = zeros(numallev,1); decay = zeros(numallev,1); amp = zeros(numallev,1);
iei = []; cviei = zeros(numtr,2);
ii = 1;
% Extrace cropped out region around events
for iTr = 1:numtr
    % Initialize und unpack varaibles
    numev = numel(evINFO(traceidx(iTr)).onsetidx);
    onsetidx = evINFO(traceidx(iTr)).onsetidx;
    baseline = evINFO(traceidx(iTr)).baselinevalues;
    for iE = 1:numev
        % Event preceding part
        if onsetidx(iE)-PREPOINTS < 1
            croptrace(ii,PREPOINTS-onsetidx(iE)+2:PREPOINTS+1) = worktraces(1:onsetidx(iE), iTr)';
        else
            croptrace(ii,1:PREPOINTS+1) = worktraces(onsetidx(iE)-PREPOINTS:onsetidx(iE), iTr)';
        end
        % Event following part
        if onsetidx(iE)+POSTPOINTS > size(worktraces,1)
            croptrace(ii,PREPOINTS+2:PREPOINTS+1+size(worktraces,1)-onsetidx(iE)) = worktraces(onsetidx(iE)+1:end,iTr);
        else
            croptrace(ii,PREPOINTS+2:end) = worktraces(onsetidx(iE)+1:onsetidx(iE)+POSTPOINTS, iTr);
        end
        % Delta F over F
        croptrace(ii,:) = (croptrace(ii,:)-baseline(iE))./baseline(iE);
        
        % Extract parameter
        rise(ii)  = evINFO(traceidx(iTr)).risetimes(iE);
        decay(ii) = evINFO(traceidx(iTr)).decaytimes(iE);
        amp(ii) = evINFO(traceidx(iTr)).amps_dFoF(iE);
        if iE < numev, iei = [iei evINFO(traceidx(iTr)).ieis(iE)]; end
        ii = ii+1;
    end
    cviei(iTr) = evINFO(traceidx(iTr)).cviei;
end

%% Plots
sp1 = subplot('Position', overlaypos);
sp2a = subplot('Position', risepos);
sp2b = subplot('Position', decaypos);
sp3 = subplot('Position', amphistpos);
sp4 = subplot('Position', ampbarpos);
sp5 = subplot('Position', ieihistpos);
sp6 = subplot('Position', ieibarpos);
sp7 = subplot('Position', cvieipos);

axes(sp2a);
boxplot(rise.*1000, 'Boxstyle', 'filled', 'colors', TRCPLOTCOL1, 'Symbol', 'ok'); % in [ms]
set(gca,'FontSize', FONTSIZE-1);set(gca,'TitleFontSizeMultiplier', 1.25);
set(gca, 'XTickLabel',''); ylabel('time [ms]'); title('Risetime');

axes(sp2a);
boxplot(decay.*1000, 'Boxstyle', 'filled', 'colors', TRCPLOTCOL1, 'Symbol', 'ok'); % in [ms]
set(gca,'FontSize', FONTSIZE-1);set(gca,'TitleFontSizeMultiplier', 1.25);
set(gca, 'XTickLabel',''); ylabel('time [ms]'); title('Decaytime');

axes(sp3);
histogram(amp,nbins, 'FaceColor', TRCPLOTCOL1);
set(gca,'FontSize', FONTSIZE-1); set(gca,'TitleFontSizeMultiplier', 1.25);
title('dFoF amplitude'); xlabel('dFoF amplitudes');

axes(sp4);
boxplot(amp, 'Boxstyle', 'filled', 'colors', TRCPLOTCOL1, 'Symbol', 'ok'); % in [ms]
set(gca,'FontSize', FONTSIZE-1); set(gca,'TitleFontSizeMultiplier', 1.25);
ylabel('dFoF');set(gca, 'XTickLabel','');title('dF/F amplitude');

axes(sp5);
histogram(iei,nbins, 'FaceColor', TRCPLOTCOL1);
set(gca,'FontSize', FONTSIZE-1); set(gca,'TitleFontSizeMultiplier', 1.25);
title('Inter-Event-Intervals'); xlabel('intervals [s]');

axes(sp6);
boxplot(iei, 'Boxstyle', 'filled', 'colors', TRCPLOTCOL1, 'Symbol', 'ok'); % in [s]
set(gca,'FontSize', FONTSIZE-1);set(gca,'TitleFontSizeMultiplier', 1.25);
ylabel('interval [s]');
set(gca, 'XTickLabel',''); title('Inter-Event-Intervals');

axes(sp7);
boxplot(cviei, 'Boxstyle', 'filled', 'colors', TRCPLOTCOL1);
set(gca,'FontSize', FONTSIZE-1); set(gca,'TitleFontSizeMultiplier', 1.25);
title('IEI CV');

axes(sp1);
for iE = 1:numallev
    p=plot(xpoints,croptrace(iE,:),'linewidth', LW1, 'color', TRCPLOTCOL1); p.Color(4) = TRCALPHA;
    hold on;
end
plot(xpoints,mean(croptrace,1),'linewidth', LW2, 'color', [0 0 0]);
xl=xline(xpoints(PREPOINTS+1), 'Linewidth',LW2+1, 'Color', [0 0 0]); xl.Color(4) = 0.9; %XLINEALPHA;
xlim([xpoints(1) xpoints(end)]);
xtckmarks = xpoints(1:5:end);
xlabs = cell(numel(xtckmarks),1);
for ii = 1:numel(xtckmarks), xlabs{ii} = num2str(xtckmarks(ii)); end
xticks(xtckmarks); set(gca, 'XTickLabel', xlabs, 'FontSize', FONTSIZE-1); xtickangle(45);
xlabel('time [ms]'); ylabel('dF/F');
yrange = [min(croptrace, [], 'all') max(croptrace, [], 'all')];
ylim([yrange(1)-diff(yrange)/10 yrange(2)+diff(yrange)/10]);
set(gca,'TitleFontSizeMultiplier', 1.25); title('Aligned events');
        
        
end