%% Enter...
% filepointer = 'C:\Users\lena_\Projects\iGluSnFR\Dev_V2\cropped_sources\1-200_20210115_7030_sl_02_20x_TTX_ts_07.tif';
% frametime = 0.12288; % [s]
% px_bg_yx = [77 80; 25 28; 94 14]; % 1
% px_stna_yx = [31 102; 27 63; 42 3]; % 2
% px_sta_yx = [69 30; 55 49; 73 49]; % 3

filepointer = 'C:\Users\lena_\Projects\iGluSnFR\Dev_V2\cropped_sources\800-1000_20210114_7035_sl_02_20x_TTX_ts_01.tif';
frametime = 0.06144; % [s]
px_bg_yx = [91 7; 22 48; 3 79]; % 1
px_stna_yx = [63 56; 30 65; 31 17]; % 2
px_sta_yx = [74 13; 72 18; 74 17]; % 3
ev = 63;

px_all_yx = [px_bg_yx; px_stna_yx; px_sta_yx];
px_all_id = [repelem(1,size(px_bg_yx,1))'; repelem(2,size(px_stna_yx,1))'; repelem(3,size(px_sta_yx,1))'];

% Add paths
wd = pwd(); p = split(wd,'\'); p = p(1:end-1); pp = []; for ii=1:numel(p), pp = fullfile(pp,p{ii}); end
addpath(genpath(pp)), clear('p','pp');

mode = 'roll';
output = 'dFoF';
winsizes = [5 10 15]; % [frames]
shift = 5;

tmpim = read_timeseries(filepointer, 'load');
tmpim = double(tmpim);

imsize = size(tmpim);
n_px = imsize(1)*imsize(2);
n_frames = imsize(3);
lintraces = reshape(tmpim,[n_px n_frames]);

lindx = sub2ind([imsize(1) imsize(2)], px_all_yx(:,1),px_all_yx(:,2));
usetraces = double(lintraces(lindx,:));

baseline_vals = cell(numel(winsizes),1);

for iWin = 1:numel(winsizes)
    tmpwin = winsizes(iWin);
    
    avtraces = usetraces;
    avsteps = NaN(n_frames,2);
    avsteps(1:tmpwin+shift,1) = 1; avsteps(1:tmpwin+shift,2) = tmpwin;
    avsteps(tmpwin+shift+1:n_frames,1) = 2:n_frames-tmpwin-shift+1;
    avsteps(tmpwin+shift+1:n_frames,2) = tmpwin+1:n_frames-shift;
    
    for iAv = 1:n_frames
        avtraces(:,iAv) = mean(usetraces(:,avsteps(iAv,1):avsteps(iAv,2)),2);
    end
    
    baseline_vals{iWin} = avtraces;
    
end

differences_1 = baseline_vals{1}-baseline_vals{2};
differences_2 = baseline_vals{1}-baseline_vals{3};
differences_3 = baseline_vals{2}-baseline_vals{3};

lims = NaN(numel(px_all_id),2);
edges = cell(numel(px_all_id),1);
for ii = 1:numel(px_all_id)
    m1 = min([min(differences_1(ii,:)) min(differences_2(ii,:)) min(differences_3(ii,:))]);
    m2 = max([max(differences_1(ii,:)) max(differences_2(ii,:)) max(differences_3(ii,:))]);
    lims(ii,:) = [m1 m2];
    edges{ii} = linspace(round(m1), round(m2), 20);
end

cnt = 1;
f_1=figure(); 
subplot(3,3,1); histogram(differences_1(cnt,:),edges{cnt}); title(sprintf('BL %i - BL %i', winsizes(1),winsizes(2)));
subplot(3,3,2); histogram(differences_2(cnt,:),edges{cnt}, 'FaceColor', 'red'); title(sprintf('Diff. BL %i - %i', winsizes(1),winsizes(3)));
subplot(3,3,3); histogram(differences_3(cnt,:),edges{cnt}, 'FaceColor', 'green'); title(sprintf('Diff. BL %i - %i', winsizes(2),winsizes(3)));
cnt = cnt+1;
subplot(3,3,4); histogram(differences_1(cnt,:),edges{cnt});
subplot(3,3,5); histogram(differences_2(cnt,:),edges{cnt}, 'FaceColor', 'red');
subplot(3,3,6); histogram(differences_3(cnt,:),edges{cnt}, 'FaceColor', 'green');
cnt = cnt+1;
subplot(3,3,7); histogram(differences_1(cnt,:),edges{cnt});
subplot(3,3,8); histogram(differences_2(cnt,:),edges{cnt}, 'FaceColor', 'red');
subplot(3,3,9); histogram(differences_3(cnt,:),edges{cnt}, 'FaceColor', 'green');
sgtitle('Baseline differences - Background');

f_2=figure(); 
cnt = cnt+1;
subplot(3,3,1); histogram(differences_1(cnt,:),edges{cnt}); title(sprintf('BL %i - BL %i', winsizes(1),winsizes(2))); ylabel(sprintf('Px %i',cnt));
subplot(3,3,2); histogram(differences_2(cnt,:),edges{cnt}, 'FaceColor', 'red'); title(sprintf('BL %i - BL %i', winsizes(1),winsizes(3)));
subplot(3,3,3); histogram(differences_3(cnt,:),edges{cnt}, 'FaceColor', 'green'); title(sprintf('BL %i - BL %i', winsizes(2),winsizes(3)));
cnt = cnt+1;
subplot(3,3,4); histogram(differences_1(cnt,:),edges{cnt}); ylabel(sprintf('Px %i',cnt));
subplot(3,3,5); histogram(differences_2(cnt,:),edges{cnt}, 'FaceColor', 'red');
subplot(3,3,6); histogram(differences_3(cnt,:),edges{cnt}, 'FaceColor', 'green');
cnt = cnt+1;
subplot(3,3,7); histogram(differences_1(cnt,:),edges{cnt}); ylabel(sprintf('Px %i',cnt));
subplot(3,3,8); histogram(differences_2(cnt,:),edges{cnt}, 'FaceColor', 'red');
subplot(3,3,9); histogram(differences_3(cnt,:),edges{cnt}, 'FaceColor', 'green');
sgtitle('Baseline differences - Structure/No Activity');

f_3=figure(); 
cnt = cnt+1;
subplot(3,3,1); histogram(differences_1(cnt,:),edges{cnt}); title(sprintf('BL %i - BL %i', winsizes(1),winsizes(2))); ylabel(sprintf('Px %i',cnt));
subplot(3,3,2); histogram(differences_2(cnt,:),edges{cnt}, 'FaceColor', 'red'); title(sprintf('BL %i - BL %i', winsizes(1),winsizes(3)));
subplot(3,3,3); histogram(differences_3(cnt,:),edges{cnt}, 'FaceColor', 'green'); title(sprintf('BL %i - BL %i', winsizes(2),winsizes(3)));
cnt = cnt+1;
subplot(3,3,4); histogram(differences_1(cnt,:),edges{cnt}); ylabel(sprintf('Px %i',cnt));
subplot(3,3,5); histogram(differences_2(cnt,:),edges{cnt}, 'FaceColor', 'red');
subplot(3,3,6); histogram(differences_3(cnt,:),edges{cnt}, 'FaceColor', 'green');
cnt = cnt+1;
subplot(3,3,7); histogram(differences_1(cnt,:),edges{cnt}); ylabel(sprintf('Px %i',cnt));
subplot(3,3,8); histogram(differences_2(cnt,:),edges{cnt}, 'FaceColor', 'red');
subplot(3,3,9); histogram(differences_3(cnt,:),edges{cnt}, 'FaceColor', 'green');
sgtitle('Baseline differences - Structure/Activity');

%% *****
lims = NaN(numel(px_all_id),2);
edges = cell(numel(px_all_id),1);
for ii = 1:numel(px_all_id)
    m1 = min([min(baseline_vals{1}(ii,:)) min(baseline_vals{2}(ii,:)) min(baseline_vals{3}(ii,:))]);
    m2 = max([max(baseline_vals{1}(ii,:)) max(baseline_vals{2}(ii,:)) max(baseline_vals{3}(ii,:))]);
    lims(ii,:) = [m1 m2];
    edges{ii} = linspace(round(m1), round(m2), 20);
end

cnt = 1;
f_bl_1=figure();
subplot(3,3,1); histogram(baseline_vals{1}(cnt,:),edges{cnt}); title(sprintf('BL %i', winsizes(1))); ylabel(sprintf('Px %i',cnt));
subplot(3,3,2); histogram(baseline_vals{2}(cnt,:),edges{cnt}, 'FaceColor', 'red'); title(sprintf('BL %i', winsizes(2)));
subplot(3,3,3); histogram(baseline_vals{3}(cnt,:),edges{cnt}, 'FaceColor', 'green'); title(sprintf('BL %i', winsizes(3)));
cnt = cnt+1;
subplot(3,3,4); histogram(baseline_vals{1}(cnt,:),edges{cnt}); ylabel(sprintf('Px %i',cnt));
subplot(3,3,5); histogram(baseline_vals{2}(cnt,:),edges{cnt}, 'FaceColor', 'red');
subplot(3,3,6); histogram(baseline_vals{3}(cnt,:),edges{cnt}, 'FaceColor', 'green');
cnt = cnt+1;
subplot(3,3,7); histogram(baseline_vals{1}(cnt,:),edges{cnt}); ylabel(sprintf('Px %i',cnt));
subplot(3,3,8); histogram(baseline_vals{2}(cnt,:),edges{cnt}, 'FaceColor', 'red');
subplot(3,3,9); histogram(baseline_vals{3}(cnt,:),edges{cnt}, 'FaceColor', 'green');
sgtitle('Baseline (BL) values - Background');

f_bl_2=figure();
cnt = cnt+1;
subplot(3,3,1); histogram(baseline_vals{1}(cnt,:),edges{cnt}); title(sprintf('BL %i', winsizes(1))); ylabel(sprintf('Px %i',cnt));
subplot(3,3,2); histogram(baseline_vals{2}(cnt,:),edges{cnt}, 'FaceColor', 'red'); title(sprintf('BL %i', winsizes(2)));
subplot(3,3,3); histogram(baseline_vals{3}(cnt,:),edges{cnt}, 'FaceColor', 'green'); title(sprintf('BL %i', winsizes(3)));
cnt = cnt+1;
subplot(3,3,4); histogram(baseline_vals{1}(cnt,:),edges{cnt}); ylabel(sprintf('Px %i',cnt));
subplot(3,3,5); histogram(baseline_vals{2}(cnt,:),edges{cnt}, 'FaceColor', 'red');
subplot(3,3,6); histogram(baseline_vals{3}(cnt,:),edges{cnt}, 'FaceColor', 'green');
cnt = cnt+1;
subplot(3,3,7); histogram(baseline_vals{1}(cnt,:),edges{cnt}); ylabel(sprintf('Px %i',cnt));
subplot(3,3,8); histogram(baseline_vals{2}(cnt,:),edges{cnt}, 'FaceColor', 'red');
subplot(3,3,9); histogram(baseline_vals{3}(cnt,:),edges{cnt}, 'FaceColor', 'green');
sgtitle('Baseline (BL) values - Structure/No Activity');

f_bl_3=figure(); 
cnt = cnt+1;
subplot(3,3,1); histogram(baseline_vals{1}(cnt,:),edges{cnt}); title(sprintf('BL %i', winsizes(1))); ylabel(sprintf('Px %i',cnt));
subplot(3,3,2); histogram(baseline_vals{2}(cnt,:),edges{cnt}, 'FaceColor', 'red'); title(sprintf('BL %i', winsizes(2)));
subplot(3,3,3); histogram(baseline_vals{3}(cnt,:),edges{cnt}, 'FaceColor', 'green'); title(sprintf('BL %i', winsizes(3)));
cnt = cnt+1;
subplot(3,3,4); histogram(baseline_vals{1}(cnt,:),edges{cnt}); ylabel(sprintf('Px %i',cnt));
subplot(3,3,5); histogram(baseline_vals{2}(cnt,:),edges{cnt}, 'FaceColor', 'red');
subplot(3,3,6); histogram(baseline_vals{3}(cnt,:),edges{cnt}, 'FaceColor', 'green');
cnt = cnt+1;
subplot(3,3,7); histogram(baseline_vals{1}(cnt,:),edges{cnt}); ylabel(sprintf('Px %i',cnt));
subplot(3,3,8); histogram(baseline_vals{2}(cnt,:),edges{cnt}, 'FaceColor', 'red');
subplot(3,3,9); histogram(baseline_vals{3}(cnt,:),edges{cnt}, 'FaceColor', 'green');
sgtitle('Baseline (BL) values - Structure/Activity');

%% *** Plot F/F traces
pxplot = [7 8 9]; %1 2 3    4 5 6   7 8 9
foftraces_1 = usetraces(pxplot,:)./baseline_vals{1}(pxplot,:);
foftraces_2 = usetraces(pxplot,:)./baseline_vals{2}(pxplot,:);
foftraces_3 = usetraces(pxplot,:)./baseline_vals{3}(pxplot,:);

figure();
px_yx = px_all_yx(pxplot,:);
cnt=1;
subplot(3,1,cnt); plot(foftraces_1(cnt,:),'k'); ylabel(sprintf('y %i, x %i',px_yx(cnt,1),px_yx(cnt,2)));
hold on, plot(foftraces_2(cnt,:),'r');
hold on, plot(foftraces_3(cnt,:),'b');
xline(ev, 'g');
legend(sprintf('size %i',winsizes(1)),sprintf('size %i',winsizes(2)),sprintf('size %i',winsizes(3)), 'Event');
xlim([1 size(foftraces_1,2)]);
cnt=cnt+1;
subplot(3,1,cnt); plot(foftraces_1(cnt,:),'k'); ylabel(sprintf('y %i, x %i',px_yx(cnt,1),px_yx(cnt,2)));
hold on, plot(foftraces_2(cnt,:),'r');
hold on, plot(foftraces_3(cnt,:),'b');
xlim([1 size(foftraces_1,2)]);
xline(ev, 'g');
cnt=cnt+1;
subplot(3,1,cnt); plot(foftraces_1(cnt,:),'k'); ylabel(sprintf('y %i, x %i',px_yx(cnt,1),px_yx(cnt,2)));
hold on, plot(foftraces_2(cnt,:),'r');
hold on, plot(foftraces_3(cnt,:),'b');
xlim([1 size(foftraces_1,2)]);
xline(ev, 'g');
sgtitle('F/F traces - Structure/Activity'); %Background  Structure/No activity  Structure/Activity

%% *** Baseline
pxplot = [7 8 9]; %1 2 3    4 5 6   7 8 9
plottraces_1 = baseline_vals{1}(pxplot,:);
plottraces_2 = baseline_vals{2}(pxplot,:);
plottraces_3 = baseline_vals{3}(pxplot,:);

figure();
px_yx = px_all_yx(pxplot,:);
cnt=1;
subplot(3,1,cnt); plot(usetraces(pxplot(cnt),:),'Color', [.5 .5 .5]); ylabel(sprintf('y %i, x %i',px_yx(cnt,1),px_yx(cnt,2)));
hold on, plot(plottraces_1(cnt,:),'k'); 
hold on, plot(plottraces_2(cnt,:),'r');
hold on, plot(plottraces_3(cnt,:),'b');
legend('raw trace', sprintf('size %i',winsizes(1)),sprintf('size %i',winsizes(2)),sprintf('size %i',winsizes(3)));
xlim([1 size(plottraces_1,2)]);
cnt=cnt+1;
subplot(3,1,cnt); plot(usetraces(pxplot(cnt),:),'Color', [.5 .5 .5]); ylabel(sprintf('y %i, x %i',px_yx(cnt,1),px_yx(cnt,2))); 
hold on, plot(plottraces_1(cnt,:),'k');
hold on, plot(plottraces_2(cnt,:),'r');
hold on, plot(plottraces_3(cnt,:),'b');
xlim([1 size(plottraces_1,2)]);
cnt=cnt+1;
subplot(3,1,cnt); plot(usetraces(pxplot(cnt),:),'Color', [.5 .5 .5]); ylabel(sprintf('y %i, x %i',px_yx(cnt,1),px_yx(cnt,2))); ylabel(sprintf('y %i, x %i',px_yx(cnt,1),px_yx(cnt,2)));
hold on, plot(plottraces_1(cnt,:),'k'); 
hold on, plot(plottraces_2(cnt,:),'r');
hold on, plot(plottraces_3(cnt,:),'b');
xlim([1 size(plottraces_1,2)]);
sgtitle('Baseline values - Structure/Activity'); %Background  Structure/No activity  Structure/Activity
