%% Enter...
% filepointer = 'C:\Users\lena_\Projects\iGluSnFR\Dev_V2\cropped_sources\1-200_20210115_7030_sl_02_20x_TTX_ts_07.tif';
% frametime = 0.12288; % [s]
% px_bg_yx = [77 80; 25 28; 94 14; 20 28; 23 19; 89 98; 96 53; 89 50; 94 9; 97 9]; % 1
% px_stna_yx = [31 102; 27 63; 42 3; 15 42; 13 29; 13 95; 14 59; 25 6; 32 38; 26 57]; % 2
% px_sta_yx = [69 30; 55 49; 73 49; 54 49; 55 49; 50 44; 50 43; 49 42; 50 40; 57 35]; % 3

% filepointer = 'C:\Users\lena_\Projects\iGluSnFR\Dev_V2\cropped_sources\800-1000_20210114_7035_sl_02_20x_TTX_ts_01.tif';
% frametime = 0.06144; % [s]
% px_bg_yx = [91 7; 22 48; 3 79; 13 80; 4 26; 94 3; 11 81; 61 69; 25 53; 13 21]; % 1
% px_stna_yx = [63 56; 30 65; 31 17;45 79; 50 37; 45 86; 47 40; 23 75; 40 14; 12 50]; % 2
% px_sta_yx = [74 13; 72 18; 74 17; 71 16; 71 15; 71 14; 70 16; 69 17; 69 16; 74 20]; % 3
% ev = 63;

% Uncropped
% filepointer = 'C:\Users\lena_\Projects\iGluSnFR\ValidationTests\1_Man_ROIs\source\20210115_7030_sl_02_20x_TTX_ts_07.lsm';
% frametime = 0.12288; % [s]
% px_bg_yx = [136 109; 177 112; 177 111; 177 110; 184 108; 184 101; 184 100; 221 14; 226 14; 244 15]; % 1
% px_stna_yx = [29 186; 26 179; 25 217; 33 152; 18 157; 18 217; 69 173; 75 147; 115 10; 104 9]; % 2
% px_sta_yx = [202 133; 201 127; 200 127; 213 139; 222 131; 222 130; 222 129; 219 119; 205 106; 211 120]; % 3

% filepointer = 'C:\Users\lena_\Projects\iGluSnFR\ValidationTests\1_Man_ROIs\source\20210114_7035_sl_03_20x_TTX_ts_01.lsm';
% frametime = 0.06144; % [s]
% px_bg_yx = [47 49; 47 50; 49 46; 46 46; 45 45; 18 4; 9 4; 11 2; 18 3; 1 14; 1 8]; % 1
% px_stna_yx = [39 44; 39 45; 36 50; 35 50; 46 58; 46 57; 96 66; 101 56; 101 55; 95 47]; % 2
% px_sta_yx = [9 57; 8 59; 10 60; 10 59; 10 58; 10 57; 11 57; 11 58; 11 59; 5 56]; % 3
% ev = 63;

% 2PM
filepointer = 'C:\Users\lena_\Projects\iGluSnFR\ValidationTests\1_Man_ROIs\source\20210105_PCP2R4WT_7029_iGluSnFR.A184S_tFrame_12a.nd2';
frametime = 0.0512; % [s]
px_bg_yx = [3 407; 3 406; 3 405; 9 320; 7 327; 96 468; 40 33; 46 28; 62 49; 61 46]; % 1
px_stna_yx = [25 479; 47 316; 24 245;26 259;117 127; 68 22; 38 213;75 187; 59 281; 44 326]; % 2
px_sta_yx = [65 302; 65 301; 62 308; 62 309; 58 301; 58 302; 56 293; 56 292; 56 294; 69 304]; % 3

px_all_yx = [px_bg_yx; px_stna_yx; px_sta_yx]; %px_sta_yx
px_all_yx = px_all_yx+1;
px_all_id = [repelem(1,size(px_bg_yx,1))'; repelem(2,size(px_stna_yx,1))'; repelem(3,size(px_sta_yx,1))']; %repelem(3,size(px_sta_yx,1))'

% Add paths
wd = pwd(); p = split(wd,'\'); p = p(1:end-1); pp = []; for ii=1:numel(p), pp = fullfile(pp,p{ii}); end
addpath(genpath(pp)), clear('p','pp');

tmpim = read_timeseries(filepointer, 'load');
tmpim = double(tmpim);

imsize = size(tmpim);
n_px = imsize(1)*imsize(2);
n_frames = imsize(3);
lintraces = reshape(tmpim,[n_px n_frames]);

lindx = sub2ind([imsize(1) imsize(2)], px_all_yx(:,1),px_all_yx(:,2));
usetraces = double(lintraces(lindx,:));
sd_traces = std(usetraces,[],2);
% av_traces = mean(usetraces,2); sd_traces = sd_traces./av_traces;
sd_px_rec = std(lintraces,[],2);
av_px_rec = mean(lintraces,2);
z_sc_px_rec = sd_px_rec./av_px_rec;
% test_prctile = prctile(sd_px_rec,40,'all');
range = [prctile(sd_px_rec,5,'all') prctile(sd_px_rec,95,'all')];
cutoff = range(1) + diff(range)*0.15
range(1) + diff(range)*0.2
range(1) + diff(range)*0.25
range(1) + diff(range)*0.3

%% Histograms of intensities
edges = linspace(min(usetraces,[],'all'), max(usetraces,[],'all'), 20);
figure();
subplot(2,3,1); histogram(usetraces(1:10,:),edges, 'FaceColor', 'red'); title('Background px');
subplot(2,3,2); histogram(usetraces(11:20,:),edges); title('Structure px');
subplot(2,3,3); histogram(usetraces(21:30,:),edges,'FaceColor', 'green'); title('Structure Activity px');
edges_sd = linspace(min(sd_traces,[],'all'), max(sd_traces,[],'all'), 5);
subplot(2,3,4); histogram(sd_traces(1:10,:),edges_sd, 'FaceColor', 'red'); title('SD Background px');
subplot(2,3,5); histogram(sd_traces(11:20,:),edges_sd); title('SD Structure px');
subplot(2,3,6); histogram(sd_traces(21:30,:),edges_sd,'FaceColor', 'green'); title('SD Structure Activity px');

%%
figure();
subplot(3,1,1); histogram(av_px_rec); title('Mean per px, whole rec.');
subplot(3,1,2); histogram(sd_px_rec); title('SD per px, whole rec.');
subplot(3,1,3); histogram(z_sc_px_rec); title('SD/Mean per px, whole rec.');

%%
figure, plot(sd_traces(1:10), mean(usetraces(1:10,:),2),'*r');
hold on, plot(sd_traces(11:20), mean(usetraces(11:20,:),2),'*b');
hold on, plot(sd_traces(21:30), mean(usetraces(21:30,:),2),'*k');
xlim([edges_sd(1)-10 edges_sd(end)+10]);xlabel('SD'); ylabel('Av. intensity');
legend('Background', 'Structure/No activity', 'Structure/Activity');

%% *** Plot
pxplot = [4 5 6]; %1 2 3    4 5 6   7 8 9
ylimits = [min(usetraces,[],'all')-20 max(usetraces,[],'all')];

figure();
px_yx = px_all_yx(pxplot,:);
cnt=1;
subplot(3,1,cnt); plot(usetraces(pxplot(cnt),:)); ylabel(sprintf('y %i, x %i',px_yx(cnt,1),px_yx(cnt,2)));
xlim([1 size(usetraces,2)]); ylim(ylimits);
legend(strcat('SD: ', num2str(round(std(usetraces(pxplot(cnt),:))))));
cnt=cnt+1;
subplot(3,1,cnt); plot(usetraces(pxplot(cnt),:)); ylabel(sprintf('y %i, x %i',px_yx(cnt,1),px_yx(cnt,2)));
xlim([1 size(usetraces,2)]); ylim(ylimits);
legend(strcat('SD: ', num2str(round(std(usetraces(pxplot(cnt),:))))));
cnt=cnt+1;
subplot(3,1,cnt); plot(usetraces(pxplot(cnt),:)); ylabel(sprintf('y %i, x %i',px_yx(cnt,1),px_yx(cnt,2)));
xlim([1 size(usetraces,2)]); ylim(ylimits);
legend(strcat('SD: ', num2str(round(std(usetraces(pxplot(cnt),:))))));
sgtitle('Raw Px traces - Structure/No activity'); %Background  Structure/No activity  Structure/Activity


