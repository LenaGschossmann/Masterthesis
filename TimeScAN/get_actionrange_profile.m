function [] = get_actionrange_profile()

addpath('C:\Users\lena_\Projects\code_extern\Matlab\gramm-master');

defpath = 'D:\Masterthesis stuff\Analysis\Confocal\20201117_A895_CAG.FLEX.iGluSnFR.A184V\20201117_A895_sl_04_63x_ts_01';
[f,p] = uigetfile({'*.csv'; '*.xlsx'},'Select file', defpath);
data = readmatrix(fullfile(p,f));
datainfo = data(2:end,1:2);
dataxprofile = data(1,3:end);
datayprofile = data(2:end,3:end);

% figure();
% for ii=1:size(datainfo,1)
%     plot(dataxprofile, yprofile(ii,:),'b-');
%     hold on;
% end

% Smooth y data
win = 5; hwin = floor(win/2);
ii = hwin+1; % startposition
ndata = size(dataxprofile,2);
smoothyprofile = datayprofile;
while ii+hwin <= ndata
    smoothyprofile(:,ii) = mean(datayprofile(:,ii-hwin:ii+hwin),2);
    ii=ii+1;
end

yprofile = datayprofile;

% Baseline profile
baseidx = find(datainfo(:,1)== -1);
baseprofile = mean(yprofile(baseidx,:),1);
% Prepeak profile
preidx = find(datainfo(:,1)== 0);
preprofile = yprofile(preidx,:);
% Peak profile
peakidx = find(datainfo(:,1)== 1);
peakprofile = yprofile(peakidx,:);
% Postpeak profile
postidx = find(datainfo(:,1)== 2);
postprofile = yprofile(postidx,:);

nf = numel(preidx) + numel(peakidx) + numel(postidx);

% Subtract baseline
dFpreprofile = preprofile./baseprofile;
dFpeakprofile = peakprofile./baseprofile;
dFpostprofile = postprofile./baseprofile;

str = struct();
str.id = [repelem({'pre'},numel(preidx))'; 'peak'; repelem({'post'},numel(postidx))'];
str.x = repmat(dataxprofile,[nf 1]);
str.xall = repmat(dataxprofile,[nf+1 1]);
str.y = yprofile([preidx; peakidx;postidx],:);
str.dFy = [dFpreprofile; dFpeakprofile; dFpostprofile];
str.basey = baseprofile;
str.yall = [baseprofile;yprofile([preidx; peakidx;postidx],:)];
str.idall = ['baseline';repelem({'pre'},numel(preidx))'; 'peak'; repelem({'post'},numel(postidx))'];

% Plotting
% Colormap
cmhex = {'#878787','#9e0142','#c2a5cf','#f46d43'};
cm=[];
for ii = 1:numel(cmhex), cm = [cm; sscanf(cmhex{ii}(2:end),'%2x%2x%2x',[1 3])/255]; end

% Raw Intensity
roi = 'ROI 1';
nm = strcat(p,'ActionRange_',roi,'_raw_smImJ.png'); % Outputname
ylab = 'Raw intensity [a.u.]';
tit = sprintf('Action range profile - %s - smoothed',roi);
f=figure();
g = gramm('x',str.xall, 'y',str.yall, 'color', str.idall);
g.geom_line();
g.set_color_options('map', cm);
g.set_names('color', 'Frames', 'x', 'Distance [um]', 'y', ylab);
g.set_title(tit);
g.draw();
saveas(f,nm);

% dF
nm = strcat(p,'ActionRange_',roi,'_dF_smImJ.png'); % Outputname
ylab = 'dF';
tit = sprintf('Action range profile - %s - smoothed',roi);
f=figure();
g = gramm('x',str.x, 'y',str.dFy, 'color', str.id);
g.geom_line();
g.set_color_options('map', cm(2:end,:));
g.set_names('color', 'Frames', 'x', 'Distance [um]', 'y', ylab);
g.set_title(tit);
g.draw();
saveas(f,nm);
end