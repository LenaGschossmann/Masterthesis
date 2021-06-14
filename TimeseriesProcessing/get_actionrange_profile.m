function [] = get_actionrange_profile()
% dF/F values are scaled to 16bit!
addpath('C:\Users\lena_\Projects\code_extern\Matlab\gramm-master');

defpath = 'D:\Masterthesis stuff\Analysis\Confocal\20201117_A895_CAG.FLEX.iGluSnFR.A184V\20201117_A895_sl_04_63x_ts_01';
[fn,p] = uigetfile({'*.csv'; '*.xlsx'},'Select file', defpath);
seglength = inputdlg('Enter length of profile line in µm'); seglength = str2num(seglength{1});
data = readmatrix(fullfile(p,fn));
datainfo = data(2:end,1:2);
dataxprofile_ori = data(1,3:end);
datayprofile = data(2:end,3:end);
scppx = seglength/size(dataxprofile_ori,2);

% figure();
% for ii=1:size(datainfo,1)
%     plot(dataxprofile, yprofile(ii,:),'b-');
%     hold on;
% end

% % Smooth y data
% win = round(size(dataxprofile,2)/10); hwin = floor(win/2);
% ii = hwin+1; % startposition
% ndata = size(dataxprofile,2);
% smoothyprofile = datayprofile;
% 
% while ii+hwin <= ndata
%     smoothyprofile(:,ii) = mean(datayprofile(:,ii-hwin:ii+hwin),2);
%     ii=ii+1;
% end
% 
% yprofile = smoothyprofile(:,hwin+1:ndata-hwin);
% dataxprofile = dataxprofile(hwin+1:ndata-hwin);

yprofile = datayprofile;
dataxprofile = dataxprofile_ori;

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

%% Fit gauss curves
for ii = 1: size(str.x,1)
    Y = str.y(ii,:);
    X = str.x(ii,:);
    fitmdl = fit(X',Y','gauss1');
    coeffs = [fitmdl.a1 fitmdl.b1 fitmdl.c1]; 
    str.fit(ii,:) = coeffs(1).*exp(-((X-coeffs(2))./coeffs(3)).^2);
end


%% Center on max and find FWHM
npts = size(str.x,2);
str.fit_ctr = NaN(size(str.x,1), npts*2);
fwhm = NaN(size(str.x,1),1);
ctridx=1:npts;
for ii = 1: size(str.x,1)
    [valmax, idxmax] = max(str.fit(ii,:));
    newidx = ctridx-idxmax;
    str.fit_ctr(ii,npts+newidx) = str.fit(ii,:);
    minval = min(str.fit(ii,:));
    hmax = ((valmax-minval)/2)+minval;
    tmp = str.fit(ii,:)>=hmax;
    tmpx1 = find(tmp==1,1,'first'); tmpx2 = find(tmp==1,1,'last');
    fwhm(ii) = (tmpx2-tmpx1)*scppx;
end

str.fit_ctr(:,all(isnan(str.fit_ctr),1)) = [];
str.x_fit_ctr = (1:size(str.fit_ctr,2)) - size(str.fit_ctr,2)/2;


%% Save fit and traces
tbl_out = table();
tbl_out.type = str.id;
tbl_out.fwhm_um = fwhm;
outputname = strcat(p,fn(1:end-4), '_fwhm.xlsx');
writetable(tbl_out, outputname)

tbl_out = table();
tbl_out.xvalues = str.x'.*scppx;
tbl_out.yvalues = str.y';
outputname = strcat(p,fn(1:end-4), '_fit.xlsx');
writetable(tbl_out, outputname, 'Sheet',1);
tbl_out = table();
tbl_out.xvalues = str.x'.*scppx;
tbl_out.fit = str.fit';
writetable(tbl_out, outputname, 'Sheet',2);

tbl_out = table();
tbl_out.xvalues = str.x_fit_ctr'.*scppx;
tbl_out.yvalues = str.fit_ctr';
outputname = strcat(p,fn(1:end-4), '_fit_ctr.xlsx');
writetable(tbl_out, outputname);

%% Plotting
% Colormap
% cmhex = {'#878787','#9e0142','#c2a5cf','#f46d43'};
cmhex = {'#878787','#9e0142','#c2a5cf','6dcd59ff'};
% cmrgb = [0 0 0; 128 222 44
cm=[];
for ii = 1:numel(cmhex), cm = [cm; sscanf(cmhex{ii}(2:end),'%2x%2x%2x',[1 3])/255]; end

% Raw Intensity
roi = 'ROI 1';
nm = strcat(p,'ActionRange_',roi,'_raw_smImJ.png'); % Outputname
ylab = 'Raw intensity [a.u.]';
tit = sprintf('Action range profile - %s - smoothed',roi);
f=figure();
g = gramm('x',str.x_fit_ctr, 'y',str.fit_ctr, 'color', str.id);
% g = gramm('x',str.xall, 'y',str.yall, 'color', str.idall);
g.geom_line();
g.set_color_options('map', cm);
g.set_names('color', 'Frames', 'x', 'Distance [um]', 'y', ylab);
g.set_title(tit);
g.draw();
% saveas(f,nm);

% % dF
% nm = strcat(p,'ActionRange_',roi,'_dF_smImJ.png'); % Outputname
% ylab = 'dF';
% tit = sprintf('Action range profile - %s - smoothed',roi);
% f=figure();
% g = gramm('x',str.x, 'y',str.dFy, 'color', str.id);
% g.geom_line();
% g.set_color_options('map', cm(2:end,:));
% g.set_names('color', 'Frames', 'x', 'Distance [um]', 'y', ylab);
% g.set_title(tit);
% g.draw();
% saveas(f,nm);
end