%% User Interfacce for ROIFlow Analysis for 2PM and Confocal data of spontaneous events %%

addpath('C:\Users\lena_\Projects\code_extern\Matlab\WriteImageJRoi\');

global ROILIST REVINFO LOADED PATH FILELIST SELFILE NFILES...
    NROIS SELROI FT CMAP roidata mainfig txtsm ddlist sp1a sp1b txtbig...
    fig3_pos FONTSIZE EVDATA DELROIS ALLROIS PARAMS DETECTED fig2_pos imsp2_pos...
    but2_po1_pos_1 trcsp2_pos but1_2 but1_3 but1_4 but1_5 but1_6 shownonev...
    TRCMODE PLOTMODE AUTO BATCH

%% Variables
ROILIST = {''}; REVINFO = {''}; FILELIST = []; DONELIST = cell(1,2); UNDONELIST = cell(1,2);
LOADED = false; PATH = [];
NFILES = 0; NROIS = 0;
roidata = [];
SELROI = []; SELFILE = []; DELROIS = []; ALLROIS=[];
FT = []; CMAP = [];
EVDATA = cell(1,12);
DETECTED = false;

shownonev = true;
PARAMS = get_predefined_params();
TRCMODE = 'dFoF';
PLOTMODE = 'aip';
AUTO = true;
BATCH = true;

%% GUI Parameters
close all;
FONTSIZE = 10;

scrsz = get(0,'ScreenSize'); scrsz=scrsz(3:4);
whfig1 = [round(4*scrsz(1)/12) scrsz(2)-50]; % Main (overview) window
whunits = 1./[7 15];
whddlist = [3 0.5].*whunits;
whbut1 = [2 0.7].*whunits;
whtxtbg = [3 2].*whunits;
whtxtsm = [5 1].*whunits;
whimsp1 = [5 8].*whunits;
whimsp1b = [0.5 4].*whunits;
hvspace1 = [0.5 0.5].*whunits;

whfig2 = [round(8*scrsz(1)/12) round(1.5*scrsz(2)/3)]; % Trace (overview) window
whunits = 1./[32 14];
whimsp2 = [4 4].*whunits;
whtrcsp2 = [23 10].*whunits;
whbut2 = [2 1].*whunits;
hvspace2 = [1 1].*whunits;

whfig3 = [round(4*scrsz(1)/12) round(1.5*scrsz(2)/3)]; % Revision window

txtbgcol = [0.2 0.2 0.2];
txtfgcol = [1 1 1];

%% GUI positions
fig1_pos = [10 10 whfig1];
but1_pos_1 = [5*hvspace1(1) hvspace1(2)*0.7 whbut1];
imsp1_pos = [hvspace1(1)*1.2 but1_pos_1(2)+but1_pos_1(4)+0.5*hvspace1(2) whimsp1];
imsp1b_pos = [imsp1_pos(1)+imsp1_pos(3) imsp1_pos(2)+whimsp1(2)/2-whimsp1b(2)/2 whimsp1b];
txtbig_pos = [hvspace1(1) imsp1_pos(2)+imsp1_pos(4)+1.5*hvspace1(2) whtxtbg];
ddlist_pos = [txtbig_pos(1) txtbig_pos(2)+txtbig_pos(4)+0.2*hvspace1(2) whddlist];
but1_pos_5 = [txtbig_pos(1)+txtbig_pos(3)+1.5*hvspace1(1) imsp1_pos(2)+imsp1_pos(4)+0.5*hvspace1(2) whbut1(1)*0.45 whbut1(2)*0.7];
but1_pos_6 = [but1_pos_5(1)+but1_pos_5(3)+0.2*hvspace1(1) but1_pos_5(2) whbut1(1)*0.45 whbut1(2)*0.7];
but1_pos_2 = [but1_pos_5(1) but1_pos_5(2)+but1_pos_5(4)+0.5*hvspace1(2) whbut1];
but1_pos_3 = [but1_pos_2(1) but1_pos_2(2)+but1_pos_2(4)+0.5*hvspace1(2) whbut1];
but1_pos_4 = [but1_pos_2(1) but1_pos_3(2)+but1_pos_3(4)+0.5*hvspace1(2) whbut1];
txtsm_pos = [hvspace1(1)*2 but1_pos_4(2)+but1_pos_4(4)+0.2*hvspace1(2) whtxtsm];
txttitle_pos = [hvspace1(1)*2.5 txtsm_pos(2)+0.5*txtsm_pos(4) whtxtsm];
rb_pos = [but1_pos_2(1)+but1_pos_2(3)+hvspace1(1)*0.2 but1_pos_2(2) whbut1(1)/4 whbut1(2)];
fig2_pos = [whfig1(1) round(scrsz(2)/2.25) whfig2];
imsp2_pos = [hvspace2(1)*1.5 5.5*hvspace2(2) whimsp2];
trcsp2_pos = [imsp2_pos(1)+imsp2_pos(3)+hvspace2(1)*2.5 hvspace2(2)*2.5 whtrcsp2];
but2_po1_pos_1 = [hvspace2(1)*2 hvspace2(2)*2.5 whbut2(1)*1.5 whbut2(2)];
fig3_pos = [round(fig2_pos(1)+whfig2(1)/10) 50 whfig3];

mainfig = figure('Name', 'ROIFlow - Automated Timeseries Analysis', 'Position', fig1_pos, 'toolbar', 'none', 'menu', 'none'); axis off;
set(mainfig,'WindowKeyPressFcn',@keyPressCallback);
sp1a = subplot('Position', imsp1_pos); axis('off');
sp1b = subplot('Position', imsp1b_pos); axis('off');
% trcfig = figure('Name', 'Trace overview', 'Position', fig2_pos, 'toolbar', 'none', 'menu', 'none'); axis off;
% sp2_1 = subplot('Position', imsp2_pos); axis('off');
% sp2_2 = subplot('Position', trcsp2_pos); axis('off');
% Add UI controls
ddlist = uicontrol('parent', mainfig, 'style', 'popupmenu','units','normalized','position', ddlist_pos,'FONTSIZE', FONTSIZE, 'string', ROILIST, 'Callback', {@cb_ddlist});
but1_1 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but1_pos_1,'string', 'Load/Select','FONTSIZE', FONTSIZE+2, 'callback', {@cb_but1_load});
but1_2 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but1_pos_2,'string', 'Accept & Save','FONTSIZE', FONTSIZE,'Enable','off', 'callback', {@cb_but1_accept});
but1_3 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but1_pos_3,'string', 'Revision (R)','FONTSIZE', FONTSIZE, 'Enable','off','callback', {@cb_but1_revision});
but1_4 =uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but1_pos_4,'string', 'Run Detection','FONTSIZE', FONTSIZE, 'Enable','off','callback', {@cb_but1_detection});
but1_5 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but1_pos_5,'string', '<<','FONTSIZE', FONTSIZE+1,'Enable','off', 'callback', {@cb_but56});
but1_6 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but1_pos_6,'string', '>>','FONTSIZE', FONTSIZE+1,'Enable','off', 'callback', {@cb_but56});
txtbig = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txtbig_pos,'string', '','BackgroundColor', txtbgcol,'ForegroundColor', txtfgcol, 'FONTSIZE', FONTSIZE);
txttitle = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txttitle_pos,'string', 'ROIFlow Event detection','FONTSIZE', round(FONTSIZE*1.3));
txtsm = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txtsm_pos,'string', '','FONTSIZE', FONTSIZE-2);
rb = uicontrol('parent', mainfig, 'style', 'radiobutton','units', 'normalized','position', rb_pos,'string', 'All','FONTSIZE', FONTSIZE-1,'Value', BATCH, 'callback', {@cb_rb_batch});

uiwait;

close all;
clear all;

%% Local Callback functions
function open_trcfig()
global  fig2_pos trcfig sp2_1 sp2_2 imsp2_pos but2_po1_pos_1 but2_po1 trcsp2_pos...
    FONTSIZE LOADED SELROI DETECTED PLOTMODE TRCMODE
if ishandle(trcfig), figure(trcfig);
else, trcfig = figure('Name', 'Trace overview', 'Position', fig2_pos, 'toolbar', 'none', 'menu', 'none'); axis off;
end
    set(trcfig,'WindowKeyPressFcn',@keyPressCallback);
sp2_1 = subplot('Position', imsp2_pos); axis('off');
sp2_2 = subplot('Position', trcsp2_pos); axis('off');
but2_po1 = uicontrol('parent', trcfig, 'style', 'pushbutton','units', 'normalized','position', but2_po1_pos_1,'string', 'Discard ROI','FONTSIZE', FONTSIZE,'BackgroundColor', [.4 .4 .4], 'callback', {@cb_but2_del_roi});

if LOADED
    plot_rois(SELROI, PLOTMODE,sp2_1,[]);
    plot_traces(SELROI, TRCMODE, DETECTED(SELROI));
end
end

function cb_but1_load(~,~)
global LOADED FILELIST SELFILE NFILES PATH FONTSIZE...
    but1_4 but1_5 but1_6 shownonev
if LOADED
    filesfig = figure('toolbar', 'none', 'menu', 'none'); axis('off');
    fl1 = uicontrol('parent', filesfig, 'style', 'listbox','units', 'normalized','position', [.05 .05 .4 .8],'string', FILELIST([FILELIST{:,3}]==0,1),'FONTSIZE', FONTSIZE, 'callback', {@cb_selfile});
    uicontrol('parent', filesfig, 'style', 'text','units', 'normalized','position', [.05 .85 .4 .1],'string', 'NOT REVISED','FONTSIZE', FONTSIZE+2);
    fl2 = uicontrol('parent', filesfig, 'style', 'listbox','units', 'normalized','position', [.55 .05 .4 .8],'string', FILELIST([FILELIST{:,3}]==1,1),'FONTSIZE', FONTSIZE, 'callback', {@cb_selfile});
    uicontrol('parent', filesfig, 'style', 'text','units', 'normalized','position', [.55 .85 .4 .1],'string', 'ACCEPTED','FONTSIZE', FONTSIZE+2);
    pter = false(NFILES,1); pter(SELFILE) = true;
    if FILELIST{SELFILE,3}==1, pter([FILELIST{:,3}]==0) = []; fl2.Value = find(pter==1);
    else, pter([FILELIST{:,3}]==1) = []; fl1.Value = find(pter==1);
    end
else
    if shownonev, [filenames, PATH] = uigetfile('*ROIDATA*.mat', 'Select processed .mat files with ROI data','Multiselect', 'on');
    else, [filenames, PATH] = uigetfile('ROIDATA*.mat', 'Select processed .mat files with ROI data','Multiselect', 'on');
    end
    if isa(filenames,'char'), filenames = {filenames}; end
    NFILES = numel(filenames);
    for i = 1:NFILES
        FILELIST{i,1} = filenames{i};
        FILELIST{i,2} = i;
        FILELIST{i,3} = false; % Codes if file has been processed
    end
    SELFILE = 1;
    LOADED = true;
    but1_4.Enable='on';
    but1_5.Enable='on';
    but1_6.Enable='on';
    load_mat();
end
end

function load_mat()
global PATH FILELIST SELFILE SELROI ROILIST EVDATA...
    NROIS DETECTED FT CMAP roidata txtsm ddlist sp1a sp1b ALLROIS PLOTMODE
load(fullfile(PATH,FILELIST{SELFILE,1}));
NROIS = roidata.n_rois;
SELROI = 1;
FT = roidata.frametime_s;
ROILIST = cell(NROIS,2);
ALLROIS = cell(NROIS,3);
DETECTED = false(NROIS,1);
for i = 1:NROIS
    ROILIST{i,1} = sprintf('ROI %i',i); ROILIST{i,2} = i;
    ALLROIS{i,3} = true;
end
ALLROIS(:,1:2) = ROILIST;
EVDATA = cell(NROIS,10);
txtsm.String = FILELIST{SELFILE,1};
ddlist.String = ROILIST(:,1);
ddlist.Value = 1;
pause(0.1);
cmap_n = max(cellfun(@max ,roidata.roi_bounds(:,2)));
CMAP = get_colormap([1 0 0],[1 1 0],[0 1 1],cmap_n);
% CMAP = get_colormap([1 0 0],[1 1 0],[0 1 1],NROIS);
plot_rois(1:NROIS, PLOTMODE, sp1a, sp1b);
open_trcfig();
end

function cb_but1_accept(~,~)
global SELFILE FILELIST EVDATA roidata...
    PATH trcfig BATCH
if ishandle(trcfig), figure(trcfig); close gcf; end
tmpfl = find([FILELIST{:,3}]== 0);
if BATCH
    % Process all files with chosen setting
    process_batch(tmpfl);
else
    [~] = save_event_info(EVDATA, roidata, PATH, FILELIST{SELFILE},false);
    FILELIST{SELFILE,3} = true;
    % Work on next file
    if ~isempty(tmpfl)
        SELFILE = FILELIST{tmpfl(1),2};
        load_mat();
    else, uiresume; close all;
    end
end
end

function cb_rb_batch(hOb,~)
global BATCH
BATCH = hOb.Value;
end

function cb_but1_revision(~,~)
global fig3_pos sp2_2 sp1a sp1b ddlist EVDATA DETECTED ROILIST...
    roidata SELROI FT TRCMODE PLOTMODE
if ~strcmp(EVDATA{SELROI,11}, 'perfect') && ~isempty(EVDATA{SELROI,6})
    goOn = false;
    while ~goOn
        [goOn, keep_rev] = revise_events(fig3_pos, EVDATA(SELROI,:), roidata.dFoF_traces(SELROI,:), FT);
    end
    if ~isempty(keep_rev)
        tmp_rev = EVDATA{SELROI,6}(1:numel(keep_rev));
        EVDATA{SELROI,7} = [EVDATA{SELROI,7} tmp_rev(keep_rev)];
        EVDATA{SELROI,6} = EVDATA{SELROI,6}(numel(keep_rev)+1:end);
        if ishandle(sp2_2), plot_traces(SELROI, TRCMODE, DETECTED(SELROI));
        else, open_trcfig();
        end
        if DETECTED(SELROI), update_main_fig();end
        clear('tmp_rev');
    end
    % Discard Roi from list if no events left
    if isempty(EVDATA{SELROI,7})
        currpter = find(([ROILIST{:,2}]==SELROI)==1);
        ROILIST(currpter,:) = [];
        if ~isempty(ROILIST), ddlist.Value = currpter; SELROI = ROILIST{currpter,2};
        else, ddlist.String = 'NO ROIs';ddlist.Enable='off'; SELROI = 0;
        end
        plot_rois([ROILIST{:,2}], PLOTMODE, sp1a, sp1b);
    end
end
end

function cb_but1_detection(~,~)
global roidata EVDATA ALLROIS PARAMS NROIS SELROI...
    sp2_2 sp1a sp1b but1_3 but1_2 trcfig ddlist DETECTED ROILIST TRCMODE PLOTMODE FT AUTO
% Detects events in all traces
if SELROI == 0, ROILIST = ALLROIS(:,1:2); SELROI = 1; end
[EVDATA, detroi, PARAMS] = run_event_detection(EVDATA, ROILIST, roidata.dFoF_traces, roidata.FoF_traces, roidata.traces, PARAMS, FT, AUTO, true);
DETECTED(detroi) = true;
but1_3.Enable='on';
but1_2.Enable='on';
% Discard Rois from list if they dont show events
if any(logical(cellfun(@isempty, EVDATA(:,7))) &...
        logical(cellfun(@isempty, EVDATA(:,6))))
    pterdel = logical(cellfun(@isempty, EVDATA(:,7)))&logical(cellfun(@isempty, EVDATA(:,6)));
    pternew = find(pterdel==0 & [ALLROIS{:,3}]');
    delrois = find(pterdel==1);
    ROILIST = ALLROIS(pternew,:);
    if ~isempty(ddlist.Value) && any(ddlist.Value == delrois) && ~isempty(pternew)
        ddlist.Value = 1;
        SELROI = ROILIST{1,2};
    end
else
    ROILIST = ALLROIS([ALLROIS{:,3}]',:);
    pternew=1:NROIS;
end
if isempty(pternew), ddlist.String = 'NO ROIs';ddlist.Enable='off'; SELROI = 0; close(trcfig);
elseif pternew~=0, ddlist.Enable='on';ddlist.String = ROILIST(:,1);
end
if SELROI ~= 0
    if ishandle(sp2_2), plot_traces(SELROI, TRCMODE, DETECTED(SELROI));
    else, open_trcfig();
    end
if DETECTED(SELROI), update_main_fig();end
end
plot_rois(pternew, PLOTMODE, sp1a, sp1b);
end

function cb_ddlist(hOb,~)
global SELROI sp2_1 sp2_2 DETECTED ROILIST TRCMODE PLOTMODE
SELROI = ROILIST{hOb.Value,2};
if ishandle(sp2_2)
    plot_traces(SELROI, TRCMODE, DETECTED(SELROI));
    plot_rois(SELROI, PLOTMODE,sp2_1,[]);
else, open_trcfig();
end
if DETECTED(SELROI), update_main_fig(); end
end

function cb_but56(hOb,~)
if strcmp(hOb.String,'<<'), swdir = -1;
else, swdir = 1; end
switch_trace(swdir);
end

function keyPressCallback(~,ev)
if strcmp(ev.Key,'rightarrow'), swdir = 1;switch_trace(swdir);
elseif strcmp(ev.Key,'leftarrow'), swdir = -1;switch_trace(swdir);
elseif strcmp(ev.Key,'r'), cb_but1_revision([],[]);
end
end

function switch_trace(swdir)
global SELROI sp2_1 sp2_2 DETECTED EVDATA txtbig ROILIST...
    TRCMODE ddlist PLOTMODE
if SELROI ~= 0
    tmpid = find(SELROI == [ROILIST{:,2}]);
    if swdir == -1
        if tmpid > 1, SELROI = ROILIST{tmpid-1,2}; tmpid=tmpid-1; end
    elseif swdir == 1
        if tmpid < size(ROILIST,1), SELROI = ROILIST{tmpid+1,2}; tmpid=tmpid+1; end
    end
    ddlist.Value = tmpid;
    if ishandle(sp2_2)
        plot_traces(SELROI, TRCMODE, DETECTED(SELROI));
        plot_rois(SELROI, PLOTMODE,sp2_1,[]);
    else, open_trcfig();
    end
    if DETECTED(SELROI)
        update_main_fig();
    end
end
end

function cb_but2_del_roi(~,~)
global ALLROIS SELROI ROILIST ddlist sp1a sp1b DETECTED TRCMODE PLOTMODE
% Discard Roi
ALLROIS{SELROI,3} = false;
currpter = find(([ROILIST{:,2}]==SELROI)==1);
ROILIST(currpter,:) = [];
if ~isempty(ROILIST)
    SELROI = ROILIST{currpter,2}; ddlist.String = ROILIST(:,1); ddlist.Value = currpter;
    plot_traces(SELROI, TRCMODE, DETECTED(SELROI));
    if DETECTED(SELROI), update_main_fig(); end
else,SELROI=0; ddlist.String = 'NO ROIs';ddlist.Enable='off';
end
plot_rois([ROILIST{2,:}], PLOTMODE, sp1a,sp1b);
end

function cb_selfile(hOb,~)
global SELFILE DETECTED FILELIST
tmpsel = hOb.Value;
if DETECTED, quest = questdlg('Discard current detection information?', 'Change file', 'Yes', 'No', 'No');
else, quest = 'Yes';
end
if strcmp(quest,'Yes')
    if FILELIST{strcmp(hOb.String(hOb.Value),FILELIST(:,1)),3} == 1
        tmpacc = find([FILELIST{:,3}]==1);
        FILELIST{SELFILE,3} = false;
        SELFILE = FILELIST{tmpacc(tmpsel),1};
    else
        tmprev = find([FILELIST{:,3}]==0);
        FILELIST{SELFILE,3} = false;
        SELFILE = FILELIST{tmprev(tmpsel),2};
    end
    close gcf;
    load_mat();
end
end

function plot_rois(plotids, mode, spa, spb)
global roidata CMAP
axes(spa); cla;
if strcmp(mode,'aip')
    uint8_im = (roidata.aip./max(roidata.aip,[],'all'))*255;
    axis([0 size(uint8_im,2) 0 size(uint8_im,1)]);
    imagesc(uint8(uint8_im));colormap('gray');axis('off');daspect([1 1 1]);
elseif strcmp(mode,'bg')
    imsize = size(roidata.aip);
    rim = zeros(imsize(1), imsize(2),1);
    gim = zeros(imsize(1), imsize(2),1);
    bim = zeros(imsize(1), imsize(2),1);
    bg_col = [0.1529 0.1882 0.6706]; % RGB
    rim(background==1) = bg_col(1);
    gim(background==1) = bg_col(2);
    bim(background==1) = bg_col(3);
    rgbim = cat(3, rim, gim,bim);
    imshow(rgbim);
end

for cc = 1:numel(plotids)
    B = roidata.roi_bounds{plotids(cc),1};
    for k = 1:length(B)
        boundary = B{k};
        hold on; plot(boundary(:,2), boundary(:,1), 'Color',CMAP(:,roidata.roi_bounds{plotids(cc),2})', 'LineWidth', 1)
    end
end

if ~isempty(spb)
    axes(spb); cla;
    v1 = min(cellfun(@min, roidata.roi_bounds(:,2)));
    v2 = max(cellfun(@max, roidata.roi_bounds(:,2)));
    colormap(spb, CMAP(:,v1:v2)');
    cb = colorbar();
    steps = ceil(v2-v1/10);
    cbtl = linspace(v1, v2, steps)';
    tmplim = get(cb,'YLim');
    cbt = linspace(tmplim(1), tmplim(2), steps);
    set(cb, 'YTick', cbt, 'YTickLabel',cbtl, 'Location', 'eastoutside');
end

end

function plot_traces(roiid, type, events)
global roidata FT FONTSIZE EVDATA SELROI sp2_2
axes(sp2_2);
cla;
trccol = [0 0 0]; %[0.6510 0.0667 0.2784];
if strcmp(type,'dFoF'), tmptraces = roidata.dFoF_traces; yaxlab = 'dFoF +1';
elseif strcmp(type,'FoF'), tmptraces = roidata.FoF_traces; yaxlab = 'FoF';
else, tmptraces = roidata.traces; yaxlab = 'Intensity [a.u.]';
end
plot(FT.*(1:size(tmptraces,2)), tmptraces(roiid,:)+1, 'Color', trccol, 'LineWidth', 1);
if ~isempty(EVDATA{SELROI,10})
    hold on, yline(1+EVDATA{SELROI,10}(2)*EVDATA{SELROI,10}(3), 'cyan');
    text(1,1+EVDATA{SELROI,10}(2)*EVDATA{SELROI,10}(3), strcat('det.thresh. x ',num2str(EVDATA{SELROI,10}(3))), 'Fontsize', 6);
end
if ~isempty(EVDATA{SELROI,10})
    yline(1+EVDATA{SELROI,10}(1)*EVDATA{SELROI,10}(4), 'red');
    text(1,1+EVDATA{SELROI,10}(1)*EVDATA{SELROI,10}(4), strcat('%-tile val. x ',num2str(EVDATA{SELROI,10}(4))), 'Fontsize', 6);
end
if ~isempty(EVDATA{SELROI,10})
    hold on, yline(1-EVDATA{SELROI,10}(1), 'Color',[.5 .5 .5]);
    yline(1+EVDATA{SELROI,10}(1), 'Color',[.5 .5 .5]);
end
yrange = [min(tmptraces(roiid,:))+1 max(tmptraces(roiid,:))+1];
ylim([yrange(1)-diff(yrange)/10 yrange(2)+diff(yrange)/10]);
xlim([FT size(tmptraces(roiid,:),2)*FT]);
xlabel('Time [s]'); ylabel(yaxlab);
set(gca, 'FontSize', FONTSIZE-1);
if events
    save_pts = EVDATA{SELROI,7}.*FT;
    rev_pts = EVDATA{SELROI,6}.*FT;
    hold on;
    for iEv = 1:numel(save_pts), xline(save_pts(iEv), 'Color',[0 0 0], 'LineWidth', 1.2, 'Alpha',.4); end
    for iEv = 1:numel(rev_pts), xline(rev_pts(iEv), 'r', 'LineWidth', 1.2, 'Alpha', .4); end
    hold off;
end
end

function update_main_fig()
global txtbig but1_3 EVDATA SELROI NROIS
txtbig.String = {'',sprintf('Accepted ROIs: %i / %i', sum(strcmp([EVDATA(:,11)],'perfect')), NROIS), '',...
    sprintf('State of ROI: %s', EVDATA{SELROI,11}),'',...
    sprintf('Accepted: %i | Revision: %i', numel(EVDATA{SELROI,7}), numel(EVDATA{SELROI,6})), '',...
    sprintf('Detect: > %ixSD | Save: > %ixSD',EVDATA{SELROI,10}(1), EVDATA{SELROI,10}(2))};
pause(0.1);
if strcmp(EVDATA{SELROI,11}, 'perfect'), but1_3.Enable = 'off';
else, but1_3.Enable = 'on';
end
end

function process_batch(file_idx)
global PATH FILELIST PARAMS
nfiles = numel(file_idx);
exp_summary = struct('Rec_Name', [], 'Num_ROIs',[],'ToT_Num_Events',[],'Amp_Mean',[],'Amp_FoF_Mean',[],'Amp_dFoF_Mean',[],...
    'Amp_SD',[],'Amp_FoF_SD',[],'Amp_dFoF_SD',[], 'IEI_Mean',[],'EvRate_Mean',[],...
    'Area_um_Mean',[],'Dend_Area_um',[],'Bg_Area_um',[],'FoV_Area_um',[],'StructNorm_Area',[],'Events_per_Dend_Area',[]);
sum_cnt = 1;
for iFile = 1:nfiles
        fprintf('Processing file %i / %i\n', iFile, nfiles);
        tmpfile = FILELIST{file_idx(iFile),2};
        evdata = cell(1,12);
        
        % Load
        load(fullfile(PATH,FILELIST{tmpfile,1}));
        nrois = roidata.n_rois;
        roilist = cell(nrois,2);
        for i = 1:nrois
            roilist{i,1} = sprintf('ROI %i',i); roilist{i,2} = i;
        end
        cmap_n = max(cellfun(@max ,roidata.roi_bounds(:,2)));
        cmap = get_colormap([1 0 0],[1 1 0],[0 1 1],cmap_n);
        % cmap = get_colormap([1 0 0],[1 1 0],[0 1 1],NROIS);
        
        % Detect
        [evdata, detroi, PARAMS] = run_event_detection(evdata, roilist, roidata.dFoF_traces, roidata.FoF_traces, roidata.traces, PARAMS, roidata.frametime_s, true, false);
        
        % Discard Rois from list if they dont show events
        if any(logical(cellfun(@isempty, evdata(:,7))) &...
                logical(cellfun(@isempty, evdata(:,6))))
            keep = find(logical(cellfun(@isempty, evdata(:,7)))&logical(cellfun(@isempty, evdata(:,6))) == 0);
            roilist = roilist(keep,:); 
        end
        
        % Save
        [tmp_summary] = save_event_info(evdata, roidata, PATH, FILELIST{tmpfile}, true);
%         FILELIST{tmpfile,3} = true;

        % Write all rercording summaries into master table
        if ~isempty(tmp_summary.Num_ROIs)
            exp_summary(sum_cnt).Rec_Name = FILELIST{tmpfile}(1:end-4);
            exp_summary(sum_cnt).Num_ROIs = tmp_summary.Num_ROIs;
            exp_summary(sum_cnt).ToT_Num_Events = tmp_summary.ToT_Num_Events;
            exp_summary(sum_cnt).Num_Events_Mean = tmp_summary.Num_Events_Mean;
            exp_summary(sum_cnt).Amp_Mean = tmp_summary.Amp_Mean;
            exp_summary(sum_cnt).Amp_FoF_Mean = tmp_summary.Amp_FoF_Mean;
            exp_summary(sum_cnt).Amp_dFoF_Mean = tmp_summary.Amp_dFoF_Mean;
            exp_summary(sum_cnt).Amp_SD = tmp_summary.Amp_SD;
            exp_summary(sum_cnt).Amp_FoF_SD = tmp_summary.Amp_FoF_SD;
            exp_summary(sum_cnt).Amp_dFoF_SD = tmp_summary.Amp_dFoF_SD;
            exp_summary(sum_cnt).Mean_SNR = tmp_summary.Mean_SNR;
            exp_summary(sum_cnt).Event_confidence = tmp_summary.Event_confidence;
            exp_summary(sum_cnt).IEI_Mean = tmp_summary.IEI_Mean;
            exp_summary(sum_cnt).EvRate_Mean = tmp_summary.EvRate_Mean;
            exp_summary(sum_cnt).REC_time_s = tmp_summary.REC_time_s;
            exp_summary(sum_cnt).Area_um_Mean = tmp_summary.Area_um_Mean;
            exp_summary(sum_cnt).StructNorm_Area = tmp_summary.StructNorm_Area;
            exp_summary(sum_cnt).Dend_Area_um = tmp_summary.Dend_Area_um;
            exp_summary(sum_cnt).Bg_Area_um = tmp_summary.Bg_Area_um;
            exp_summary(sum_cnt).FoV_Area_um = tmp_summary.FoV_Area_um;
            exp_summary(sum_cnt).Events_per_Dend_Area = tmp_summary.Events_per_Dend_Area;
            sum_cnt = sum_cnt+1;
        end
end

% Write master table to file
savepath = fullfile(PATH,FILELIST{tmpfile}(1:8));
exp_summary_xls = strcat(savepath, '_Master.xlsx'); if isa(exp_summary_xls,'cell'), exp_summary_xls=exp_summary_xls{1};end
writetable(struct2table(exp_summary), exp_summary_xls);

msgbox('All files processed!');
uiresume; close all;
end