%% User Interfacce for ROIFlow Analysis for 2PM and Confocal data of spontaneous events %%

global ROILIST REVINFO LOADED PATH FILELIST SELFILE NFILES...
    NROIS SELROI FT CMAP roidata mainfig txtsm ddlist sp1 txtbig...
    fig3_pos FONTSIZE EVDATA DELROIS ALLROIS PARAMS DETECTED fig2_pos imsp2_pos...
    but2_po1_pos_1 trcsp2_pos but1_2 but1_3 but1_4 but1_5 but1_6 shownonev...
    TRCMODE PLOTMODE

%% Variables
ROILIST = {''}; REVINFO = {''}; FILELIST = []; DONELIST = cell(1,2); UNDONELIST = cell(1,2);
LOADED = false; PATH = [];
NFILES = 0; NROIS = 0;
roidata = [];
SELROI = []; SELFILE = []; DELROIS = []; ALLROIS=[];
FT = []; CMAP = [];
EVDATA = cell(1,10);
DETECTED = false;

shownonev = true;
PARAMS = get_predefined_params();
TRCMODE = 'dFoF';
PLOTMODE = 'aip';

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
imsp1_pos = [hvspace1(1)*1.5 but1_pos_1(2)+but1_pos_1(4)+0.5*hvspace1(2) whimsp1];
txtbig_pos = [hvspace1(1) imsp1_pos(2)+imsp1_pos(4)+1.5*hvspace1(2) whtxtbg];
ddlist_pos = [txtbig_pos(1) txtbig_pos(2)+txtbig_pos(4)+0.2*hvspace1(2) whddlist];
but1_pos_5 = [but1_pos_2(1) imsp1_pos(2)+imsp1_pos(4)+0.5*hvspace1(2) whbut1(1)*0.45 whbut1(2)*0.7];
but1_pos_6 = [but1_pos_5(1)+but1_pos_5(3)+0.2*hvspace1(1) but1_pos_5(2) whbut1(1)*0.45 whbut1(2)*0.7];
but1_pos_2 = [txtbig_pos(1)+txtbig_pos(3)+2*hvspace1(1) but1_pos_5(2)+but1_pos_5(4)+0.5*hvspace1(2) whbut1];
but1_pos_3 = [but1_pos_2(1) but1_pos_2(2)+but1_pos_2(4)+0.5*hvspace1(2) whbut1];
but1_pos_4 = [but1_pos_2(1) but1_pos_3(2)+but1_pos_3(4)+0.5*hvspace1(2) whbut1];
txtsm_pos = [hvspace1(1)*2 but1_pos_4(2)+but1_pos_4(4)+0.2*hvspace1(2) whtxtsm];
txttitle_pos = [hvspace1(1)*2.5 txtsm_pos(2)+0.5*txtsm_pos(4) whtxtsm];
fig2_pos = [whfig1(1) round(scrsz(2)/2.25) whfig2];
imsp2_pos = [hvspace2(1)*1.5 5.5*hvspace2(2) whimsp2];
trcsp2_pos = [imsp2_pos(1)+imsp2_pos(3)+hvspace2(1)*2.5 hvspace2(2)*2.5 whtrcsp2];
but2_po1_pos_1 = [hvspace2(1)*2 hvspace2(2)*2.5 whbut2(1)*1.5 whbut2(2)];
fig3_pos = [round(fig2_pos(1)+whfig2(1)/10) 50 whfig3];

mainfig = figure('Name', 'ROIFlow - Automated Timeseries Analysis', 'Position', fig1_pos, 'toolbar', 'none', 'menu', 'none'); axis off;
set(mainfig,'WindowKeyPressFcn',@keyPressCallback);
sp1 = subplot('Position', imsp1_pos); axis('off');
% trcfig = figure('Name', 'Trace overview', 'Position', fig2_pos, 'toolbar', 'none', 'menu', 'none'); axis off;
% sp2_1 = subplot('Position', imsp2_pos); axis('off');
% sp2_2 = subplot('Position', trcsp2_pos); axis('off');
% Add UI controls
ddlist = uicontrol('parent', mainfig, 'style', 'popupmenu','units','normalized','position', ddlist_pos,'FONTSIZE', FONTSIZE, 'string', ROILIST, 'Callback', {@cb_ddlist});
but1_1 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but1_pos_1,'string', 'Load/Select','FONTSIZE', FONTSIZE+2, 'callback', {@cb_but1_load});
but1_2 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but1_pos_2,'string', 'Accept & Save','FONTSIZE', FONTSIZE,'Enable','off', 'callback', {@cb_but1_accept});
but1_3 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but1_pos_3,'string', 'Revision','FONTSIZE', FONTSIZE, 'Enable','off','callback', {@cb_but1_revision});
but1_4 =uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but1_pos_4,'string', 'Run Detection','FONTSIZE', FONTSIZE, 'Enable','off','callback', {@cb_but1_detection});
but1_5 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but1_pos_5,'string', '<<','FONTSIZE', FONTSIZE+1,'Enable','off', 'callback', {@cb_but56});
but1_6 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but1_pos_6,'string', '>>','FONTSIZE', FONTSIZE+1,'Enable','off', 'callback', {@cb_but56});
txtbig = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txtbig_pos,'string', '','BackgroundColor', txtbgcol,'ForegroundColor', txtfgcol, 'FONTSIZE', FONTSIZE);
txttitle = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txttitle_pos,'string', 'ROIFlow Event detection','FONTSIZE', round(FONTSIZE*1.3));
txtsm = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txtsm_pos,'string', '','FONTSIZE', FONTSIZE-2);

uiwait;


%% Local Callback functions
function open_trcfig()
global  fig2_pos trcfig sp2_1 sp2_2 imsp2_pos but2_po1_pos_1 but2_po1 trcsp2_pos...
    FONTSIZE LOADED SELROI DETECTED PLOTMODE
if ishandle(trcfig), figure(trcfig);
else, trcfig = figure('Name', 'Trace overview', 'Position', fig2_pos, 'toolbar', 'none', 'menu', 'none'); axis off;
end
    set(trcfig,'WindowKeyPressFcn',@keyPressCallback);
sp2_1 = subplot('Position', imsp2_pos); axis('off');
sp2_2 = subplot('Position', trcsp2_pos); axis('off');
but2_po1 = uicontrol('parent', trcfig, 'style', 'pushbutton','units', 'normalized','position', but2_po1_pos_1,'string', 'Discard ROI','FONTSIZE', FONTSIZE,'BackgroundColor', [.4 .4 .4], 'callback', {@cb_but2_del_roi});

if LOADED
    plot_rois(SELROI, PLOTMODE,sp2_1);
    plot_traces(SELROI, 'dFoF', DETECTED(SELROI));
end
end

function cb_but1_load(~,~)
global LOADED FILELIST SELFILE NFILES PATH FONTSIZE...
    but1_2 but1_3 but1_4 but1_5 but1_6 shownonev
if LOADED
    filesfig = figure('toolbar', 'none', 'menu', 'none'); axis('off');
    uicontrol('parent', filesfig, 'style', 'listbox','units', 'normalized','position', [.05 .05 .4 .8],'string', FILELIST([FILELIST{:,3}]==0,1),'FONTSIZE', FONTSIZE, 'callback', {@cb_selfile});
    uicontrol('parent', filesfig, 'style', 'text','units', 'normalized','position', [.05 .85 .4 .1],'string', 'NOT REVISED','FONTSIZE', FONTSIZE+2);
    uicontrol('parent', filesfig, 'style', 'listbox','units', 'normalized','position', [.55 .05 .4 .8],'string', FILELIST([FILELIST{:,3}]==1,1),'FONTSIZE', FONTSIZE, 'callback', {@cb_selfile});
    uicontrol('parent', filesfig, 'style', 'text','units', 'normalized','position', [.55 .85 .4 .1],'string', 'ACCEPTED','FONTSIZE', FONTSIZE+2);
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
    but1_2.Enable='on';
    but1_3.Enable='on';
    but1_4.Enable='on';
    but1_5.Enable='on';
    but1_6.Enable='on';
    load_mat();
end
end

function load_mat()
global PATH FILELIST SELFILE SELROI ROILIST...
    NROIS DETECTED FT CMAP roidata txtsm ddlist sp1 ALLROIS PLOTMODE
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
txtsm.String = FILELIST{SELFILE,1};
ddlist.String = ROILIST(:,1);
ddlist.Value = 1;
pause(0.1);
CMAP = get_colormap([1 0 0],[1 1 0],[0 1 1],NROIS);
plot_rois(1:NROIS, PLOTMODE, sp1);
open_trcfig();
end

function cb_but1_accept(~,~)
global SELFILE FILELIST EVDATA roidata...
    PATH FILELIST trcfig
if ishandle(trcfig), figure(trcfig); close gcf; end
save_event_info(EVDATA, roidata, PATH, FILELIST{SELFILE})
FILELIST{SELFILE,3} = true;

% Work on next file
tmpfl = find([FILELIST{:,3}]== 0);
if ~isempty(tmpfl)
    SELFILE = FILELIST{tmpfl(1),2};
    load_mat();
else, uiresume; close all;
end
end

function cb_but1_revision(~,~)
global fig3_pos txtbig sp2_2 ddlist EVDATA DETECTED ROILIST...
    roidata SELROI FT TRCMODE PLOTMODE
if ~isempty(EVDATA{SELROI,6})
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
        if DETECTED(SELROI)
            txtbig.String = {'',sprintf('Total events: %i', numel(EVDATA{SELROI,3})), '',...
                sprintf('Events for revision: %i', numel(EVDATA{SELROI,6})),'',...
                sprintf('Detect: > %ixSD | Save: > %ixSD',EVDATA{SELROI,10}(1), EVDATA{SELROI,10}(2))};
            pause(0.1);
        end
        clear('tmp_rev');
    end
    % Discard Roi from list if no events left
    if isempty(EVDATA{SELROI,7})
        currpter = find(([ROILIST{:,2}]==SELROI)==1);
        ROILIST(currpter,:) = [];
        if ~isempty(ROILIST), ddlist.Value = currpter; SELROI = ROILIST{currpter,2};
        else, ddlist.String = 'NO ROIs';ddlist.Enable='off'; SELROI = 0;
        end
        plot_rois([ROILIST{2,:}], PLOTMODE, sp1);
    end
end
end

function cb_but1_detection(~,~)
global roidata EVDATA ALLROIS PARAMS NROIS SELROI...
    sp2_2 sp1 txtbig ddlist DETECTED ROILIST TRCMODE PLOTMODE
% Detects events in all traces
% DFOFTHRESH = inputdlg('dF/F threshold for event detection:', 'Event detection', [60 10], DFOFTHRESH);
% DFOFTHRESH=str2num(DFOFTHRESH{1});
[EVDATA, detroi, PARAMS] = run_event_detection(EVDATA, ROILIST, roidata.dFoF_traces, roidata.FoF_traces, roidata.traces, PARAMS);
DETECTED(detroi) = true;
% Discard Rois from list if they dont show events
if any(cellfun(@isempty, EVDATA(:,7)))
    pterdel = cellfun(@isempty, EVDATA(:,7));
    pternew = find(pterdel==0 & [ALLROIS{:,3}]');
    delrois = find(pterdel==1);
    ROILIST = ALLROIS(pternew,:);
    if ~isempty(ddlist.Value) && any(ddlist.Value == delrois)
        ddlist.Value = 1;
        SELROI = ROILIST{1,2};
    end
else
    ROILIST = ALLROIS([ALLROIS{:,3}]');
    pternew=1:NROIS;
end
if isempty(pternew), ddlist.String = 'NO ROIs';ddlist.Enable='off'; SELROI = 0;
elseif pternew~=0, ddlist.Enable='on';ddlist.String = ROILIST(:,1);
end
plot_rois(pternew, PLOTMODE, sp1);
if ishandle(sp2_2), plot_traces(SELROI, TRCMODE, DETECTED(SELROI));
else, open_trcfig();
end
if DETECTED(SELROI)
    txtbig.String = {'',sprintf('Total events: %i', numel(EVDATA{SELROI,3})), '',...
        sprintf('Events for revision: %i', numel(EVDATA{SELROI,6})),'',...
        sprintf('Detect: > %ixSD | Save: > %ixSD',EVDATA{SELROI,10}(1), EVDATA{SELROI,10}(2))};
    pause(0.1);
end
end

function cb_ddlist(hOb,~)
global SELROI sp2_1 sp2_2 DETECTED EVDATA txtbig ROILIST TRCMODE PLOTMODE
SELROI = ROILIST{hOb.Value,2};
if ishandle(sp2_2)
    plot_traces(SELROI, TRCMODE, DETECTED(SELROI));
    plot_rois(SELROI, PLOTMODE,sp2_1);
else, open_trcfig();
end
if DETECTED(SELROI)
    txtbig.String = {'',sprintf('Total events: %i', numel(EVDATA{SELROI,3})), '',...
        sprintf('Events for revision: %i', numel(EVDATA{SELROI,6})),'',...
        sprintf('Detect: > %ixSD | Save: > %ixSD',EVDATA{SELROI,10}(1), EVDATA{SELROI,10}(2))};
    pause(0.1);
end
end

function cb_but56(hOb,~)
if strcmp(hOb.String,'<<'), swdir = -1;
else, swdir = 1; end
switch_trace(swdir);
end

function keyPressCallback(~,ev)
if strcmp(ev.Key,'rightarrow'), swdir = 1;switch_trace(swdir);
elseif strcmp(ev.Key,'leftarrow'), swdir = -1;switch_trace(swdir);
end
end

function switch_trace(swdir)
global SELROI sp2_1 sp2_2 DETECTED EVDATA txtbig ROILIST...
    TRCMODE ddlist PLOTMODE
tmpid = find(SELROI == [ROILIST{:,2}]);
if swdir == -1
    if tmpid > 1, SELROI = ROILIST{tmpid-1,2}; tmpid=tmpid-1; end
elseif swdir == 1
    if tmpid < size(ROILIST,1), SELROI = ROILIST{tmpid+1,2}; tmpid=tmpid+1; end
end
ddlist.Value = tmpid;
if ishandle(sp2_2)
    plot_traces(SELROI, TRCMODE, DETECTED(SELROI));
    plot_rois(SELROI, PLOTMODE,sp2_1);
else, open_trcfig();
end
if DETECTED(SELROI)
    txtbig.String = {'',sprintf('Total events: %i', numel(EVDATA{SELROI,3})), '',...
        sprintf('Events for revision: %i', numel(EVDATA{SELROI,6})),'',...
        sprintf('Detect: > %ixSD | Save: > %ixSD',EVDATA{SELROI,10}(1), EVDATA{SELROI,10}(2))};
    pause(0.1);
end
end

function cb_but2_del_roi(~,~)
global ALLROIS SELROI ROILIST EVDATA ddlist sp1 txtbig DETECTED TRCMODE PLOTMODE
% Discard Roi
ALLROIS{SELROI,3} = false;
currpter = find(([ROILIST{:,2}]==SELROI)==1);
ROILIST(currpter,:) = [];
if ~isempty(ROILIST)
    SELROI = ROILIST{currpter,2}; ddlist.String = ROILIST(:,1); ddlist.Value = currpter;
    plot_traces(SELROI, TRCMODE, DETECTED(SELROI));
    if DETECTED(SELROI)
        txtbig.String = {'',sprintf('Total events: %i', numel(EVDATA{SELROI,3})), '',...
            sprintf('Events for revision: %i', numel(EVDATA{SELROI,6})),'',...
            sprintf('Detect: > %ixSD | Save: > %ixSD',EVDATA{SELROI,10}(1), EVDATA{SELROI,10}(2))};
        pause(0.1);
    end
else,SELROI=0; ddlist.String = 'NO ROIs';ddlist.Enable='off';
end
plot_rois([ROILIST{2,:}], PLOTMODE, sp1);
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

function plot_rois(plotids, mode, sp)
global roidata CMAP
axes(sp);
if strcmp(mode,'aip')
    uint8_im = (roidata.aip./max(roidata.aip,[],'all'))*255;
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
        hold on, plot(boundary(:,2), boundary(:,1), 'Color',CMAP(:,roidata.roi_bounds{plotids(cc),2})', 'LineWidth', 1)
    end
end
end

function plot_traces(roiid, type, events)
global roidata FT FONTSIZE EVDATA SELROI sp2_2
axes(sp2_2);
trccol = [0 0 0]; %[0.6510 0.0667 0.2784];
if strcmp(type,'dFoF'), tmptraces = roidata.dFoF_traces; yaxlab = 'dFoF';
elseif strcmp(type,'FoF'), tmptraces = roidata.FoF_traces; yaxlab = 'FoF';
else, tmptraces = roidata.traces; yaxlab = 'Intensity [a.u.]';
end
plot(FT.*(1:size(tmptraces,2)), tmptraces(roiid,:), 'Color', trccol, 'LineWidth', 1);
yrange = [min(tmptraces(roiid,:)) max(tmptraces(roiid,:))];
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

function plot_raster()
end