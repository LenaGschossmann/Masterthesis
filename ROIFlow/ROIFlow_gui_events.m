%% User Interfacce for ROIFlow Analysis for 2PM and Confocal data of spontaneous events %%

addpath('C:\Users\lena_\Projects\code_extern\Matlab\WriteImageJRoi\');

global ROILIST REVINFO LOADED PATH FILELIST SELFILE NFILES FILTKERNEL...
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
FILTKERNEL = []; %[1 2 1]; [1 3 3 1]

shownonev = true;
PARAMS = get_predefined_params();
TRCMODE = 'FoF';
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
whtxtbg = [3 2.3].*whunits;
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
rb_pos_1 = [but1_pos_2(1)+but1_pos_2(3)+hvspace1(1)*0.2 but1_pos_2(2) whbut1(1)/4 whbut1(2)];
rb_pos_2 = [rb_pos_1(1) but1_pos_4(2) whbut1(1)/2 whbut1(2)];
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
rb_1 = uicontrol('parent', mainfig, 'style', 'radiobutton','units', 'normalized','position', rb_pos_1,'string', 'All','FONTSIZE', FONTSIZE-1,'Value', BATCH, 'callback', {@cb_rb_batch});
rb_2 = uicontrol('parent', mainfig, 'style', 'radiobutton','units', 'normalized','position', rb_pos_2,'string', 'Auto','FONTSIZE', FONTSIZE-1,'Value', AUTO, 'callback', {@cb_rb_auto});

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
global PATH FILELIST SELFILE SELROI ROILIST EVDATA FILTKERNEL PARAMS...
    NROIS DETECTED FT CMAP roidata txtsm ddlist sp1a sp1b ALLROIS PLOTMODE
load(fullfile(PATH,FILELIST{SELFILE,1}));
NROIS = roidata.n_rois;
SELROI = 1;
FT = roidata.frametime_s;
ROILIST = cell(NROIS,2);
ALLROIS = cell(NROIS,3);
DETECTED = false(NROIS,1);
EVDATA = struct();
% Filter traces
if ~isempty(FILTKERNEL), trc_filtrd = smooth_data(roidata.traces, FILTKERNEL);
else, trc_filtrd = roidata.traces;
end
% Calculate FoF
[~,FoF_filtrd, dFoF_filtrd] = rollBase_dFoF(trc_filtrd,roidata.baseline_frames,PARAMS.dFoF_baseline_shift, 'roll');
for iRoi = 1:NROIS
    ROILIST{iRoi,1} = sprintf('ROI %i',iRoi); ROILIST{iRoi,2} = iRoi;
    ALLROIS{iRoi,3} = true;
    EVDATA(iRoi).size_FiltKernel = numel(FILTKERNEL);
    EVDATA(iRoi).filtrd_trace = trc_filtrd(iRoi,:);
    EVDATA(iRoi).filtrd_FoFtrace = FoF_filtrd(iRoi,:);
    EVDATA(iRoi).filtrd_dFoFtrace = dFoF_filtrd(iRoi,:);
    EVDATA(iRoi).FoF_SDnoise = [];
end
ALLROIS(:,1:2) = ROILIST;
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
global SELFILE FILELIST EVDATA ALLROIS roidata...
    PATH trcfig BATCH PARAMS
if ishandle(trcfig), figure(trcfig); close gcf; end
tmpfl = find([FILELIST{:,3}]== 0);
if BATCH && size(FILELIST,1) > 1
    % Process all files with chosen setting
    process_batch(tmpfl);
else
    roiselection = [ALLROIS{:,3}];
    [~,~,~] = save_event_info(EVDATA, roiselection, roidata, PARAMS, PATH, FILELIST{SELFILE},false);
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

function cb_rb_auto(hOb,~)
global AUTO
AUTO = hOb.Value;
end

function cb_but1_revision(~,~)
global fig3_pos sp2_2 sp1a sp1b ddlist EVDATA DETECTED ROILIST...
    SELROI FT TRCMODE PLOTMODE PARAMS AUTO
goOn = false;
while ~goOn
    [goOn, keep_ev] = revise_events(AUTO, fig3_pos, PARAMS, EVDATA, SELROI, FT);
end
if ~AUTO && ~isempty(keep_ev)
    EVDATA(SELROI).onset_idx(~keep_ev) = [];
    EVDATA(SELROI).peak_idx(~keep_ev) = [];
%     EVDATA(SELROI).peak_amp(~keep_ev) = [];
%     EVDATA(SELROI).FoF_peak_amp(~keep_ev) = [];
    if ishandle(sp2_2), plot_traces(SELROI, TRCMODE, DETECTED(SELROI));
    else, open_trcfig();
    end
    if DETECTED(SELROI), update_main_fig();end
    clear('tmp_rev');
end
% Discard Roi from list if no events left
if ~AUTO && isempty(keep_ev) && all(~keep_ev)
    currpter = find(([ROILIST{:,2}]==SELROI)==1);
    ROILIST(currpter,:) = [];
    if ~isempty(ROILIST), ddlist.Value = currpter; SELROI = ROILIST{currpter,2};
    else, ddlist.String = 'NO ROIs';ddlist.Enable='off'; SELROI = 0;
    end
    plot_rois([ROILIST{:,2}], PLOTMODE, sp1a, sp1b);
end
end

function cb_but1_detection(~,~)
global EVDATA ALLROIS PARAMS NROIS SELROI...
    sp2_2 sp1a sp1b but1_3 but1_2 trcfig ddlist DETECTED ROILIST TRCMODE PLOTMODE FT AUTO
% Detects events in all traces
if SELROI == 0, ROILIST = ALLROIS(:,1:2); SELROI = 1; end
[EVDATA, detroi, PARAMS] = run_event_detection(EVDATA, ROILIST, PARAMS, FT, AUTO, true);
DETECTED(detroi) = true;
but1_3.Enable='on';
but1_2.Enable='on';
% Discard Rois from list if they dont show events

if any(~([EVDATA(:).events_detected]))
    pternew = find([EVDATA(:).events_detected] == 1 & [ALLROIS{:,3}]);
    delrois = find([EVDATA(:).events_detected] == 0);
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
global SELROI sp2_1 sp2_2 DETECTED ROILIST...
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
global ALLROIS SELROI ROILIST ddlist sp1a sp1b sp2_1 trcfig DETECTED TRCMODE PLOTMODE
% Discard Roi
ALLROIS{SELROI,3} = false;
currpter = find(([ROILIST{:,2}]==SELROI)==1);
ROILIST(currpter,:) = [];
if ~isempty(ROILIST)
    if currpter > size(ROILIST,1), currpter = 1; end
    SELROI = ROILIST{currpter,2}; ddlist.String = ROILIST(:,1); ddlist.Value = currpter;
    plot_traces(SELROI, TRCMODE, DETECTED(SELROI));
    if DETECTED(SELROI), update_main_fig(); end
    plot_rois([ROILIST{:,2}], PLOTMODE, sp1a,sp1b);
    plot_rois(SELROI, PLOTMODE, sp2_1,[]);
else,SELROI=0; ddlist.String = 'NO ROIs';ddlist.Enable='off';close(trcfig);
end
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
    boundary = roidata.roi_bounds{plotids(cc),1}; boundary = boundary{1};
    hold on
    plot(boundary(:,2), boundary(:,1), 'Color',CMAP(:,roidata.roi_bounds{plotids(cc),2})', 'LineWidth', 1)
end

if ~isempty(spb)
    axes(spb); cla;
    %     v1 = min(cellfun(@min, roidata.roi_bounds(:,2)));
    v2 = max(cellfun(@max, roidata.roi_bounds(:,2)));
    if v2 > 12 && v2 < 40
        stepsize = 5;
        steps = ceil(v2/stepsize);
        v2 = ceil(v2/steps)*steps;
    elseif v2 >=40
        stepsize = 10;
        steps = ceil(v2/stepsize);
        v2 = ceil(v2/steps)*steps;
    else
        steps = v2;
    end
    colormap(spb, CMAP(:,1:v2)');
    cb = colorbar();
    cbtl = linspace(1, v2, steps)';
    tmplim = get(cb,'YLim');
    cbt = linspace(tmplim(1), tmplim(2), steps);
    set(cb, 'YTick', cbt, 'YTickLabel',cbtl, 'Location', 'eastoutside');
end

end

function plot_traces(roiid, type, events)
global FT FONTSIZE EVDATA SELROI sp2_2
axes(sp2_2);
cla;
trccol = [0 0 0]; %[0.6510 0.0667 0.2784];
if strcmp(type,'dFoF'), tmptrace = EVDATA(roiid).filtrd_FoFtrace; yaxlab = 'dFoF'; %tmptraces = roidata.dFoF_traces; 
elseif strcmp(type,'FoF'), tmptrace = EVDATA(roiid).filtrd_FoFtrace; yaxlab = 'FoF'; %tmptraces = roidata.FoF_traces
else, tmptrace = EVDATA(roiid).filtrd_trace; yaxlab = 'Intensity [a.u.]'; % tmptraces = roidata.traces;
end
plot(FT.*(1:size(tmptrace,2)), tmptrace, 'Color', trccol, 'LineWidth', 1);
if ~isempty(EVDATA(SELROI).FoF_SDnoise) && strcmp(type,'FoF')
    hold on, yline(1-EVDATA(SELROI).FoF_SDnoise, 'Color',[.5 .5 .5]); yline(1+EVDATA(SELROI).FoF_SDnoise, 'Color',[.5 .5 .5]);
    
    hold on, yline(EVDATA(SELROI).FoF_threshold, 'cyan');
    text(1,EVDATA(SELROI).FoF_threshold, strcat('det.thresh. ',num2str(EVDATA(SELROI).FoF_threshold)), 'Fontsize', 6);
    
    if ~isempty(EVDATA(SELROI).ctrl_av_peak_FoF)
        hold on, yline(EVDATA(SELROI).ctrl_av_peak_FoF, 'blue');
        text(1,EVDATA(SELROI).ctrl_av_peak_FoF, 'ctrl av. peak amp ', 'Fontsize', 6);
    end
end
yrange = [min(tmptrace) max(tmptrace)];
ylim([yrange(1)-diff(yrange)/10 yrange(2)+diff(yrange)/10]);
xlim([FT size(tmptrace,2)*FT]);
xlabel('Time [s]'); ylabel(yaxlab);
set(gca, 'FontSize', FONTSIZE-1);
if events
    rev_pts = EVDATA(SELROI).peak_idx .*FT;
    hold on;
    for iEv = 1:numel(rev_pts), xline(rev_pts(iEv), 'r', 'LineWidth', 1.2, 'Alpha', .4); end
    hold off;
end
end

function update_main_fig()
global txtbig EVDATA SELROI NROIS
tot_ev = 0; all_peak_amp = [];
for iEv = 1:size(EVDATA,2)
    tot_ev = tot_ev + numel(EVDATA(SELROI).onset_idx);
    all_peak = [all_peak_amp EVDATA(SELROI).filtrd_FoFtrace(EVDATA(SELROI).peak_idx)'];
end
min_peak = round(min(all_peak),2);
av_peak = round(mean(all_peak),2);
max_peak = round(max(all_peak),2);

txtbig.String = {'',sprintf('ROIs with events: %i / %i', sum([EVDATA.events_detected]), NROIS), '',...
    sprintf('# Events curr. ROI: %i', numel(EVDATA(SELROI).onset_idx)),'',...
    sprintf('# Events total: %i', tot_ev), '',...
    strcat('Min/Av/Max Peak: ',num2str(min_peak),'/',num2str(av_peak),'/',num2str(max_peak))};
pause(0.1);
end

function process_batch(file_idx)
global PATH FILELIST PARAMS FILTKERNEL
nfiles = numel(file_idx);
rec_master = struct();
roi_master = struct();
synchronicity_master = struct();
sync_cnt = 0;
roi_cnt = 0;
sum_cnt = 1;
for iFile = 1:nfiles
    fprintf('Processing file %i / %i\n', iFile, nfiles);
    tmpfile = FILELIST{file_idx(iFile),2};   
    evdata = struct();
    
    % Load
    load(fullfile(PATH,FILELIST{tmpfile,1}));
    nrois = roidata.n_rois;
    roilist = cell(nrois,3);
    % Filter traces
    if ~isempty(FILTKERNEL), trc_filtrd = smooth_data(roidata.traces, FILTKERNEL);
    else, trc_filtrd = roidata.traces;
    end
    % Calculate FoF
    [~,FoF_filtrd, dFoF_filtrd] = rollBase_dFoF(trc_filtrd,roidata.baseline_frames,PARAMS.dFoF_baseline_shift, 'roll');
    for iRoi = 1:nrois
        roilist{iRoi,1} = sprintf('ROI %i',iRoi); roilist{iRoi,2} = iRoi; roilist{iRoi,3} = true;
        evdata(iRoi).size_FiltKernel = numel(FILTKERNEL);
        evdata(iRoi).filtrd_trace = trc_filtrd(iRoi,:);
        evdata(iRoi).filtrd_FoFtrace = FoF_filtrd(iRoi,:);
        evdata(iRoi).filtrd_dFoFtrace = dFoF_filtrd(iRoi,:);
        evdata(iRoi).FoF_SDnoise = [];
    end
    cmap_n = max(cellfun(@max ,roidata.roi_bounds(:,2)));
    cmap = get_colormap([1 0 0],[1 1 0],[0 1 1],cmap_n);
    % cmap = get_colormap([1 0 0],[1 1 0],[0 1 1],NROIS);
    
    % Detect
    [evdata, ~, PARAMS] = run_event_detection(evdata, roilist, PARAMS, roidata.frametime_s, true, false);
    
    % Discard Rois from list if they dont show events
    for iRoi = 1:size(evdata,1)
        if isempty(evdata(iRoi).onset_idx)
            roilist{iRoi,3} = false;
        end
    end
    
    % Save
    roiselection = [roilist{:,3}];
    [tmp_summary, tmp_roimaster, tmp_sync_dist] = save_event_info(evdata, roiselection, roidata, PARAMS, PATH, FILELIST{tmpfile}, true);
    %         FILELIST{tmpfile,3} = true;
    
    % Write all rercording summaries into master table
    if ~isempty(tmp_summary.Num_ROIs)
        rec_master(sum_cnt).Rec_Name = FILELIST{tmpfile}(1:end-4);
        rec_master(sum_cnt).Num_ROIs = tmp_summary.Num_ROIs;
        rec_master(sum_cnt).ToT_Num_Events = tmp_summary.ToT_Num_Events;
        rec_master(sum_cnt).Num_Events_Mean = tmp_summary.Num_Events_Mean;
        rec_master(sum_cnt).Amp_Mean = tmp_summary.Amp_Mean;
        rec_master(sum_cnt).Amp_FoF_Mean = tmp_summary.Amp_FoF_Mean;
        rec_master(sum_cnt).Peak_dFoF_Mean = tmp_summary.Peak_dFoF_Mean;
        rec_master(sum_cnt).Amp_SD = tmp_summary.Amp_SD;
        rec_master(sum_cnt).Amp_FoF_SD = tmp_summary.Amp_FoF_SD;
        rec_master(sum_cnt).Peak_dFoF_SD = tmp_summary.Peak_dFoF_SD;
        rec_master(sum_cnt).Mean_SNR = tmp_summary.Mean_SNR;
        rec_master(sum_cnt).IEI_Mean = tmp_summary.IEI_Mean;
        rec_master(sum_cnt).EvRate_Hz_Mean = tmp_summary.EvRate_Hz_Mean;
        rec_master(sum_cnt).REC_time_s = tmp_summary.REC_time_s;
        rec_master(sum_cnt).Area_um_Mean = tmp_summary.Area_um_Mean;
        rec_master(sum_cnt).Area_ROI_PercOf_Struct = tmp_summary.Area_ROI_PercOf_Struct;
        rec_master(sum_cnt).Dend_Area_um = tmp_summary.Dend_Area_um;
        rec_master(sum_cnt).Bg_Area_um = tmp_summary.Bg_Area_um;
        rec_master(sum_cnt).FoV_Area_um = tmp_summary.FoV_Area_um;
        rec_master(sum_cnt).Events_per_Dend_Area = tmp_summary.Events_per_Dend_Area;
        sum_cnt = sum_cnt+1;
        
        % ROI master table
        roi_n = size(tmp_roimaster,2);
        for iRoi = 1:roi_n
            roi_master(roi_cnt+iRoi).Recording = FILELIST{tmpfile,1}(1:end-4);
            roi_master(roi_cnt+iRoi).ROI = tmp_roimaster(iRoi).ROI;
            roi_master(roi_cnt+iRoi).Num_Events = tmp_roimaster(iRoi).Num_Events;
            roi_master(roi_cnt+iRoi).Amp_Mean = tmp_roimaster(iRoi).Amp_Mean;
            roi_master(roi_cnt+iRoi).Amp_FoF_Mean = tmp_roimaster(iRoi).Amp_FoF_Mean;
            roi_master(roi_cnt+iRoi).Peak_dFoF_Mean = tmp_roimaster(iRoi).Peak_dFoF_Mean;
            roi_master(roi_cnt+iRoi).Amp_SD = tmp_roimaster(iRoi).Amp_SD;
            roi_master(roi_cnt+iRoi).Amp_FoF_SD = tmp_roimaster(iRoi).Amp_FoF_SD;
            roi_master(roi_cnt+iRoi).Peak_dFoF_SD = tmp_roimaster(iRoi).Peak_dFoF_SD;
            roi_master(roi_cnt+iRoi).Trc_Mean = tmp_roimaster(iRoi).Trc_Mean;
            roi_master(roi_cnt+iRoi).Trc_FoF_Mean = tmp_roimaster(iRoi).Trc_FoF_Mean;
            roi_master(roi_cnt+iRoi).Trc_dFoF_Mean = tmp_roimaster(iRoi).Trc_dFoF_Mean;
            roi_master(roi_cnt+iRoi).Trc_SD = tmp_roimaster(iRoi).Trc_SD;
            roi_master(roi_cnt+iRoi).Trc_FoF_SD = tmp_roimaster(iRoi).Trc_FoF_SD;
            roi_master(roi_cnt+iRoi).Trc_dFoF_SD = tmp_roimaster(iRoi).Trc_dFoF_SD;
            roi_master(roi_cnt+iRoi).Threshold_Slope = tmp_roimaster(iRoi).Threshold_Slope;
            roi_master(roi_cnt+iRoi).Threshold_Amp = tmp_roimaster(iRoi).Threshold_Amp;
            roi_master(roi_cnt+iRoi).Trc_SNR = tmp_roimaster(iRoi).Trc_SNR;
            roi_master(roi_cnt+iRoi).IEI_Mean = tmp_roimaster(iRoi).IEI_Mean;
            roi_master(roi_cnt+iRoi).CV_IEI = tmp_roimaster(iRoi).CV_IEI;
            roi_master(roi_cnt+iRoi).EvRate_Hz = tmp_roimaster(iRoi).EvRate_Hz;
            roi_master(roi_cnt+iRoi).Area_um = tmp_roimaster(iRoi).Area_um;
            roi_master(roi_cnt+iRoi).Area_ROI_PercOf_Struct = tmp_roimaster(iRoi).Area_ROI_PercOf_Struct;
            roi_master(roi_cnt+iRoi).Ctr_X_um = tmp_roimaster(iRoi).Ctr_X_um;
            roi_master(roi_cnt+iRoi).Ctr_Y_um = tmp_roimaster.Ctr_Y_um;
            roi_master(roi_cnt+iRoi).Subtract_value = tmp_roimaster(iRoi).Subtract_value;
            roi_master(roi_cnt+iRoi).Offset = tmp_roimaster(iRoi).Offset;
            roi_master(roi_cnt+iRoi).Pmt_gain = tmp_roimaster(iRoi).Pmt_gain;
        end
        roi_cnt = roi_cnt + roi_n;
        
        % Synchronicity and Distance master table
        if tmp_summary.Num_ROIs > 1
            sync_n = size(tmp_sync_dist,1);
            for iSync = 1:sync_n
                synchronicity_master(sync_cnt+iSync).Recording = FILELIST{tmpfile,1}(1:end-4);
                synchronicity_master(sync_cnt+iSync).ROI_1 = tmp_sync_dist.ROI_1(iSync);
                synchronicity_master(sync_cnt+iSync).ROI_2 = tmp_sync_dist.ROI_2(iSync);
                synchronicity_master(sync_cnt+iSync).Dist_um = tmp_sync_dist.Dist_um(iSync);
                synchronicity_master(sync_cnt+iSync).PearsonR = tmp_sync_dist.PearsonR(iSync);
                synchronicity_master(sync_cnt+iSync).EvRate_Hz_ROI_1 = tmp_sync_dist.EvRate_Hz_ROI_1(iSync);
                synchronicity_master(sync_cnt+iSync).EvRate_Hz_ROI_2 = tmp_sync_dist.EvRate_Hz_ROI_2(iSync);
            end
            sync_cnt = sync_cnt + sync_n;
        end
    end
end

% Write master table to file
savepath = fullfile(PATH,FILELIST{tmpfile}(1:8));
rec_master_xls = strcat(savepath, '_Rec_master.xlsx'); if isa(rec_master_xls,'cell'), rec_master_xls=rec_master_xls{1};end
roi_master_xls = strcat(savepath, '_ROI_master.xlsx'); if isa(roi_master_xls,'cell'), roi_master_xls=roi_master_xls{1};end
sync_master_xls = strcat(savepath, '_Sync_master.xlsx'); if isa(sync_master_xls,'cell'), sync_master_xls=sync_master_xls{1};end
writetable(struct2table(rec_master), rec_master_xls);
writetable(struct2table(roi_master), roi_master_xls);
writetable(struct2table(synchronicity_master), sync_master_xls);

msgbox('All files processed!');
uiresume; close all;
end