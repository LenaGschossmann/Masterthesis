%% User Interfacce for ROIFlow Analysis for 2PM and Confocal data of spontaneous events %%

clear all;
addpath(genpath('C:\Users\lena_\Projects\code\Masterthesis\OpenBFiles'));
addpath('C:\Users\lena_\Projects\code\Masterthesis\ROIFlow\test_scripts');

%% Parameters
global FONTSIZE SELECTED FILELIST NFILES INFOSTRUCT FTINFO SCALEINFO BGINFO AIPLIST...
    but_2 but_4 but_5  list_1 txt_4 txtlist_3 txtlist_4  txtlist_5 txtlist_6...
    PARAMS setcol unsetcol

FONTSIZE = 10;
SELECTED = 0;
FILELIST = {''};
NFILES = 0;
INFOSTRUCT = struct('name',[],'path',[],'frametime_s',[], 'bg_threshold',[],'mocorr',[]);
FTINFO = struct('pointer', [],'ft_s',[], 'txt',[]);
SCALEINFO = struct('pointer', [],'xy_scale',[], 'txt',[]);
BGINFO = struct('pointer', [],'bg_val',[], 'txt',[]);
AIPLIST = [];
PARAMS = get_predefined_params('subtract_perc', 0.15,...
    'connc_px_thresh',10,'perc_thresh',95, 'perc_fac', 1.15, 'overlap_thresh',0.5);

%% GUI Parameters
scrsz = get(0,'ScreenSize'); scrsz=scrsz(3:4);
whfig = [1000 600];
whunits = 1./[18 20];
% Units normalized
whlist1 = [5 18].*whunits;
whlist2 = [3 14].*whunits;
whlist3 = [1 14].*whunits;
whbut = [2 1.5].*whunits;
whbutsm = [2 0.8].*whunits;
whtxt = [3 0.7].*whunits;
hvspacer = [0.7 1].*whunits;

txtbgcol = [0.2 0.2 0.2];
txtfgcol = [1 1 1];
setcol = [.2 .5 .8];
unsetcol = [.8 .8 .8];

% GUI positions
fig_pos = [round((scrsz(1)-whfig(1))/2) round((scrsz(2)-whfig(2))/2) whfig];
list1_pos = [hvspacer whlist1];
but_pos_0 = [list1_pos(1)+whlist1(1)+hvspacer(1) hvspacer(2)*2.3 whbut(1) whbut(2)*0.75];
but_pos_1 = [list1_pos(1)+whlist1(1)+hvspacer(1) hvspacer(2)*2 whbut(1) whbut(2)*2];
but_pos_2 = [list1_pos(1)+whlist1(1)+hvspacer(1) but_pos_1(2)+but_pos_1(4)+hvspacer(2)*0.7 whbut(1) whbut(2)*0.8];
but_pos_3 = [list1_pos(1)+whlist1(1)+hvspacer(1) but_pos_2(2)+but_pos_2(4)+hvspacer(2)*0.2 whbut(1) whbut(2)*0.8];
but_pos_4 = [list1_pos(1)+whlist1(1)+hvspacer(1) but_pos_3(2)+but_pos_3(4)+hvspacer(2)*0.5 whbut];
but_pos_5 = [list1_pos(1)+whlist1(1)+hvspacer(1) but_pos_4(2)+whbut(2)+hvspacer(2)*0.5 whbut];
but_pos_7 = [list1_pos(1)+whlist1(1)+hvspacer(1) but_pos_5(2)+whbut(2)+hvspacer(2)*0.7 whbut(1) whbut(2)*2];
txtlist1_pos_2 = [but_pos_1(1)+whbut(1)+hvspacer(1) hvspacer(2)*1.2+whbutsm(2) whlist2];
txtlist2_pos_2 = [txtlist1_pos_2(1)+whlist2(1)+hvspacer(1)/20 txtlist1_pos_2(2) whlist3];
txtlist1_pos_3 = [txtlist2_pos_2(1)+whlist3(1)+hvspacer(1)*0.75 hvspacer(2)*1.2+whbutsm(2) whlist2];
txtlist2_pos_3 = [txtlist1_pos_3(1)+whlist2(1)+hvspacer(1)/20 txtlist1_pos_2(2) whlist3];
txt_pos_2 = [txtlist1_pos_2(1) txtlist1_pos_2(2)+whlist2(2)+0.1*hvspacer(2) whtxt];
txt_pos_3 = [txtlist1_pos_3(1) txtlist1_pos_2(2)+whlist2(2)+0.1*hvspacer(2) whtxt];
txt_pos_4 = [list1_pos(1)+whlist1(1)+hvspacer(1)*3 list1_pos(2)+whlist1(2)-whtxt(2)*2 whtxt(1)*1.5 whtxt(2)];
butsm_pos_2 = [txtlist1_pos_2(1)+hvspacer(1)*0.3 hvspacer(2) whbutsm];
butsm_pos_3 = [txtlist1_pos_3(1)+hvspacer(1)*0.3 hvspacer(2) whbutsm];

mainfig = figure('Name', 'ROIFlow - Automated Timeseries Analysis', 'Position', fig_pos, 'toolbar', 'none', 'menu', 'none');
set(gca, 'Color', 'none'); axis off;
title('ROIFlow Batch Processing', 'FontSize', round(FONTSIZE*1.3));

% Add UI controls
list_1 = uicontrol('parent', mainfig, 'style', 'listbox','units', 'normalized','position', list1_pos,'FONTSIZE', FONTSIZE-1, 'string',FILELIST,'Callback', {@cb_list1});
txt_2 = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txt_pos_2,'string', 'Frame time set','FONTSIZE', FONTSIZE);
txt_3 = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txt_pos_3,'string', 'XY Scale set','FONTSIZE', FONTSIZE);
txt_4 = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txt_pos_4,'string', '','FONTSIZE', FONTSIZE);
txtlist_3 = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txtlist1_pos_2,'FONTSIZE', FONTSIZE-2, 'string','','BackgroundColor', txtbgcol, 'ForegroundColor', txtfgcol);
txtlist_4 = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txtlist2_pos_2,'FONTSIZE', FONTSIZE-2, 'string','','BackgroundColor', txtbgcol, 'ForegroundColor', txtfgcol);
txtlist_5 = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txtlist1_pos_3,'FONTSIZE', FONTSIZE-2, 'string','','BackgroundColor', txtbgcol, 'ForegroundColor', txtfgcol);
txtlist_6 = uicontrol('parent', mainfig, 'style', 'text','units', 'normalized','position', txtlist2_pos_3,'FONTSIZE', FONTSIZE-2, 'string','','BackgroundColor', txtbgcol, 'ForegroundColor', txtfgcol);
but_1 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but_pos_1,'string', 'Start','FONTSIZE', FONTSIZE+4, 'BackgroundColor', [.6 .6 .6],'callback', {@cb_but1_start});
but_2 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but_pos_2,'string', 'MoCorr','FONTSIZE', FONTSIZE,'BackgroundColor',unsetcol, 'callback', {@cb_but2_mocorr});
but_3 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but_pos_3,'string', '2PM DarkCount','FONTSIZE', FONTSIZE,'BackgroundColor',unsetcol, 'callback', {@cb_but3_2pm_corr});
but_4 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but_pos_4,'string', 'Set FrameTime','FONTSIZE', FONTSIZE,'BackgroundColor',unsetcol, 'callback', {@cb_but4_frametime});
but_5 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but_pos_5,'string', 'Set Scale','FONTSIZE', FONTSIZE, 'BackgroundColor',unsetcol,'callback', {@cb_but5_setscale});
but_7 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', but_pos_7,'string', 'Open','FONTSIZE', FONTSIZE+2, 'BackgroundColor', [.6 .6 .6],'callback', {@cb_but7_open});
butsm_2 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', butsm_pos_2,'string', 'Check FT','FONTSIZE', FONTSIZE, 'callback', {@cb_butsm});
butsm_3 = uicontrol('parent', mainfig, 'style', 'pushbutton','units', 'normalized','position', butsm_pos_3,'string', 'Check Scale','FONTSIZE', FONTSIZE, 'callback', {@cb_butsm});
uiwait;

%% Local Callback functions
function cb_list1(hObj,~)
global SELECTED INFOSTRUCT but_2 but_4 but_5 setcol unsetcol
SELECTED = hObj.Value;
if numel(SELECTED) == 1
    if ~isnan(INFOSTRUCT(SELECTED).frametime_s), but_4.BackgroundColor = setcol; else, but_4.BackgroundColor = unsetcol; end
    if ~isnan(INFOSTRUCT(SELECTED).xy_scale), but_5.BackgroundColor = setcol; else, but_5.BackgroundColor = unsetcol; end
    if INFOSTRUCT(SELECTED).mocorr, but_2.BackgroundColor = setcol; else, but_2.BackgroundColor = unsetcol; end
else
    but_2.BackgroundColor = unsetcol; but_4.BackgroundColor = unsetcol;
    but_5.BackgroundColor = unsetcol;
end
end

function cb_but7_open(~,~)
global NFILES FILELIST INFOSTRUCT AIPLIST list_1 txt_4
[filenames, path] = uigetfile({'*.tif*'; '*.lsm'; '*.nd2'}, 'Select files for processing','Multiselect', 'on');
if isa(filenames,'char'), filenames = {filenames}; end
NFILES = numel(filenames);
FILELIST = cell(NFILES,1);
AIPLIST = cell(NFILES,1);
t_load = 30;
for iF = 1:NFILES
    tic;
    txt_4.String = sprintf('Loading %i/%i | Expected time: %i min',iF, NFILES, round(t_load*NFILES/60));
    pause(0.1);
    AIPLIST{iF} = read_timeseries(fullfile(path, filenames{iF}), 'aip');
    INFOSTRUCT(iF).path = fullfile(path, filenames{iF});
    FILELIST{iF,1} = filenames{iF}(1:end-4);
    INFOSTRUCT(iF).name = FILELIST{iF};
    INFOSTRUCT(iF).mocorr = false;
    INFOSTRUCT(iF).corr_2pm = '\';
    INFOSTRUCT(iF).pmt_gain = NaN;
    INFOSTRUCT(iF).offset = NaN;
    INFOSTRUCT(iF).refset = NaN;
    INFOSTRUCT(iF).frametime_s = NaN;
    INFOSTRUCT(iF).xy_scale = NaN;
    list_1.String = FILELIST;
    list_1.Max = NFILES;
    t_load = toc;
end
txt_4.String = 'Finished Loading';
pause(0.1);
end

function cb_but2_mocorr(~,~)
global INFOSTRUCT SELECTED FONTSIZE but_2 setcol unsetcol
if ~isempty(SELECTED) && SELECTED(1)  ~= 0
    for iSel = 1:numel(SELECTED)
        if ~INFOSTRUCT(SELECTED(iSel)).mocorr
            INFOSTRUCT(SELECTED(iSel)).mocorr = true;
            but_2.BackgroundColor = setcol;
        else
            INFOSTRUCT(SELECTED(iSel)).mocorr = false;
            but_2.BackgroundColor = unsetcol;
        end
    end
    tmppter = find([INFOSTRUCT(:).mocorr] == true);
    mocorrlist = cell(numel(tmppter),1);
    for iF = 1:numel(tmppter), mocorrlist{iF,1} = INFOSTRUCT(tmppter(iF)).name; end
    mocorrfig = figure('toolbar', 'none', 'menu', 'none');
    title('Motion Correction');
    uicontrol('parent', mocorrfig, 'style', 'text','units', 'normalized','position', [.05 .05 .9 .8],'string', mocorrlist,'FONTSIZE', FONTSIZE);
end
end

function cb_but3_2pm_corr(~,~)
global INFOSTRUCT SELECTED but_3 txt_4 setcol
if ~isempty(SELECTED) && SELECTED(1)  ~= 0
    % Read excel file with the appropriate value to correct for the offset
    % during the 2P recording. This value is the difference of a background
    % recording with offset 0 and with the offset used for the recording
    % which is analyzed
    [f_corrfile, p_corrfile] = uigetfile({'*.csv*'; '*.xlsx'; '*.xls'}, 'Select file for background correction','Multiselect', 'off');
    f_corrfile = fullfile(p_corrfile,f_corrfile);
    if ~isa(f_corrfile,'char'), f_corrfile = f_corrfile{1}; end
    for iSel = 1:numel(SELECTED)
        tmpf = readmatrix(f_corrfile);
        INFOSTRUCT(SELECTED(iSel)).pmt_gain = tmpf(1);
        INFOSTRUCT(SELECTED(iSel)).offset = tmpf(2);
        INFOSTRUCT(SELECTED(iSel)).corr_2pm = tmpf(3);
    end
    clear('tmpf');
    but_3.BackgroundColor = setcol;
    txt_4.String = ''; pause(0.1);
end
end

function cb_but4_frametime(~,~)
global INFOSTRUCT SELECTED FTINFO FILELIST but_4 txt_4 txtlist_3 txtlist_4 setcol
if ~isempty(SELECTED) && SELECTED(1) ~= 0
    txt_4.String = 'Set frame time...';
    pause(0.1);
    tmp_ft = inputdlg('Enter frame time of selected files in ms:');
    if ~isempty(tmp_ft)
        ft_s = str2num(tmp_ft{1})/1000;
        for iSel = 1:numel(SELECTED)
            INFOSTRUCT(SELECTED(iSel)).frametime_s = ft_s;
            idx = SELECTED(iSel) == [FTINFO.pointer];
            if all(~idx)
                FTINFO.pointer = [FTINFO.pointer SELECTED(iSel)];
                FTINFO.ft_s = [FTINFO.ft_s ft_s];
                nchar = numel(FILELIST{SELECTED(iSel)});
                FTINFO.txt = [FTINFO.txt; strcat(FILELIST{SELECTED(iSel)}(1:10),'..',FILELIST{SELECTED(iSel)}(nchar-5:nchar))];
            else
                FTINFO.ft_s(idx) = ft_s;
            end
        end
        txtlist_3.String = FTINFO.txt;
        txtlist_4.String = FTINFO.ft_s;
        pause(0.1);
        but_4.BackgroundColor = setcol;
    end
else
    txt_4.String = 'SELECT FILE!';
    pause(0.1);
end
txt_4.String = ''; pause(0.1);
end

function cb_but5_setscale(~,~)
global INFOSTRUCT SELECTED SCALEINFO FILELIST txt_4 txtlist_5 txtlist_6 but_5 setcol
if ~isempty(SELECTED) && SELECTED(1) ~= 0
    txt_4.String = 'Set XY Scale...';
    pause(0.1);
    tmp_scale = inputdlg('Enter XY scaling in µm:');
    if ~isempty(tmp_scale)
        tmp_scale = str2num(tmp_scale{1});
        for iSel = 1:numel(SELECTED)
            INFOSTRUCT(SELECTED(iSel)).xy_scale = tmp_scale;
            idx = SELECTED(iSel) == [SCALEINFO.pointer];
            if all(~idx)
                SCALEINFO.pointer = [SCALEINFO.pointer SELECTED(iSel)];
                SCALEINFO.xy_scale = [SCALEINFO.xy_scale tmp_scale];
                nchar = numel(FILELIST{SELECTED(iSel)});
                SCALEINFO.txt = [SCALEINFO.txt; strcat(FILELIST{SELECTED(iSel)}(1:10),'..',FILELIST{SELECTED(iSel)}(nchar-5:nchar))];
            else
                SCALEINFO.xy_scale(idx) = tmp_scale;
            end
        end
        txtlist_5.String = SCALEINFO.txt;
        txtlist_6.String = SCALEINFO.xy_scale;
        pause(0.1);
        but_5.BackgroundColor = setcol;
    end
else
    txt_4.String = 'SELECT FILE!';
    pause(0.1);
end
txt_4.String = ''; pause(0.1);
end

% function cb_but6_setbg(~,~)
% global INFOSTRUCT SELECTED BGINFO FILELIST AIPLIST txt_4...
%     PARAMS setcol
% if ~isempty(SELECTED) && SELECTED(1) ~= 0
%     txt_4.String = 'Set Frametime...';
%     pause(0.1);
%     [val_mtrx, ~]= define_bg_threshold(AIPLIST,SELECTED, PARAMS.bg_perc);
%     for iSel = 1:numel(SELECTED)
%         INFOSTRUCT(SELECTED(iSel)).bg_threshold = val_mtrx(iSel,1);
%         idx = SELECTED(iSel) == [BGINFO.pointer];
%         if all(~idx)
%             BGINFO.pointer = [BGINFO.pointer SELECTED(iSel)];
%             BGINFO.bg_val = [BGINFO.bg_val val_mtrx(iSel,1)];
%             nchar = numel(FILELIST{SELECTED(iSel)});
%             BGINFO.txt = [BGINFO.txt; strcat(FILELIST{SELECTED(iSel)}(1:10),'..',FILELIST{SELECTED(iSel)}(nchar-5:nchar))];
%         else
%             BGINFO.bg_val(idx) = val_mtrx(iSel,1);
%         end
%     end
%     pause(0.1);
% else
%     txt_4.String = 'SELECT FILE!';
%     pause(0.1);
% end
% txt_4.String = ''; pause(0.1);
% end


function cb_butsm(hOb,~)
global INFOSTRUCT FONTSIZE SELECTED
if SELECTED(1) == 0
    msgbox('Please load files first!');
else
    checkfig = figure('toolbar', 'none', 'menu', 'none'); axis('off');
    if strcmp(hOb.String,'Check FT')
        title('Frametime - MISSING FILES');
        tmppter = find(isnan([INFOSTRUCT(:).frametime_s])==1);
        checklist = cell(numel(tmppter),1);
        for iF = 1:numel(tmppter), checklist{iF} = INFOSTRUCT(tmppter(iF)).name; end
    elseif strcmp(hOb.String,'Check Scale')
        title('XY Scale - MISSING FILES');
        tmppter = find(isnan([INFOSTRUCT(:).xy_scale])==1);
        checklist = cell(numel(tmppter),1);
        for iF = 1:numel(tmppter), checklist{iF} = INFOSTRUCT(tmppter(iF)).name; end
    else
        % Check files with associated file for correction
        bgcorrlist = [];
        for iF = 1:size(INFOSTRUCT,2)
            if strcmp(INFOSTRUCT(iF).corr_2pm, '\')
                bgcorrlist = [bgcorrlist; INFOSTRUCT(iF).name];
            end
        end        
        if ~isempty(bgcorrlist)
            bgcorrfig = figure('toolbar', 'none', 'menu', 'none');
            title('Missinng Bg Correction File');
            axis('off')
            uicontrol('parent', bgcorrfig, 'style', 'text','units', 'normalized','position', [.05 .05 .9 .8],'string', bgcorrlist,'FONTSIZE', FONTSIZE);
        end
    end
    uicontrol('parent', checkfig, 'style', 'text','units', 'normalized','position', [.05 .05 .9 .8],'string', checklist,'FONTSIZE', FONTSIZE);
end
end

function cb_but1_start(~,~)
global INFOSTRUCT SELECTED
if SELECTED == 0
    msgbox('Please load files first!');
elseif any(isnan([INFOSTRUCT(:).frametime_s]))
    msgbox('There are files without assigned frame time!');
% elseif any(isnan([INFOSTRUCT(:).bg_threshold]))
%     msgbox('There are files without assigned background threshold!');
elseif any(isnan([INFOSTRUCT(:).xy_scale]))
    msgbox('There are files without assigned XY scale!');
else
    close gcf;  
    run_roi_detection(INFOSTRUCT, []);
%     mtrx = Threshold_test_run_roi_detection(INFOSTRUCT, []);
end

end

