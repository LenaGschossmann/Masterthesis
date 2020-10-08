function ini_ctrl_box_timescan()

% Declare globally shared variables
global SCRSZ FNAME FONTSIZE  SCMAP  WINSZ...
    PLOTRANGE CLIMUI THRESHOLDTYPE THRESHOLD PEAKTHRESHOLD SMTHWIN...
    figINFO roiINFO traceINFO FTIME ...
    CURRFILE NUMFILES SAVEPARENT FULLFILENAMES ROILIST IMHW FIGCOUNTER ROICNTID...
    RANGE TOTP TRACEID ROILIST TRACEDATA ISBGSUB evINFO

% Initiate (display) variables
TRACEID = [];
showevents = false;
showalltr = false;
singlesave = true;
detectall = false;
importdata = [];
importroi = [];
evINFO = create_evInfo();

whctrl = [410 580];
vspace = [10 20 30]; % [small bigg bigger]
hspace = [5 15]; % [ctr-aligned rest]
hctr = round(whctrl(1)/2);
whin = [70 20];
whinbg = [90 20];
whbutsm = [80 20];
whbutbg = [120 25];
whdd = [200 15];
whtxt = [200 20];
whtxtbg = [whctrl(1)-2*hspace(2) 20];
whrb = [50 20];

% Display positions
positionctrl = [SCRSZ(1)-round(1.5*whctrl(1)) round(SCRSZ(2)/2-0.5*whctrl(2)) whctrl];
savetxtpos = [hctr-hspace(1)-whtxt(1)/2 vspace(2) whtxt(1)/2 whtxt(2)];
saveselbutpos = [hctr+hspace(1) savetxtpos(2) round(whbutsm(1)*0.65) whbutbg(2)];
saveallbutpos = [saveselbutpos(1)+saveselbutpos(3)+hspace(1) savetxtpos(2) round(whbutsm(1)*0.65) whbutbg(2)];
detevbutpos = [hctr-whbutbg(1)-hspace(1) savetxtpos(2)+savetxtpos(4)+vspace(2) whbutbg];
detevallrbpos = [hctr+hspace(1) detevbutpos(2) whrb];
detevsaverbpos = [detevallrbpos(1)+detevallrbpos(3)+hspace(1) detevbutpos(2) round(whrb(1)*2) whrb(2)];
disptrbutpos = [hctr-whbutbg(1)-hspace(1) detevbutpos(2)+detevbutpos(4)+vspace(2) whbutbg];
dispalltrrbpos = [hctr+hspace(1) disptrbutpos(2) whrb];
dispevtrrbpos = [dispalltrrbpos(1)+dispalltrrbpos(3)+hspace(1) disptrbutpos(2) round(whrb(1)*2) whrb(2)];
threshtxtpos = [hspace(2) disptrbutpos(2)+whbutbg(2)+vspace(2) whtxt];
threshinpos = [hctr+hspace(1) threshtxtpos(2) whin];
peakinpos = [threshinpos(1)+whin(1)+hspace(1) threshtxtpos(2) whin];
threshtypetxtpos = [hspace(2) peakinpos(2)+whin(2)+vspace(1) whtxt];
threshtyperbpos = [hctr+hspace(1) threshtypetxtpos(2) whrb;...
    hctr+hspace(1)+whrb(1)+hspace(1) threshtypetxtpos(2) whrb];
trcsmthtxtpos = [hspace(2) threshtypetxtpos(2)+whtxt(2)+vspace(1) whtxt];
trcsmthinpos = [hctr+hspace(1) trcsmthtxtpos(2) whin];
trcsmthbutpos = [trcsmthinpos(1)+trcsmthinpos(3)+hspace(1) trcsmthinpos(2) whbutsm(1)/3 whbutsm(2)];
rangetxtpos = [hspace(2) trcsmthinpos(2)+trcsmthinpos(4)+vspace(1) whtxt];
rangein1pos = [hctr+hspace(1) rangetxtpos(2) whin];
rangein2pos = [hctr+hspace(1)*2+whin(1) rangetxtpos(2) whin];
backsubutpos = [hctr+hspace(1) rangein1pos(2)+rangein1pos(4)+vspace(2) whbutbg];
backsubtxtpos = [hctr-hspace(1)-whtxt(1)*0.8 backsubutpos(2) whtxt(1)*0.8 whtxt(2)*1.5];
roilistlbpos = [hctr-round(whdd(1)/3) backsubtxtpos(2)+backsubtxtpos(4)+vspace(2) whdd(1) whdd(2)*8];
fttxtpos = [hspace(2)*2 roilistlbpos(2)+roilistlbpos(4)-whin(2) whin];
ftinpos = [fttxtpos(1) fttxtpos(2)-whin(2) whin];
openbutpos = [hctr-whbutbg(1)/2 fttxtpos(2)+fttxtpos(4)+vspace(2) whbutbg];
nametxtpos = [hspace(2) openbutpos(2)+openbutpos(4)+vspace(1) whtxtbg];

% Main window
ctrlfig = figure('Name', 'TimeSeries Analysis', 'Position', positionctrl, 'toolbar', 'none', 'menu', 'none');
set(gca, 'Color', 'none'); axis off;
title('Control Box', 'FontSize', round(FONTSIZE*1.3));

% Add UI controls
savetxt = uicontrol('parent', ctrlfig, 'style', 'text','position', savetxtpos,'string', 'Save combined:','FONTSIZE', FONTSIZE);
saveselbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', saveselbutpos,'string', 'Selected','FONTSIZE', FONTSIZE, 'callback', {@cb_savebut});
saveallbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', saveallbutpos,'string', 'All','FONTSIZE', FONTSIZE, 'callback', {@cb_savebut});

detevbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton','position', detevbutpos,'string', 'Events selected traces','FONTSIZE', FONTSIZE, 'callback', {@cb_detevbut});
detevallrb = uicontrol('parent', ctrlfig, 'style', 'radiobutton', 'position', detevallrbpos,'string', 'All', 'Value', 0,'FONTSIZE', FONTSIZE, 'callback', {@cb_detevallrb});
detevsaverb = uicontrol('parent', ctrlfig, 'style', 'radiobutton', 'position', detevsaverbpos,'string', 'Save single', 'Value', 1,'FONTSIZE', FONTSIZE, 'callback', {@cb_detevsaverb});

disptrbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', disptrbutpos,'string', 'Show selected traces','FONTSIZE', FONTSIZE, 'callback', {@cb_disptrbut});
dispalltrrb = uicontrol('parent', ctrlfig, 'style', 'radiobutton', 'position', dispalltrrbpos,'string', 'All', 'Value', 0,'FONTSIZE', FONTSIZE, 'callback', {@cb_dispalltrrb});
dispevtrrb = uicontrol('parent', ctrlfig, 'style', 'radiobutton', 'position', dispevtrrbpos,'string', 'Show events', 'Value', 0,'FONTSIZE', FONTSIZE, 'FOREGROUNDCOLOR', [0.5 0.5 0.5],'Enable', 'inactive', 'callback', {@cb_dispevtrrb});

threshtxt = uicontrol('parent', ctrlfig, 'style', 'text','position', threshtxtpos,'string', 'Threshold [s.d.]: start | peak','FONTSIZE', FONTSIZE);
threshin = uicontrol('parent', ctrlfig, 'style', 'edit','position', threshinpos,'FONTSIZE', FONTSIZE, 'String', num2str(THRESHOLD),'Callback', {@cb_threshin});
peakin = uicontrol('parent', ctrlfig, 'style', 'edit','position', peakinpos, 'Enable', 'off', 'FONTSIZE', FONTSIZE, 'String', num2str(PEAKTHRESHOLD),'Callback', {@cb_peakthreshin});

threshtypetxt1 = uicontrol('parent', ctrlfig, 'style', 'text','position', threshtypetxtpos,'string', 'Threshold type:','FONTSIZE', FONTSIZE);
threshtyperb1 = uicontrol('parent', ctrlfig, 'style', 'radiobutton', 'position', threshtyperbpos(1,:),'string', 'delta F', 'Value', 1,'FONTSIZE', FONTSIZE, 'callback', {@cb_threshtyperb});
threshtyperb2 = uicontrol('parent', ctrlfig, 'style', 'radiobutton', 'position', threshtyperbpos(2,:),'string', 'sd', 'Value', 0,'FONTSIZE', FONTSIZE, 'callback', {@cb_threshtyperb});

trcsmthtxt = uicontrol('parent', ctrlfig, 'style', 'text','position', trcsmthtxtpos,'string', 'Boxcar smoothing window:','FONTSIZE', FONTSIZE);
trcsmthin = uicontrol('parent', ctrlfig, 'style', 'edit','position', trcsmthinpos,'FONTSIZE', FONTSIZE, 'String', num2str(SMTHWIN),'Callback', {@cb_trcsmthin});
trcsmthbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', trcsmthbutpos,'string', 's','FONTSIZE', FONTSIZE, 'callback', {@cb_trcsmthbut});

rangetxt = uicontrol('parent', ctrlfig, 'style', 'text','position', rangetxtpos,'string', 'Data point range:','FONTSIZE', FONTSIZE);
rangein1 = uicontrol('parent', ctrlfig, 'style', 'edit','position', rangein1pos,'FONTSIZE', FONTSIZE, 'String', '','Callback', {@cb_rangein1});
rangein2 = uicontrol('parent', ctrlfig, 'style', 'edit','position', rangein2pos,'FONTSIZE', FONTSIZE, 'String', '','Callback', {@cb_rangein2});

backsubut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', backsubutpos,'string', 'RAW > Subtract','FONTSIZE', FONTSIZE, 'callback', {@cb_backsubut});
backsubtxt = uicontrol('parent', ctrlfig, 'style', 'text','position', backsubtxtpos,'string', {'Subtract background:', 'select ROI as background!'},'FONTSIZE', FONTSIZE);

ftin = uicontrol('parent', ctrlfig, 'style', 'edit', 'position', ftinpos,'string', num2str(FTIME),'FONTSIZE', FONTSIZE, 'backgroundcolor', [0.85 0.85 0.85],'callback', {@cb_ftin});
fttxt = uicontrol('parent', ctrlfig, 'style', 'text','position', fttxtpos,'string', 'ms/frame:','FONTSIZE', FONTSIZE);
roilistlb = uicontrol('parent', ctrlfig, 'style', 'listbox','position', roilistlbpos,'FONTSIZE', FONTSIZE, 'string',ROILIST,'Callback', {@cb_roilistlb});
openbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', openbutpos,'string', 'Open file','FONTSIZE', FONTSIZE, 'callback', {@cb_openbut});
nametxt = uicontrol('parent', ctrlfig, 'style', 'text','position', nametxtpos,'string', FNAME,'FONTSIZE', FONTSIZE);

uiwait;

%% Local callbacks
    function cb_openbut(~,~)
        [filenames, pathnames] = uigetfile({'*.csv', '*.xlsx'}, 'Multiselect', 'off'); % Select (multiple) files for processing
        FNAME = filenames;
        if isa(filenames,'cell') || isa(filenames,'char')
            if ~isa(filenames,'cell'), filenames = {filenames}; end
            if isempty(SAVEPARENT), SAVEPARENT = uigetdir(pathnames, 'Parent folder for saving'); end % Select path for parent folder of saved files
            FULLFILENAMES = fullfile(pathnames, filenames);
            NUMFILES = size(filenames,2);
            filepointer = FULLFILENAMES{CURRFILE};
            datatable = readtable(filepointer);
            importroi = datatable.Properties.VariableNames;                   % ROI references for TRACEDATA
            importdata = table2array(datatable);                             % Data to work with (is subject to changes)
            if unique(diff(importdata(:,1))) == 1                            % Delete first column if it carries only row ids
                importdata = importdata(:,2:end);
                importroi = importroi(2:end);
            end 
            for iTr = 1:numel(importroi), importroi{iTr}= strtrim(strrep(importroi{iTr},'_','')); end
            TRACEDATA = importdata; ROILIST = importroi;
%             ft = inputdlg('Frame time in ms:', 'Frame Time', [1 40]); ft = str2double(ft{1}); 
%             if ~isempty(ft), FTIME = ft/1000; end % Convert to sec
            RANGE = [1 size(TRACEDATA,1)]; TOTP = RANGE(2);
            rangein1.String = num2str(RANGE(1)); rangein2.String = num2str(RANGE(2)); 
            roilistlb.String = ROILIST;
            TRACEID = 1;
            nametxt.String = FNAME;
            evINFO = create_evInfo();
        end
    end

    function cb_ftin(hObj,~)
        ft = str2double(hObj.String);
        if ~isempty(ft), FTIME = ft/1000; end
    end

    function cb_roilistlb(hObj,~)
        TRACEID = hObj.Value;
        roilistlb.Max = numel(ROILIST);
    end

    function cb_backsubut(hObj,~)
        if isempty(TRACEID)
            msgbox('Select a ROI for background subtraction');
        elseif numel(TRACEID) > 1
            msgbox('Select 1 ROI for background subtraction');
        else
            if strncmp(hObj.String,'RAW',3)
                % Subtract background
                ROILIST(TRACEID) = [];
                bgdata = TRACEDATA(:,TRACEID);
                TRACEDATA(:,TRACEID) = [];
                TRACEDATA = TRACEDATA-bgdata;
                
                hObj.String = 'Bg. subtracted';
                hObj.BackgroundColor = '#77AC30';
                roilistlb.Value = 1;
                roilistlb.String = ROILIST;
                ISBGSUB = true;
            else
                hObj.String = 'RAW > Subtract';
                hObj.BackgroundColor = [0.94 0.94 0.94];
                ROILIST = importroi;                                        % Set back to roilist with background trace included
                TRACEDATA = importdata;                                     % Set back to data without background subtraction
                roilistlb.Value = 1;
                roilistlb.String = ROILIST;
                ISBGSUB = false;
            end
        end
    end

    function cb_rangein1(hObj,~)
        RANGE(1) = str2double(get(hObj,'String'));
        if RANGE(1) > TOTP || RANGE(1) >= RANGE(2), RANGE(1) = 1; set(hObj,'String', num2str(RANGE(1))); end
    end

    function cb_rangein2(hObj,~)
        RANGE(2) = str2double(get(hObj,'String'));
        if RANGE(2) > TOTP, RANGE(2) = TOTP; set(hObj,'String', num2str(RANGE(2))); end
        if RANGE(2) <= RANGE(1), RANGE(2) = RANGE(1); set(hObj,'String', num2str(RANGE(2))); end
    end

    function cb_threshtyperb(hObj,~)
        if strcmp(hObj.String, 'delta F')
            THRESHOLDTYPE = 'dF';
            threshtyperb2.Value = 0;
            peakin.Enable = 'off';
        else
            THRESHOLDTYPE = 'sd';
            threshtyperb1.Value = 0;
            peakin.Enable = 'on';
        end
    end

    function cb_trcsmthin(hObj,~)
        if strcmp(trcsmthbut.String,'s'), SMTHWIN = str2double(get(hObj,'String'));
        else, SMTHWIN = str2double(get(hObj,'String'))*FTIME;
        end
    end

    function cb_trcsmthbut(hObj,~)
        if strcmp(hObj.String, 's'), hObj.String = 'f'; SMTHWIN = str2double(trcsmthin.String)*FTIME;
        else, hObj.String = 's'; SMTHWIN = str2double(trcsmthin.String);
        end
    end

    function cb_threshin(hObj,~)
        THRESHOLD = str2double(get(hObj,'String'));
    end

    function cb_peakthreshin(hObj,~)
        PEAKTHRESHOLD = str2double(get(hObj,'String'));
    end

    function cb_disptrbut(~,~)
        if isempty(TRACEID) && ~showalltr
            msgbox('Please select a ROI first!');
        else
            traces_overview(showalltr, showevents);
        end
    end

    function cb_dispalltrrb(hObj,~)
        showalltr = logical(hObj.Value);
    end

    function cb_dispevtrrb(hObj,~)
        showevents = logical(hObj.Value);
    end

    function cb_detevbut(~,~)
        if isempty(TRACEID) && ~showalltr
            msgbox('Please select a ROI first!');
        else
            ev_detection_overview(detectall); %detect_events(detectall)
            dispevtrrb.Enable = 'on'; dispevtrrb.ForegroundColor = [0 0 0];
        end
    end

    function cb_detevallrb(hObj,~)
        detectall = logical(hObj.Value);
    end

    function cb_detevsaverb(hObj,~)
        singlesave = logical(hObj.Value);
    end

    function cb_savebut()
        
    end

    function cb_savengobut(~,~)
        figures = get(groot,'Children');
        if isempty(roiINFO(1).position)
            msgbox('You have to select a ROI first!');
        else
            findfig = strncmp({figures(:).Name}, 'Fig',3);
            currfig = figures(find(findfig,1, 'first'));
            figid = currfig.UserData;
            [~,figidx] = find([figINFO(:).IDs] == figid);
            % Save
            savedir = strcat(SAVEPARENT, '\', FNAME, '\');
            if ~exist(savedir,'dir'), mkdir(savedir); end
            savepointer = strcat(savedir,FNAME,'_');
            for iE = 1: size(traceINFO,2), traceINFO(iE).frametime = FTIME; end
            save(strcat(savepointer, '_figure_info', '.mat'), 'figINFO');
            save(strcat(savepointer, '_roi_info', '.mat'), 'roiINFO');
            save(strcat(savepointer, '_traces_info', '.mat'), 'traceINFO');
            save_files(savepointer, figidx, 'all');
            
            % Open next file
            CURRFILE = CURRFILE +1;
            if CURRFILE > NUMFILES
                disp('All files have been handled.');
                uiresume(); close all;
            else
                figures = get(groot,'Children');
                tmpfig = figures(~strncmp({figures(:).Name}, 'Line',4));
                for iF = 1:numel(tmpfig), close(tmpfig(iF)); end
                filepointer = FULLFILENAMES{CURRFILE};
                open_file(filepointer, rangein1, rangein2);
            end
        end
        roiINFO = struct(); roiINFO(1).name = []; roiINFO(1).mask = []; roiINFO(1).position = []; roiINFO(1).ID = []; roiINFO(1).selected = []; roiINFO(1).saved = []; roiINFO(1).mode = []; roiINFO(1).PLOTRANGE = [];
        figINFO = struct('IDs', [], 'name', [], 'PLOTRANGE',[], 'cSCMAP',[], 'csclimits', [], 'avwinsize',[], 'saved',[]);
        traceINFO = struct('figID',[], 'fig_params',[],'roiID',[], 'binned_roi_av',[],'dFoF_roi_av',[], 'timestamp',[], 'save',[], 'currmode',[], 'showtot', []);
        FIGCOUNTER = 0;
        ROICNTID = 1;
    end

end