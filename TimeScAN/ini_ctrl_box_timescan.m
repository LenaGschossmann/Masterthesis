function ini_ctrl_box_timescan()

%% Declare globally shared variables
global SCRSZ FNAME FONTSIZE SMTHWIN SAVEPATH...
    FTIME FTIMEPREV SAVEPARENT FULLFILENAMES ROILIST...
    RANGE TOTP TRACEID TRACEDATA evINFO FTIMEVEC XTYPE...
    BASEPOINTS PREPOINTS POSTPOINTS LOCAREAINFO SOMAINFO PXCLASSCNT...

%% Initiate (display) variables
TRACEID = [];
showevents = false;
showalltr = false;
singlesave = false;
detectall = false;
importgraydata = [];
importgrayrois = cell(0);
if strcmp(XTYPE, 's'), txtcutin1 = num2str(PREPOINTS*FTIME); txtcutin2 = num2str(POSTPOINTS*FTIME);
else, txtcutin1 = num2str(PREPOINTS); txtcutin2 = num2str(POSTPOINTS);
end
% checkicon = imread('check_icon.png');

whctrl = [410 560];
vspace = [10 20 30]; % [small bigg bigger]
hspace = [5 15]; % [ctr-aligned rest]
hctr = round(whctrl(1)/2);
whin = [70 20];
whinbg = [90 20];
whbutsm = [80 20];
whbutbg = [100 25];
whdd = [200 15];
whtxt = [200 20];
whtxtbg = [whctrl(1)-2*hspace(2) 20];
whrb = [50 20];
% whicon = [20 20];

% Display positions
positionctrl = [SCRSZ(1)-round(1.5*whctrl(1)) round(SCRSZ(2)/2-0.5*whctrl(2)) whctrl];
savetxtpos = [hspace(2) vspace(2) whtxt(1)/2 whtxt(2)];
saveselbutpos = [savetxtpos(1)+savetxtpos(3)+hspace(1) savetxtpos(2) round(whbutsm(1)*0.65) whbutbg(2)];
saveallbutpos = [saveselbutpos(1)+saveselbutpos(3)+hspace(1) savetxtpos(2) round(whbutsm(1)*0.65) whbutbg(2)];
summbutpos = [saveallbutpos(1)+saveallbutpos(3)+2*hspace(2) savetxtpos(2) whbutbg];
detevbutpos = [hctr-whbutbg(1)-hspace(1) savetxtpos(2)+savetxtpos(4)+vspace(2) whbutbg];
detevallrbpos = [hctr+hspace(1) detevbutpos(2) whrb];
detevsaverbpos = [detevallrbpos(1)+detevallrbpos(3)+hspace(1) detevbutpos(2) round(whrb(1)*2) whrb(2)];
disptrbutpos = [hctr-whbutbg(1)-hspace(1) detevbutpos(2)+detevbutpos(4)+vspace(1) whbutbg];
dispalltrrbpos = [hctr+hspace(1) disptrbutpos(2) whrb];
dispevtrrbpos = [dispalltrrbpos(1)+dispalltrrbpos(3)+hspace(1) disptrbutpos(2) round(whrb(1)*2) whrb(2)];
baseptxtpos = [hspace(2) disptrbutpos(2)+disptrbutpos(4)+vspace(1) whtxt];
basepinpos = [hctr+hspace(1) baseptxtpos(2) whin];
cuttxtpos = [hspace(2) baseptxtpos(2)+baseptxtpos(4)+vspace(1) whtxt];
cut1inpos = [hctr+hspace(1) cuttxtpos(2) whin];
cut2inpos = [cut1inpos(1)+whin(1)+hspace(1) cuttxtpos(2) whin];
trcsmthtxtpos = [hspace(2) cuttxtpos(2)+whtxt(2)+vspace(1) whtxt];
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
tsmetabutpos = [fttxtpos(1) ftinpos(2)-whbutsm(2)-vspace(1) whin(1) whbutsm(2)];
% tsmetacheckpos = [fttxtpos(1)+fttxtpos(3)+hspace(1)/3 fttxtpos(2) whicon];
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

summbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton','position', summbutpos,'string', 'Summary','FONTSIZE', FONTSIZE, 'callback', {@cb_summbut});

detevbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton','position', detevbutpos,'string', 'Detect events','FONTSIZE', FONTSIZE, 'callback', {@cb_detevbut});
detevallrb = uicontrol('parent', ctrlfig, 'style', 'radiobutton', 'position', detevallrbpos,'string', 'All', 'Value', detectall,'FONTSIZE', FONTSIZE, 'callback', {@cb_detevallrb});
% detevsaverb = uicontrol('parent', ctrlfig, 'style', 'radiobutton', 'position', detevsaverbpos,'string', 'Save single', 'Value', singlesave,'FONTSIZE', FONTSIZE, 'callback', {@cb_detevsaverb});

disptrbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', disptrbutpos,'string', 'Show traces','FONTSIZE', FONTSIZE, 'callback', {@cb_disptrbut});
dispalltrrb = uicontrol('parent', ctrlfig, 'style', 'radiobutton', 'position', dispalltrrbpos,'string', 'All', 'Value', showalltr,'FONTSIZE', FONTSIZE, 'callback', {@cb_dispalltrrb});
dispevtrrb = uicontrol('parent', ctrlfig, 'style', 'radiobutton', 'position', dispevtrrbpos,'string', 'Show events', 'Value', 0,'FONTSIZE', FONTSIZE, 'FOREGROUNDCOLOR', [0.5 0.5 0.5],'Enable', 'inactive', 'callback', {@cb_dispevtrrb});

baseptxt = uicontrol('parent', ctrlfig, 'style', 'text','position', baseptxtpos,...
    'string', sprintf('# of baseline points [%s]:', XTYPE),'FONTSIZE', FONTSIZE);
basepin = uicontrol('parent', ctrlfig, 'style', 'edit','position', basepinpos,'FONTSIZE', FONTSIZE, 'String', num2str(BASEPOINTS*FTIME),'Callback', {@cb_basepin});

cuttxt = uicontrol('parent', ctrlfig, 'style', 'text','position', cuttxtpos,'string', sprintf('Cut-out before & after event [%s]:', XTYPE),'FONTSIZE', FONTSIZE);
cut1in = uicontrol('parent', ctrlfig, 'style', 'edit','position', cut1inpos,'FONTSIZE', FONTSIZE, 'String', txtcutin1,'Callback', {@cb_cut1in});
cut2in = uicontrol('parent', ctrlfig, 'style', 'edit','position', cut2inpos,'FONTSIZE', FONTSIZE, 'String', txtcutin2,'Callback', {@cb_cut2in});

trcsmthtxt = uicontrol('parent', ctrlfig, 'style', 'text','position', trcsmthtxtpos,'string', 'Boxcar smoothing window:','FONTSIZE', FONTSIZE);
trcsmthin = uicontrol('parent', ctrlfig, 'style', 'edit','position', trcsmthinpos,'FONTSIZE', FONTSIZE, 'String', num2str(SMTHWIN),'Callback', {@cb_trcsmthin});
trcsmthbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', trcsmthbutpos,'string', 's','FONTSIZE', FONTSIZE, 'callback', {@cb_trcsmthbut});

rangetxt = uicontrol('parent', ctrlfig, 'style', 'text','position', rangetxtpos,'string', 'Data point range:','FONTSIZE', FONTSIZE);
rangein1 = uicontrol('parent', ctrlfig, 'style', 'edit','position', rangein1pos,'FONTSIZE', FONTSIZE, 'String', '','Callback', {@cb_rangein1});
rangein2 = uicontrol('parent', ctrlfig, 'style', 'edit','position', rangein2pos,'FONTSIZE', FONTSIZE, 'String', '','Callback', {@cb_rangein2});

backsubut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', backsubutpos,'string', 'RAW > Subtract','FONTSIZE', FONTSIZE, 'callback', {@cb_backsubut});
backsubtxt = uicontrol('parent', ctrlfig, 'style', 'text','position', backsubtxtpos,'string', {'Subtract background:', 'select ROI as background!'},'FONTSIZE', FONTSIZE);

tsmetabut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', tsmetabutpos,'string', 'Import time','FONTSIZE', FONTSIZE, 'callback', {@cb_tsmetabut});
ftin = uicontrol('parent', ctrlfig, 'style', 'edit', 'position', ftinpos,'string', num2str(FTIME*1000),'Enable','on','FONTSIZE', FONTSIZE, 'backgroundcolor', [0.85 0.85 0.85],'callback', {@cb_ftin});
fttxt = uicontrol('parent', ctrlfig, 'style', 'text','position', fttxtpos,'string', 'ms/frame:','FONTSIZE', FONTSIZE);
roilistlb = uicontrol('parent', ctrlfig, 'style', 'listbox','position', roilistlbpos,'FONTSIZE', FONTSIZE, 'string',ROILIST,'Callback', {@cb_roilistlb});
openbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', openbutpos,'string', 'Open file','FONTSIZE', FONTSIZE, 'callback', {@cb_openbut});
nametxt = uicontrol('parent', ctrlfig, 'style', 'text','position', nametxtpos,'string', FNAME,'FONTSIZE', FONTSIZE);
% sp1 = subplot('position',tsmetacheckpos, 'parent', ctrlfig);
% imshow(checkicon, 'parent', sp1, 'magnification', 0.1);
uiwait;

%% Local callbacks
    function cb_openbut(~,~)
        [filenames, pathnames] = uigetfile({'*.csv'; '*.xlsx'}, 'Select intensity, location, dendritic count data','Multiselect', 'on');
        if isa(filenames,'char'), filenames = {filenames}; end
        if isa(filenames,'cell')
            if numel(filenames) > 1 % Location file was selected
                matchloc = cellfun(@isempty, regexp(filenames, 'loc_ar', 'match'));
                matchcnt = cellfun(@isempty, regexp(filenames, 'dend_px', 'match'));
                idLoc = find(matchloc==0);
                idCnt = find(matchcnt==0);
                idObs = find(matchcnt==1 & matchloc == 1);
                FNAME = filenames{idObs};
                % Open file with location and area information(col1: area in �m, col2-4: X, Y, Z
                if ~isempty(idLoc)
                    nameloc = fullfile(pathnames, filenames{idLoc});
                    if isa(nameloc, 'cell'), nameloc = nameloc{1}; end
                    importlocarea = readtable(nameloc);
                    importlocarea = table2array(importlocarea);
                    % Ask for table with location of somata and load
                    defpath = split(pathnames, '\');
                    defpath = fullfile(defpath{1:end-2}, '*.csv; *.xlsx');
                    [somafile, somapath] = uigetfile(defpath, 'Select table with location of somata');
                    if ~isempty(somafile) && isa(somafile, 'char')
                        somaname = fullfile(somapath, somafile);
                        if isa(somaname, 'cell'), somaname  = somaname{1}; end
                        SOMAINFO = readtable(somaname);
                        SOMAINFO = table2array(SOMAINFO);
                    else, SOMAINFO = [];
                    end
                    clear('locname', 'somaname');
                else
                    importlocarea = [];
                end
                if ~isempty(idCnt)
                    namecnt = fullfile(pathnames, filenames{idCnt});
                    if isa(namecnt, 'cell'), namecnt = namecnt{1}; end
                    PXCLASSCNT = readtable(namecnt);
                    PXCLASSCNT = table2array(PXCLASSCNT); % col1: non-dendritic, col2: dendritic
                end
            else
                FNAME = filenames{1};
                SOMAINFO =  [];
                importlocarea = [];
            end
            if isempty(SAVEPARENT), SAVEPARENT = uigetdir(pathnames, 'Parent folder for saving'); end % Select path for parent folder of saved files
            FULLFILENAMES = fullfile(pathnames, FNAME);
            if isa(FULLFILENAMES, 'cell'), FULLFILENAMES = FULLFILENAMES{1}; end
            [~,f,~] = fileparts(FULLFILENAMES);
            SAVEPATH = strcat(SAVEPARENT,'\',f);
            if ~exist(SAVEPATH,'dir'), mkdir(SAVEPATH); end
            
            % Read intensity data
            datatable = readtable(FULLFILENAMES);
            datacolnames = datatable.Properties.VariableNames;              % ROI references for TRACEDATA
            datarray = table2array(datatable); clear('datatable');
            
            for ii = 1:numel(datacolnames)                                  % Remove any column with information that is not the mean intensity
                datacolnames{ii} = upper(datacolnames{ii});
                if strncmp(datacolnames{ii},'MEAN',4)
                    importgraydata = [importgraydata datarray(:,ii)];
                    importgrayrois = [importgrayrois datacolnames{ii}(6:end-1)]; % Strip name from "Mean()" (works for data saved from ImageJ ROI measurement)
                end
            end
            clear('datarray','datacolnames');
            
            % Sort by number
            sortidx = roi_sort(importgrayrois);
            importgrayrois = importgrayrois(sortidx);
            importgraydata = importgraydata(:,sortidx);
            if ~isempty(importlocarea), importlocarea = importlocarea(sortidx,:); end
            TRACEDATA = importgraydata; ROILIST = importgrayrois; LOCAREAINFO = importlocarea;
%             ft = inputdlg('Frame time in ms:', 'Frame Time', [1 40]); ft = str2double(ft{1}); 
%             if ~isempty(ft), FTIME = ft/1000; end % Convert to sec
            RANGE = [1 size(TRACEDATA,1)]; TOTP = RANGE(2);
            if strcmp(XTYPE,'s'),rangein1.String = num2str(RANGE(1))*FTIME; rangein2.String = num2str(RANGE(2))*FTIME;
            else, rangein1.String = num2str(RANGE(1)); rangein2.String = num2str(RANGE(2));
            end
            roilistlb.String = ROILIST;
            TRACEID = 1;
            nametxt.String = FNAME;
            evINFO = create_evInfo();
            FTIMEVEC = (1:TOTP).*FTIME;
        end
    end

    function cb_summbut(~,~)
        if any(strcmp({evINFO(:).accepted},'accepted')), create_summary(); end
    end

    function cb_ftin(hObj,~)
        ft = str2double(hObj.String);
        if ~isempty(ft)
            FTIME = ft/1000;
            FTIMEVEC = (1:TOTP).*FTIME;
        end
        if strcmp(XTYPE,'s')
            BASEPOINTS = round(str2num(basepin.String)/FTIMEPREV); basepin.String = BASEPOINTS*FTIME;
            PREPOINTS = round(str2num(cut1in.String)/FTIMEPREV); cut1in.String = PREPOINTS*FTIME;
            POSTPOINTS = round(str2num(cut2in.String)/FTIMEPREV); cut2in.String = POSTPOINTS*FTIME;
            SMTHWIN = round(str2num(trcsmthin.String)/FTIMEPREV); trcsmthin.String = SMTHWIN*FTIME;
            if ~strcmp(rangein1.String,'')
                RANGE(1) = round(str2num(rangein1.String)/FTIMEPREV);
                if RANGE(1) < 1, RANGE(1) = 1; end
                rangein1.String = RANGE(1)*FTIME;
            end
            if ~strcmp(rangein2.String,'')
                RANGE(2) = round(str2num(rangein2.String)/FTIMEPREV);
                if RANGE(2) > TOTP, RANGE(2) = TOTP; end
                rangein2.String = RANGE(2)*FTIME;
            end
        elseif strcmp(XTYPE,'f')
            BASEPOINTS = round(str2num(basepin.String)*FTIMEPREV/FTIME); basepin.String = BASEPOINTS;
            PREPOINTS = round(str2num(cut1in.String)*FTIMEPREV/FTIME); cut1in.String = PREPOINTS;
            POSTPOINTS = round(str2num(cut2in.String)*FTIMEPREV/FTIME); cut2in.String = POSTPOINTS;
            SMTHWIN = round(str2num(trcsmthin.String)*FTIMEPREV/FTIME); trcsmthin.String = SMTHWIN;
            if ~strcmp(rangein1.String,'')
                RANGE(1) = round(str2num(rangein1.String)*FTIMEPREV/FTIME);
                if RANGE(1) < 1, RANGE(1) = 1; end
                rangein1.String = RANGE(1);
            end
            if ~strcmp(rangein2.String,'')
                RANGE(2) = round(str2num(rangein2.String)*FTIMEPREV/FTIME);
                if RANGE(2) > TOTP, RANGE(2) = TOTP; end
                rangein2.String = RANGE(2);
            end
        end
        FTIMEPREV = FTIME;
    end

    function cb_tsmetabut(~,~)
        [metafn, metapn] = uigetfile({'*.nd2'}, 'Multiselect', 'off');
        metapointer = fullfile(metapn, metafn);
        [~, imInfo] = load_bf_file(metapointer, true);
        immeta =imInfo.metadata;
        timestrings = immeta(strncmp(immeta, 'timestamp',9));
        timesplit = cellfun(@(x) regexp(x,'= ','split'), timestrings,'UniformOutput', false);
        sortsplit = cellfun(@(x) x{1}, timesplit, 'UniformOutput', false);
        timesplit = cellfun(@(x) str2double(x{2}), timesplit, 'UniformOutput', false);
        sortsplit = cellfun(@(x) regexp(x,'#','split'), sortsplit,'UniformOutput', false);
        sortsplit = cellfun(@(x) str2double(x{2}), sortsplit, 'UniformOutput', false);
        sortsplit = cell2mat(sortsplit);
        [~,sortidx] = sort(sortsplit, 'ascend');
        timeseries = cell2mat(timesplit); timeseries = timeseries(sortidx);
        format long;
        FTIMEVEC = timeseries;
        tsmetabut.BackgroundColor = '#77AC30';
        ftin.Enable = 'Off';
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
                if ~isempty(LOCAREAINFO), LOCAREAINFO(TRACEID,:) = []; end
                bgdata = TRACEDATA(:,TRACEID);
                bgdata = mean(bgdata);
                TRACEDATA(:,TRACEID) = [];
                TRACEDATA = TRACEDATA-bgdata;
                
                hObj.String = 'Bg. subtracted';
                hObj.BackgroundColor = '#77AC30';
                roilistlb.Value = 1;
                roilistlb.String = ROILIST;
                ISBGSUB = true;
                evINFO = create_evInfo();
            else
                hObj.String = 'RAW > Subtract';
                hObj.BackgroundColor = [0.94 0.94 0.94];
                ROILIST = importgrayrois;                                        % Set back to roilist with background trace included
                if ~isempty(LOCAREAINFO), LOCAREAINFO = importlocarea; end
                TRACEDATA = importgraydata;                                     % Set back to data without background subtraction
                roilistlb.Value = 1;
                roilistlb.String = ROILIST;
                ISBGSUB = false;
                evINFO = create_evInfo();
            end
        end
    end

    function cb_rangein1(hObj,~)
        RANGE(1) = str2double(get(hObj,'String'));
        if RANGE(1) > TOTP || RANGE(1) >= RANGE(2), RANGE(1) = 1; set(hObj,'String', num2str(RANGE(1))); end
    end

    function cb_rangein2(hObj,~)
        RANGE(2) = str2double(get(hObj,'String'));
        if RANGE(2) > TOTP, RANGE(2) = TOTP; end
        if RANGE(2) <= RANGE(1), RANGE(2) = RANGE(1); end
        if strcmp(XTYPE,'s'), RANGE(2) = round(str2double(get(hObj,'String'))/FTIME);
        else, RANGE(2) = str2double(get(hObj,'String'));
        end
        set(hObj,'String', num2str(RANGE(2)));
    end

    function cb_cut1in(hObj,~)
        if ~isempty(hObj.String), PREPOINTS = str2double(hObj.String); end
        if strcmp(XTYPE, 's')
            PREPOINTS = round(PREPOINTS/FTIME);
            cut1in.String = num2str(PREPOINTS*FTIME);
        end
    end

    function cb_cut2in(hObj,~)
        if ~isempty(hObj.String), POSTPOINTS = str2double(hObj.String); end
        if strcmp(XTYPE, 's')
            POSTPOINTS = round(POSTPOINTS/FTIME);
            cut2in.String = num2str(POSTPOINTS*FTIME);
        end
    end

    function cb_basepin(hObj,~)
        if ~isempty(hObj.String), BASEPOINTS = str2double(hObj.String); end
        if strcmp(XTYPE, 's'), BASEPOINTS = round(BASEPOINTS/FTIME); end
    end

    function cb_trcsmthin(hObj,~)
        if strcmp(XTYPE,'s'), SMTHWIN = round(str2double(get(hObj,'String'))/FTIME);
        else, SMTHWIN = str2double(get(hObj,'String'));
        end
    end

    function cb_trcsmthbut(hObj,~)
        if strcmp(hObj.String, 's')
            hObj.String = 'f';
            SMTHWIN = str2double(trcsmthin.String);
            XTYPE = 'f';
            TRCXLABEL = 'frames';
            basepin.String = num2str(str2double(basepin.String)/FTIME);
            trcsmthin.String = num2str(str2double(trcsmthin.String)/FTIME);
            if ~strcmp(rangein1.String,''), rangein1.String = num2str(RANGE(1)); end
            if ~strcmp(rangein2.String,''), rangein2.String = num2str(RANGE(2)); end
            cut1in.String = num2str(PREPOINTS); cut2in.String = num2str(POSTPOINTS);
        else
            hObj.String = 's';
            SMTHWIN = str2double(trcsmthin.String)/FTIME;
            XTYPE = 's';
            TRCXLABEL = 'time [s]';
            basepin.String = num2str(str2double(basepin.String)*FTIME);
            trcsmthin.String = num2str(str2double(trcsmthin.String)*FTIME);
            if ~strcmp(rangein1.String,''), rangein1.String = num2str(RANGE(1)*FTIME); end
            if ~strcmp(rangein2.String,''), rangein2.String = num2str(RANGE(2)*FTIME); end
            cut1in.String = num2str(PREPOINTS*FTIME); cut2in.String = num2str(POSTPOINTS*FTIME);
        end
        baseptxt.String = sprintf('# of baseline points [%s]:', XTYPE);
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
        if any(strcmp({evINFO(:).accepted},'accepted')), showevents = logical(hObj.Value);
        else, msgbox('Run event detection first, review and accept events!');
        end
    end

    function cb_detevbut(~,~)
        if isempty(TRACEID) && ~detectall
            msgbox('Please select a ROI first!');
        else
            ev_detection_overview(detectall, singlesave); %detect_events(detectall)
            dispevtrrb.Enable = 'on'; dispevtrrb.ForegroundColor = [0 0 0];
        end
    end

    function cb_detevallrb(hObj,~)
        detectall = logical(hObj.Value);
    end

    function cb_detevsaverb(hObj,~)
        singlesave = logical(hObj.Value);
    end

    function cb_savebut(hObj,~)
        if strcmp(hObj.String, 'All'), saveidx = 1:numel(ROILIST);
        else, saveidx = TRACEID;
        end
        save_files_2(saveidx);
    end

    function sortidx = roi_sort(namelist)
        numlist = zeros(numel(namelist),1);
        ii=1;
        while ii <= numel(namelist)
            strnumbers = regexp(namelist{ii},{'1', '2', '3', '4', '5', '6', '7', '8', '9', '0'});
            idx = [];
            for iN = 1:numel(strnumbers), if ~isempty(strnumbers{iN}), idx = [idx strnumbers{iN}]; end, end
            if ~isempty(idx), numlist(ii) = str2num(namelist{ii}(idx)); else, numlist(ii) = NaN; end
            ii= ii+1;
        end
        if ~all(isnan(numlist)), [~,sortidx] = sort(numlist,'ascend', 'MissingPlacement', 'first');
        else, sortidx = 1:numel(namelist);
        end
    end

end