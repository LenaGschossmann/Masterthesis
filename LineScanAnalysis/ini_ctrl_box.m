function ini_ctrl_box()

% Declare globally shared variables
global SCRSZ FNAME FONTSIZE SCMAPDDITEMS SCMAP WINSZDDITEMS WINSZ...
    PLOTRANGE CLIMUI THRESHOLD PEAKTHRESHOLD SMTHWIN...
    figINFO roiINFO traceINFO FTIME...
    CURRFILE NUMFILES SAVEPARENT FULLFILENAMES IMMETA

% Initiate (display) variables
whctrl = [430 450];
vspace = [10 20 30]; % [small bigg bigger]
hspace = [5 15]; % [ctr-aligned rest]
hctr = round(whctrl(1)/2);
whin = [40 20];
whbutsm = [80 20];
whbutbg = [120 20];
whdd = [100 15];
whtxt = [200 20];
whtxtbg = [whctrl(1)-2*hspace(2) 20];

% Display positions
positionctrl = [SCRSZ(1)-round(1.5*whctrl(1)) round(SCRSZ(2)/2-0.5*whctrl(2)) whctrl];
deltrcbutpos = [hctr-whbutsm(1)*2-hspace(1)*2 vspace(2) whbutsm];
showtrcbutpos = [hctr-whbutsm(1)-hspace(1) deltrcbutpos(2) whbutsm];
savengobutpos = [hctr+hspace(1) deltrcbutpos(2) whbutsm];
nosavengobutpos = [hctr+whbutsm(1)+hspace(1)*2 deltrcbutpos(2) whbutsm];
threshtxtpos = [hspace(2) deltrcbutpos(2)+whbutbg(2)+vspace(1) whtxt];
threshinpos = [hctr+hspace(2) threshtxtpos(2) whin];
peakinpos = [threshinpos(1)+whin(1)+hspace(1) threshtxtpos(2) whin];
evsmthtxtpos = [hspace(2) threshtxtpos(2)+whtxt(2)+vspace(1) whtxt];
evsmthinpos = [hctr+hspace(2) evsmthtxtpos(2) whin];
roilbpos = [hctr+hspace(2) evsmthtxtpos(2)+whtxt(2)+vspace(1) whdd(1) whdd(2)*5];
roilbtxtpos = [hspace(2) roilbpos(2)+roilbpos(4)-whtxt(2) whtxt];
newbutpos = [hctr-whbutsm(1)-hspace(1) roilbtxtpos(2)+whtxt(2)+vspace(2) whbutsm];
updatebutpos = [hctr+hspace(1) newbutpos(2) whbutsm];
sclimtxtpos = [hspace(2) newbutpos(2)+whbutsm(2)+vspace(1) whtxt];
scliminpos = [hctr+hspace(2) sclimtxtpos(2) whin;...
    hctr+hspace(2)+whin(1)+hspace(1) sclimtxtpos(2) whin];
scmapddpos = [hctr+hspace(2) sclimtxtpos(2)+whtxt(2)+vspace(1) whdd(1) whdd(2)];
scmaptxtpos = [hspace(2) scmapddpos(2)+ scmapddpos(4)-whtxt(2) whtxt];
winszddpos = [hctr+hspace(2) scmaptxtpos(2)+whtxt(2)+vspace(1) whdd(1) whdd(2)];
winsztxtpos = [hspace(2) winszddpos(2)+ winszddpos(4)-whtxt(2) whtxt];
rangetxtpos = [hspace(2) winsztxtpos(2)+whtxt(2)+vspace(1) whtxt];
rangeinpos = [hctr+hspace(2) rangetxtpos(2) whin;...
    hctr+hspace(2)+whin(1)+hspace(1) rangetxtpos(2) whin];
openbutpos = [hctr-whbutsm(1)-hspace(1) rangetxtpos(2)+rangetxtpos(4)+vspace(1) whbutsm];
metabutpos = [hctr+hspace(1) openbutpos(2) whbutsm];
nametxtpos = [hspace(2) metabutpos(2)+metabutpos(4)+vspace(2) whtxtbg];

% Main window
ctrlfig = figure('Name', 'Line Scan Analysis', 'Position', positionctrl, 'toolbar', 'none', 'menu', 'none');
set(gca, 'Color', 'none'); axis off;
title('Control Box', 'FontSize', round(FONTSIZE*1.3));

% Add UI controls
deltrcbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', deltrcbutpos,'string', 'Delete ROI','FONTSIZE', FONTSIZE, 'callback', {@cb_deltrcbut});
showtrcbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', showtrcbutpos,'string', 'Show traces','FONTSIZE', FONTSIZE, 'callback', {@cb_showtrcbut});
savengobut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', savengobutpos,'string', 'Save & Go','FONTSIZE', FONTSIZE, 'callback', {@cb_savengobut});
nosavengobut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', nosavengobutpos,'string', 'No Save & Go','FONTSIZE', FONTSIZE, 'callback', {@cb_nosavengobut});

threshtxt = uicontrol('parent', ctrlfig, 'style', 'text','position', threshtxtpos,'string', 'Threshold [s.d.]: start | peak','FONTSIZE', FONTSIZE);
threshin = uicontrol('parent', ctrlfig, 'style', 'edit','position', threshinpos,'FONTSIZE', FONTSIZE, 'String', num2str(THRESHOLD),'Callback', {@cb_threshin});
peakin = uicontrol('parent', ctrlfig, 'style', 'edit','position', peakinpos,'FONTSIZE', FONTSIZE, 'String', num2str(PEAKTHRESHOLD),'Callback', {@cb_peakthreshin});

evsmthtxt = uicontrol('parent', ctrlfig, 'style', 'text','position', evsmthtxtpos,'string', 'Smoothing window [s]','FONTSIZE', FONTSIZE);
evsmthin = uicontrol('parent', ctrlfig, 'style', 'edit','position', evsmthinpos,'FONTSIZE', FONTSIZE, 'String', num2str(SMTHWIN),'Callback', {@cb_evsmthin});

roilbtxt  = uicontrol('parent', ctrlfig, 'style', 'text','position', roilbtxtpos,'string', 'Select ROI for analysis:','FONTSIZE', FONTSIZE);
roilb = uicontrol('parent', ctrlfig, 'style', 'listbox','position', roilbpos,'FONTSIZE', FONTSIZE, 'string', '', 'max', 20,'Callback', {@cb_roilb});

newbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', newbutpos,'string', 'New figure','FONTSIZE', FONTSIZE, 'callback', {@cb_newbut});
updatebut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', updatebutpos,'string', 'Update figure', 'FONTSIZE',FONTSIZE, 'callback', {@cb_updatebut});
sclimtxt = uicontrol('parent', ctrlfig, 'style', 'text','position', sclimtxtpos,'string', 'Set colorscale limits:','FONTSIZE', FONTSIZE);
sclimin1 = uicontrol('parent', ctrlfig, 'style', 'edit','position', scliminpos(1,:),'FONTSIZE', FONTSIZE, 'String', num2str(CLIMUI(1)),'Callback', {@cb_sclimin1});
sclimin2 = uicontrol('parent', ctrlfig, 'style', 'edit','position', scliminpos(2,:),'FONTSIZE', FONTSIZE, 'String', num2str(CLIMUI(2)),'Callback', {@cb_sclimin2});

scmaptxt = uicontrol('parent', ctrlfig, 'style', 'text','position', scmaptxtpos,'string', 'Choose LUT:','FONTSIZE', FONTSIZE);
scmapdd = uicontrol('parent', ctrlfig, 'style', 'popupmenu','position', scmapddpos,'FONTSIZE', FONTSIZE, 'string', SCMAPDDITEMS, 'Callback', {@cb_scmapdd});

winsztxt = uicontrol('parent', ctrlfig, 'style', 'text','position', winsztxtpos,'string', 'Select averaging window size:','FONTSIZE', FONTSIZE);
winszdd = uicontrol('parent', ctrlfig, 'style', 'popupmenu','position', winszddpos,'FONTSIZE', FONTSIZE, 'string', WINSZDDITEMS, 'Callback', {@cb_winszdd});

rangetxt = uicontrol('parent', ctrlfig, 'style', 'text','position', rangetxtpos,'string', 'Set plot samplepoint range:','FONTSIZE', FONTSIZE);
rangein1 = uicontrol('parent', ctrlfig, 'style', 'edit','position', rangeinpos(1,:),'FONTSIZE', FONTSIZE, 'String', num2str(PLOTRANGE(1)),'Callback', {@cb_rangein1});
rangein2 = uicontrol('parent', ctrlfig, 'style', 'edit','position', rangeinpos(2,:),'FONTSIZE', FONTSIZE, 'String', num2str(PLOTRANGE(2)),'Callback', {@cb_rangein2});

nametxt = uicontrol('parent', ctrlfig, 'style', 'text','position', nametxtpos,'string', FNAME,'FONTSIZE', FONTSIZE);
openbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', openbutpos,'string', 'Open files','FONTSIZE', FONTSIZE, 'callback', {@cb_openbut});
metabut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', metabutpos,'string', 'Metadata','FONTSIZE', FONTSIZE, 'callback', {@cb_metabut});

uiwait;

%% Local callbacks

    function cb_openbut(~,~)
        [filenames, pathnames] = uigetfile({'*.nd2'; '*.tif'}, 'Multiselect', 'on'); % Select (multiple) files for processing
        if isa(filenames,'cell') || isa(filenames,'char')
            if ~isa(filenames,'cell'), filenames = {filenames}; end
            SAVEPARENT = uigetdir(pathnames, 'Parent folder for saving');   % Select path for parent folder of saved files
            FULLFILENAMES = fullfile(pathnames, filenames);
            NUMFILES = size(filenames,2);
            filepointer = FULLFILENAMES{CURRFILE};
            open_file(filepointer);
        end
    end

    function cb_metabut(~,~)
        if isempty(IMMETA)
            msgbox('No metadata has been extracted!');
        else
            % Save metadata to .txt file
            disp('Metadata is written to .txt file...');
            metapath = strcat(spath,'\', FNAME,'_metadata.txt');
            fid = fopen(metapath,'w');
            fprintf(fid,'%s\n',IMMETA);
            fclose(fid);
            open(metapath);
        end
    end

    function cb_updatebut(~,~)
        update_fig();
    end

    function cb_newbut(~,~)
        gen_new_fig();
    end

    function cb_scmapdd(hObj,~)
        SCMAPlist = get(hObj,'String');
        selected = get(hObj,'Value');
        SCMAP = SCMAPlist{selected};
    end

    function cb_winszdd(hObj,~)
        winszlist = get(hObj,'String');
        selected = get(hObj,'Value');
        WINSZ = str2double(winszlist{selected});
    end

    function cb_sclimin1(hObj,~)
        CLIMUI(1) = str2double(get(hObj,'String'));
    end

    function cb_sclimin2(hObj,~)
        CLIMUI(2) = str2double(get(hObj,'String'));
    end

    function cb_rangein1(hObj,~)
        PLOTRANGE(1) = str2double(get(hObj,'String'));
    end

    function cb_rangein2(hObj,~)
        PLOTRANGE(2) = str2double(get(hObj,'String'));
    end

    function cb_roilb(hObj,~)
        roilist = get(hObj,'String');
        selected = get(hObj,'Value');
        for iF = 1:size(roiINFO,2),roiINFO(iF).selected = false; end
        for iS = 1:numel(selected)
            try roiINFO(strcmp({roiINFO(:).name}, roilist{selected(iS)})).selected = true; end
        end
    end

    function cb_evsmthin(hObj,~)
        SMTHWIN = str2double(get(hObj,'String'));
    end

    function cb_threshin(hObj,~)
        THRESHOLD = str2double(get(hObj,'String'));
    end

    function cb_peakthreshin(hObj,~)
        PEAKTHRESHOLD = str2double(get(hObj,'String'));
    end


    function cb_deltrcbut(~,~)
        if sum([roiINFO(:).selected]) == 0 || isempty(roiINFO(1).position)
            msgbox('You have to select a ROI first!');
        else
            del = find([roiINFO(:).selected] == true);
            nroi = size(roiINFO,2);
            if isempty(roiINFO(end).position), nroi = nroi-1; end
            if ~isempty(traceINFO(1).figID)
                deltrc = zeros(numel(del),1);
                for iE = 1:numel(del)
                    deltrc(iE) = find([traceINFO(:).roiID] == roiINFO(del(iE)).ID);
                end
                traceINFO = traceINFO(~deltrc);
            end
            if numel(del) == nroi
                roiINFO = roiINFO(1); roiINFO(1).mask = []; roiINFO(1).position = []; roiINFO(1).name = []; roiINFO(1).ID = []; roiINFO(1).selected = []; roiINFO(1).saved = []; roiINFO(1).mode = []; roiINFO(1).PLOTRANGE = [];
            else
                roiINFO(del) = [];
            end
            figures = get(groot,'Children');
            tmpfig = figures(strncmp({figures(:).Name}, 'Line',4));
            roilb = findobj(tmpfig, 'type', 'uicontrol', 'style', 'listbox');
            if isempty(roiINFO(1).name), set(roilb, 'Value', 1); set(roilb, 'string', '');
            else, set(roilb, 'Value',1); set(roilb, 'string', {roiINFO(:).name});
            end
            update_fig();
            [~,closeidx] = find(strncmp({figures(:).Name},'ROI',3)==1);
            if ~isempty(closeidx), for iRoi = 1:numel(closeidx), close(figures(closeidx(iRoi))); end; end
        end
    end

    function cb_showtrcbut(~,~)
        if sum([roiINFO(:).selected]) == 0 || isempty(roiINFO(1).position)
            msgbox('Please select a ROI first!');
        else
            traces_overview();
        end
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
                uiresume();
            else
                figures = get(groot,'Children');
                tmpfig = figures(~strncmp({figures(:).Name}, 'Line',4));
                for iF = 1:numel(tmpfig), close(tmpfig(iF)); end
                filepointer = FULLFILENAMES{CURRFILE};
                open_file(filepointer);
            end
        end
    end

    function cb_nosavengobut(~,~)
        reassure = questdlg('Continue without saving?', '', 'Yes', 'Save structs', 'Return', 'Save structs');
        if ~strcmp(reassure, 'Return')
            if strcmp(reassure, 'Save structs')
                % Save structs anyway
                for iE = 1: size(traceINFO,2), traceINFO(iE).frametime = FTIME; end
                disp(strcat('Save struct files...  @ ', spath, FNAME));
                savedir = strcat(SAVEPARENT, '\', FNAME, '\');
                if ~exist(savedir,'dir'), mkdir(savedir); end
                savepointer = strcat(savedir,FNAME,'_');
                save(strcat(savepointer, '_figure_info', '.mat'), 'figINFO');
                save(strcat(savepointer, '_roi_info', '.mat'), 'roiINFO');
                save(strcat(savepointer, '_traces_info', '.mat'), 'traceINFO');
                save_files(savepointer, figidx, 'all');
            end
            
            % Open next file
            CURRFILE = CURRFILE +1;
            if CURRFILE > NUMFILES
                disp('All files have been handled.');
                uiresume();
            else
                figures = get(groot,'Children');
                tmpfig = figures(~strncmp({figures(:).Name}, 'Line',4));
                for iF = 1:numel(tmpfig), close(tmpfig(iF)); end
                filepointer = FULLFILENAMES{CURRFILE};
                open_file(filepointer);
            end
        end
    end

end