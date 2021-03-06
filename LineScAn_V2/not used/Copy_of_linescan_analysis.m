%% ** Script for line scan analysis **
function main_linescan_analysis()

%% Select file(s)
[filenames, pathnames] = uigetfile({'*.nd2', '*.tif'}, 'Multiselect', 'on');

if isa(filenames,'cell') || isa(filenames,'char')
    %% Initialize variables
    if ~isa(filenames,'cell'), filenames = {filenames}; pathnames = {pathnames}; end
    fullfilenames = fullfile(pathnames, filenames);
    numfiles = size(filenames,2);
    iFile = 1;
    scrsz = get(0, 'Screensize'); scrsz = scrsz(3:4);
    whfig = [800 600];
    whtrcfig = [1000 600];
    positionfig = [100 scrsz(2)-whfig(2)-100 whfig];
    positiontrcfig = [100 scrsz(2)-whfig(2)-100 whtrcfig];
    positionroiselect = [0.05 0.1 0.9 0.85];
    fontsize = 8;
    tfszmpl = 1.3;
    scmapdditems = {'gray', 'hot', 'summer', 'winter', 'parula'};
    scmap = scmapdditems{1};
    winszdditems = {'4', '8', '16', '32'};
    winsz = str2double(winszdditems{1}); % Line average
    plotrange = [1 1000];
    climraw = [0 255];
    climui = [-3 3];
    trcylabel = 'Av. Intensity [a.u.]';
    trcylabeldFoF = 'DeltaF/F';
    trcxlabel = 'time [s]';
    savepointer = '';
    trcplotcol1 = [0.6510 0.0667 0.2784];
    trcplotcol2 = [0.0471 0.3137 0.4588];
    trcalpha = 0.25;
    threshlw = 1;
    threshla = 0.7;
    evplotcol = [0 0.1882 0.1059];
    plotPeaks = true; % otherwise plot threshold crossings
    risetime = inputdlg('Upper limit sensor rise time [s]:','Sensor info', [1 60], {'0.2'});
    risetime = str2double(risetime{1});
    threshold = 2;
    peakthreshold = threshold+0.5;
    win_ms = 0.1; % smoothing for peak detection [s]
        
    while iFile <= numfiles
        % Load file
        if numfiles > 1,filepointer = fullfilenames{iFile};
        else, filepointer = fullfilenames{1};
        end
        [spath, fname, ~] = fileparts(filepointer);
        [tmpstack, imInfo] = load_bf_file(filepointer, true);

        % New parameters
        savepointer = '';
        figcounter = 0;
        roiInfo = struct(); roiInfo(1).name = []; roiInfo(1).mask = []; roiInfo(1).position = []; roiInfo(1).ID = []; roiInfo(1).selected = []; roiInfo(1).saved = []; roiInfo(1).mode = []; roiInfo(1).plotrange = [];
        figInfo = struct('IDs', [], 'name', [], 'plotrange',[], 'cscmap',[], 'csclimits', [], 'avwinsize',[], 'saved',[]);
        traceInfo = struct('figID',[], 'fig_params',[],'roiID',[], 'binned_roi_av',[],'dFoF_roi_av',[], 'timestamp',[], 'save',[], 'currmode',[], 'showtot', []);
        
        % Collapse 3rd dimension
        dims = size(tmpstack);
        composite2D = zeros(dims(1)*dims(3), dims(2));
        catrange = 1:dims(1);
        for frame = 1:dims(3)
            composite2D(catrange,:) = tmpstack(:,:,frame);
            catrange = catrange+dims(1);
        end
        imhw = size(composite2D);
        
        % Frametime from metadata
        immeta =imInfo.metadata;
        timestrings = flipud(immeta(strncmp(immeta, 'timestamp',9)));
        ftimeseries = cellfun(@(x) regexp(x,'= ','split'), timestrings,'UniformOutput', false);
        ftimeseries = cellfun(@(x) str2double(x{2}), ftimeseries, 'UniformOutput', false);
        ftime = mean(diff(cell2mat(ftimeseries)))/dims(1); %[s]

        % GUI & plotting
        close all;
        ini_ctrl_box();

        uiwait;
        iFile = iFile+1;
        close all;
        clear('immeta','timstampidx','roiInfo', 'traceInfo', 'figInfo', 'tmpstack', 'composite2D', 'dims', 'catrange', 'imhw');
    end
end

%% Main analysis functions
    function [binned] = average_linescan(input, wsz) % Bin input matrix
        tmphw = size(input);
        hw = round(wsz/2);
        binned = zeros(tmphw);
        iAv = hw;
        while iAv <= tmphw(2)-hw
            binned(:,iAv-hw+1:iAv+hw) = repmat(mean(input(:,iAv-hw+1:iAv+hw), 2), [1 wsz]);
            iAv = iAv+wsz;
        end
        binned(:,1:hw) = repmat(mean(input(:,1:hw), 2), [1 hw]);
        if mod(tmphw(2),wsz) ~= 0, binned(:,end-(iAv-wsz):end) = repmat(mean(input(:,end-(iAv-wsz):end), 2), [1 hw]);
        else, binned(:,end-hw+1:end) = repmat(mean(input(:,end-hw+1:end), 2), [1 hw]);
        end
    end

    function [deltaFovF] = deltaFovF_linescan(input)
        fix_avF = mean(input,1);
        deltaFovF = input-fix_avF;
        deltaFovF = deltaFovF./fix_avF;
    end

    function [scaledinput] = scale_data(input, newscale)
        newdelta = diff(newscale);
        oldscale = [min(input,[],'all') max(input,[],'all')];
        input = input-oldscale(1); % Shift and set minimum to zero
        oldscale = [0 max(input,[],'all')];
        scaledinput = input.*(100/oldscale(2)); %...express as % of max value
        scaledinput = newscale(1) + scaledinput.*(newdelta/100);
    end

%% GUI functions
    function ini_ctrl_box()
        % Spacing parameteres
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
        
        % Positions
        positionctrl = [scrsz(1)-round(1.5*whctrl(1)) round(scrsz(2)/2-0.5*whctrl(2)) whctrl];
        deltrcbutpos = [hctr-whbutsm(1)*2-hspace(1)*2 vspace(2) whbutsm];
        showtrcbutpos = [hctr-whbutsm(1)-hspace(1) deltrcbutpos(2) whbutsm];
        saveallbutpos = [hctr+hspace(1) deltrcbutpos(2) whbutsm];
        continuebutpos = [hctr+whbutsm(1)+hspace(1)*2 deltrcbutpos(2) whbutsm];
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
        scmaptxtpos = [hspace(2) scmapddpos(2)+scmapddpos(4)-whtxt(2) whtxt];
        winszddpos = [hctr+hspace(2) scmaptxtpos(2)+whtxt(2)+vspace(1) whdd(1) whdd(2)];
        winsztxtpos = [hspace(2) winszddpos(2)+winszddpos(4)-whtxt(2) whtxt];
        rangetxtpos = [hspace(2) winsztxtpos(2)+whtxt(2)+vspace(1) whtxt];
        rangeinpos = [hctr+hspace(2) rangetxtpos(2) whin;...
            hctr+hspace(2)+whin(1)+hspace(1) rangetxtpos(2) whin];
        metabutpos = [hctr-whbutsm(1)/2 rangetxtpos(2)+rangetxtpos(4)+vspace(1) whbutsm];
        nametxtpos = [hspace(2) metabutpos(2)+metabutpos(4)+vspace(2) whtxtbg];
        
        % Main window
        ctrlfig = figure('Name', 'Line Scan Analysis', 'Position', positionctrl, 'toolbar', 'none', 'menu', 'none');
        set(gca, 'Color', 'none'); axis off;
        title('Control Box', 'FontSize', round(fontsize*tfszmpl));
        
        % Add UI controls
        deltrcbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', deltrcbutpos,'string', 'Delete ROI','fontsize', fontsize, 'callback', {@cb_deltrcbut});
        showtrcbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', showtrcbutpos,'string', 'Show traces','fontsize', fontsize, 'callback', {@cb_showtrcbut});
        saveallbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', saveallbutpos,'string', 'Save all','fontsize', fontsize, 'callback', {@cb_saveallbut});
        continuebut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', continuebutpos,'string', 'Next file','fontsize', fontsize, 'callback', {@cb_continuebut});
        
        threshtxt = uicontrol('parent', ctrlfig, 'style', 'text','position', threshtxtpos,'string', 'Threshold [s.d.]: start | peak','fontsize', fontsize);
        threshin = uicontrol('parent', ctrlfig, 'style', 'edit','position', threshinpos,'fontsize', fontsize, 'String', num2str(threshold),'Callback', {@cb_threshin});
        peakin = uicontrol('parent', ctrlfig, 'style', 'edit','position', peakinpos,'fontsize', fontsize, 'String', num2str(peakthreshold),'Callback', {@cb_peakthreshin});
                
        evsmthtxt = uicontrol('parent', ctrlfig, 'style', 'text','position', evsmthtxtpos,'string', 'Smoothing window [s]','fontsize', fontsize);
        evsmthin = uicontrol('parent', ctrlfig, 'style', 'edit','position', evsmthinpos,'fontsize', fontsize, 'String', num2str(win_ms),'Callback', {@cb_evsmthin});
        
        roilbtxt  = uicontrol('parent', ctrlfig, 'style', 'text','position', roilbtxtpos,'string', 'Select ROI for analysis:','fontsize', fontsize);
        roilb = uicontrol('parent', ctrlfig, 'style', 'listbox','position', roilbpos,'fontsize', fontsize, 'string', '', 'max', 20,'Callback', {@cb_roilb});
        
        newbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', newbutpos,'string', 'New figure','fontsize', fontsize, 'callback', {@cb_newbut});
        updatebut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', updatebutpos,'string', 'Update figure', 'fontsize',fontsize, 'callback', {@cb_updatebut});
        sclimtxt = uicontrol('parent', ctrlfig, 'style', 'text','position', sclimtxtpos,'string', 'Set colorscale limits:','fontsize', fontsize);
        sclimin1 = uicontrol('parent', ctrlfig, 'style', 'edit','position', scliminpos(1,:),'fontsize', fontsize, 'String', num2str(climui(1)),'Callback', {@cb_sclimin1});
        sclimin2 = uicontrol('parent', ctrlfig, 'style', 'edit','position', scliminpos(2,:),'fontsize', fontsize, 'String', num2str(climui(2)),'Callback', {@cb_sclimin2});
        
        scmaptxt = uicontrol('parent', ctrlfig, 'style', 'text','position', scmaptxtpos,'string', 'Choose LUT:','fontsize', fontsize);
        scmapdd = uicontrol('parent', ctrlfig, 'style', 'popupmenu','position', scmapddpos,'fontsize', fontsize, 'string', scmapdditems, 'Callback', {@cb_scmapdd});
        
        winsztxt = uicontrol('parent', ctrlfig, 'style', 'text','position', winsztxtpos,'string', 'Select averaging window size:','fontsize', fontsize);
        winszdd = uicontrol('parent', ctrlfig, 'style', 'popupmenu','position', winszddpos,'fontsize', fontsize, 'string', winszdditems, 'Callback', {@cb_winszdd});
        
        rangetxt = uicontrol('parent', ctrlfig, 'style', 'text','position', rangetxtpos,'string', 'Set plot samplepoint range:','fontsize', fontsize);
        rangein1 = uicontrol('parent', ctrlfig, 'style', 'edit','position', rangeinpos(1,:),'fontsize', fontsize, 'String', num2str(plotrange(1)),'Callback', {@cb_rangein1});
        rangein2 = uicontrol('parent', ctrlfig, 'style', 'edit','position', rangeinpos(2,:),'fontsize', fontsize, 'String', num2str(plotrange(2)),'Callback', {@cb_rangein2});
    
        nametxt = uicontrol('parent', ctrlfig, 'style', 'text','position', nametxtpos,'string', fname,'fontsize', fontsize);
        metabut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', metabutpos,'string', 'Metadata','fontsize', fontsize, 'callback', {@cb_metabut});
    end

    function gen_new_fig()
        figcounter = figcounter+1;
        figure('Position', positionfig)
        set(gcf, 'UserData', figcounter);
        set(gcf, 'Name', sprintf('Figure #%i', figcounter));
        figInfo(figcounter).IDs = figcounter;
        figInfo(figcounter).name = sprintf('Figure #%i', figcounter);
        figInfo(figcounter).saved = false;
        update_fig();
    end

    function update_fig()
        figures = get(groot,'Children');
        if numel(figures) == 1
            figcounter = figcounter+1;
            currfig = figure('Position', positionfig);
            set(gcf, 'UserData', figcounter);
            set(gcf, 'Name', sprintf('Figure #%i', figcounter));
            currcnt = figcounter;
            figInfo(currcnt).IDs = figcounter;
            figInfo(currcnt).name = sprintf('Figure #%i', figcounter);
            figInfo(currcnt).saved = false;
        else
            findfig = strncmp({figures(:).Name}, 'Fig',3);
            currfig = figures(find(findfig,1, 'first'));
            currcnt = get(currfig, 'UserData');
        end
        figInfo(currcnt).plotrange = plotrange;
        figInfo(currcnt).avwinsize = winsz;
        figInfo(currcnt).cscmap = scmap;
        figInfo(currcnt).csclimits = climui;
        
        % Parameters
        spcol = 2; sprow = 2;
        hspace = 30;
        vspace = 60;
        whsp = [floor((whfig(1)-3*hspace)/spcol) floor((whfig(2)-3*vspace)/sprow)];
        whsp = whsp./whfig; hspace = hspace/whfig(1); vspace = vspace/whfig(2);
        positionspbl = [hspace vspace whsp];
        positionspbr = [positionspbl(1)+positionspbl(3)+hspace vspace whsp];
        positionspul = [hspace positionspbl(2)+positionspbl(4)+vspace whsp];
        positionspur = [positionspul(1)+positionspul(3)+hspace positionspbl(2)+positionspbl(4)+vspace whsp];
        
        % Call Analysis functions
        if plotrange(2)>imhw(1), plotrange=[1 imhw(1)]; end % Make sure the plot range doesnt exceed sampling timepoints
        rawmtrx = composite2D(plotrange(1):plotrange(2),:);
        averaged = average_linescan(rawmtrx, winsz);
        dFoF = deltaFovF_linescan(averaged);
        
        % Plot        
        set(0,'currentfig',currfig);
        figID = get(currfig,'UserData');
        colormap(scmap); 
        % Raw (scaled to 255)
%         scraw = scale_data(raw, climraw);
        scraw = uint8(rawmtrx);
        if ~isempty(roiInfo(1).ID), scrawroi = mark_rois(scraw,plotrange,'all'); else, scrawroi=scraw; end
        sp1 = subplot('Position', positionspul);
        im1 = imagesc(scrawroi, climraw);
        axis off;
        spname = 'Raw';
        title(spname);
        colorbar;
        set(im1, 'ButtonDownFcn', {@cb_subplot, scraw, figID, false});
        
        % Averaged (scaled to 255)
%         scav = scale_data(averaged, climraw);
        scav = uint8(averaged);
        if ~isempty(roiInfo(1).ID), scavroi = mark_rois(scav,plotrange,'all'); else, scavroi=scav; end
        sp2 = subplot('Position', positionspur, 'parent', currfig);
        im2 = imagesc(scavroi, climraw);
        axis off;
        spname = 'Averaged';
        title(spname);
        colorbar;
        set(im2, 'ButtonDownFcn', {@cb_subplot, scav, figID, false});
        
        % DeltaFovF (scaled to user-defined scale)
        dFoF(dFoF<climui(1)) = climui(1); dFoF(dFoF>climui(2)) = climui(2); 
        dFoF = scale_data(dFoF, climui);
        if ~isempty(roiInfo(1).ID), dFoFroi = mark_rois(dFoF, plotrange,'all'); else, dFoFroi=dFoF; end
        sp3 = subplot('Position', positionspbl, 'parent', currfig);
        im3 = imagesc(dFoFroi, climui);
        axis off;
        spname = 'Delta F over F';
        title(spname);
        colorbar;
        set(im3, 'ButtonDownFcn', {@cb_subplot,dFoF, figID, true});
        
        % Info
        currinfo = {'Information:', '', sprintf('Plotrange: %i - %i', figInfo(figID).plotrange),...
            '', sprintf('Line averaging window size: %i',figInfo(figID).avwinsize),...
            '', sprintf('Colormap: %s', figInfo(figID).cscmap),...
            '', sprintf('Color scale limits: %i - %i', figInfo(figID).csclimits),...
            '', sprintf('Time per frame: %s s', num2str(round(ftime,3)))};
        positionspbr = positionspbr.*[whfig whfig]; positionspbr(2) = round(positionspbr(2)*0.75);
        infotxt = uicontrol('style','text', 'parent', currfig, 'position', positionspbr, 'String', currinfo,'fontsize', fontsize);
    end
    
    function tmpmtrx = mark_rois(tmpmtrx, pr, rois)
        if strcmp(rois,'all')
            roiidx = 1:size(roiInfo,2);
            if isempty(roiInfo(end).name), roiidx= roiidx(1:end-1); end
        else
            roiidx = rois;
        end
        tmpmarkval = max(tmpmtrx,[],'all')*0.3;
        for iRoi = 1:numel(roiidx)
            roimask = roiInfo(roiidx(iRoi)).mask;
            if roiInfo(roiidx(iRoi)).mode == 1
                newmask = false(size(tmpmtrx));
                x1 = find(roimask(1,:),1,'first'); x2 = find(roimask(1,:),1,'last');
                newmask(:,x1:x2) = true;
            else
                % Check if ROI still lies in plotted range
                roipr = roiInfo(roiidx(iRoi)).plotrange;
                if any(roipr ~= pr)
                    if pr(1) <= roipr(1) && pr(2) <= roipr(1), newmask = false(size(roimask));
                    elseif pr(1) <= roipr(1)&& pr(2) <= roipr(2), ycut = 1:pr(2)-roipr(1)+1; ypaste = roipr(1)-pr(1)+1:pr(2)-pr(1)+1;
                    elseif pr(1) <= roipr(1) && pr(2) >= roipr(2), ycut = 1:roipr(2)-roipr(1)+1; ypaste = roipr(1)-pr(1)+1:roipr(2)-pr(1)+1;
                    elseif pr(1) >= roipr(1) && pr(2) <= roipr(2), ycut = pr(1)-roipr(1)+1:pr(2)-roipr(1)+1; ypaste = pr(1)-pr(1)+1:pr(2)-pr(1)+1;
                    elseif pr(1) >= roipr(1) && pr(2) >= roipr(2), ycut = pr(1)-roipr(1)+1:roipr(2)-roipr(1)+1; ypaste = pr(1)-pr(1)+1:roipr(2)-pr(1)+1;
                    end
                    newmask = false(size(tmpmtrx));
                    newmask(ypaste,:) = roimask(ycut,:);
                else
                    newmask = roimask;
                end
            end
            tmpmtrx(newmask) = tmpmtrx(newmask)+tmpmarkval;
        end
    end

    function create_rois(tmpmtrx, figID, cscui)
        fp = find([figInfo(:).IDs] == figID);
        if isempty(roiInfo(1).ID), roitot = 0; else, roitot = size(roiInfo,2)-1; end
        pr = figInfo(fp).plotrange;
        % Spacing parameteres
        vspace = [20 20]; % [small bigg bigger]
        hspace = [10 40 25]; % [ctr-aligned rest]
        whbut = [80 20];
        whtxt = [100 15];
        
        % Positions
        lineroipos = [hspace(2) vspace(1) whbut];
        rectroipos = [lineroipos(1)+lineroipos(3)+hspace(1) vspace(1) whbut];
        closeroipos = [whfig(1)-whbut(1)-hspace(1)*4 vspace(1) whbut];
        delroipos = [closeroipos(1)-hspace(1)-whbut(1) vspace(1) whbut];
        addroipos = [delroipos(1)-hspace(1)-whbut(1) vspace(1) whbut];
        roistxtpos = [addroipos(1)-hspace(1)-whtxt(1) vspace(1) whtxt];
        
        % Main window
        roifig = figure('Position', positionfig, 'Name', 'ROI selection');
        set(gca, 'Color', 'none'); axis off;
        roiax = subplot('Position', positionroiselect, 'parent', roifig);
        colormap(figInfo(fp).cscmap);
        if cscui, cscui = figInfo(fp).csclimits; else, cscui = climraw; end
        if ~isempty(roiInfo(1).ID), tmpmtrxroi = mark_rois(tmpmtrx, pr, 'all'); else, tmpmtrxroi=tmpmtrx; end
        imagesc(tmpmtrxroi, cscui);
        axis off;
        
        % Add UI controls
        lineroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', lineroipos,'string', 'Line ROI','fontsize', fontsize, 'callback', {@cb_drawroi,roiax, pr});
        rectroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', rectroipos,'string', 'Rectangle ROI','fontsize', fontsize, 'callback', {@cb_drawroi,roiax, pr});
        addroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', addroipos,'string', 'Add ROI','fontsize', fontsize, 'callback', {@cb_addroi, fp, tmpmtrx, pr, cscui});
        delroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', delroipos,'string', 'Delete all','fontsize', fontsize, 'callback', {@cb_delroi,tmpmtrx, fp, cscui});
        closeroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', closeroipos,'string', 'Close','fontsize', fontsize, 'callback', {@cb_closeroi});
        roistxt  = uicontrol('parent', roifig, 'style', 'text','position', roistxtpos,'string', sprintf('ROIs Total # %i', roitot),'fontsize', fontsize);
    end

    function traces_overview()
        % Image Processing Parameters
        figures = get(groot,'Children');
        findfig = strncmp({figures(:).Name}, 'Fig',3);
        currfig = figures(find(findfig,1, 'first'));
        if ~isempty(currfig)
            figid = currfig.UserData;
        else
            warning('As no figure is open, the parameters of the last created figure are used');
            figid = max([figInfo(:).IDs]);
        end
        [~,figidx] = find([figInfo(:).IDs] == figid);
        pr = figInfo(figidx).plotrange;
        trcxrange = [pr(1)*ftime pr(2)*ftime]; %[s]
        wsz = figInfo(figidx).avwinsize;
        cscm = figInfo(figidx).cscmap;
        cscl = figInfo(figidx).csclimits;
        
        % ROI Parameters
        [~,roiidx] = find([roiInfo(:).selected] == 1);
        numrois = numel(roiidx);
        
        traceidx = [];
        for iRoi = 1:numrois
            trcexists = false;
            existidx1 = [traceInfo(:).figID] == figid;
            existidx2 = [traceInfo(:).roiID] == roiInfo(roiidx(iRoi)).ID;
            if any(existidx1 & existidx2)
                existidx = find(existidx1 & existidx2);
                iEx = 1;
                while iEx <= numel(existidx)
                    traceInfo(existidx(iEx)).save = 0;
                    if all(traceInfo(existidx(iEx)).fig_params{1,1} == pr) &&...
                            traceInfo(existidx(iEx)).fig_params{2,1} == wsz &&...
                            traceInfo(existidx(iEx)).fig_params{3,1} == win_ms &&...
                            traceInfo(existidx(iEx)).fig_params{4,1} == threshold % check plotrange & binning
                        trcexists = true;
                        traceidx = [traceidx existidx(iEx)];
                        traceInfo(existidx(iEx)).events = get_trc_params(traceInfo(existidx(iEx)).binned_roi_av{1}, [], traceInfo(existidx(iEx)).events.average, traceInfo(existidx(iEx)).events.sd);
                        if roiInfo(roiidx(iRoi)).mode == 1
                            traceInfo(existidx(iEx)).tot_events = get_trc_params(traceInfo(existidx(iEx)).tot_binned_roi_av{1}, [], traceInfo(existidx(iEx)).tot_events.average, traceInfo(existidx(iEx)).tot_events.sd);
                        end
                        break;
                    end
                    iEx = iEx+1;
                end
            end
            if ~trcexists
                traceidx = create_trc(figidx, roiidx, iRoi, traceidx, pr, wsz);
            end
        end
        numrois = numel(traceidx);
        
        if ~isempty(traceidx)
            % Prepare display
            spcol = 12; sprow = 3;
            hspace = 20;
            vspace = 30;
            whbut = [100 25];
            whsp = [floor((whtrcfig(1)-7*hspace)/spcol) floor((whtrcfig(2)-(sprow-1+4)*vspace-whbut(2))/sprow)];
            whrbg = [whsp(1)*2 round(whsp(2)/2)];
            whrb = [round(0.8*whrbg(1)) round(whrbg(2)/3)];
            whan = [120 50];
            savebutpos = [round(whtrcfig(1)/2)-round(hspace/2)-whbut(1) round(vspace/2) whbut];
            closebutpos = [round(whtrcfig(1)/2)+round(hspace/2) round(vspace/2) whbut];
            rbbotpos = [0.15*whrbg(1) round(whrb(2)*0.5) whrb];
            rbtoppos = [0.15*whrbg(1) round(whrb(2)*1.5) whrb];
            whsp = whsp./whtrcfig; whrbg = whrbg./whtrcfig; hspace = hspace/whtrcfig(1); vspace = vspace/whtrcfig(2); whbut=whbut./whtrcfig; whrb=whrb./whtrcfig; whan=whan./whtrcfig; 
            spposleft = zeros(sprow,4);
            spposmid = zeros(sprow,4);
            spposright = zeros(sprow,4);
            bgpos = zeros(sprow,4);
            rbsavetrcpos = zeros(sprow,4);
            rbshowallpos = zeros(sprow,4);
            ypos = vspace*2+whbut(2);
            for iR = 1:sprow
                spposleft(iR,:) = [hspace ypos+vspace whsp(1)*2 whsp(2)];
                spposmid(iR,:) = [spposleft(iR,1)+spposleft(iR,3)+hspace*2 ypos+vspace whsp(1)*6 whsp(2)];
                spposright(iR,:) = [spposmid(iR,1)+spposmid(iR,3)+hspace ypos+vspace whsp(1)*2 whsp(2)];
                bgpos(iR,:) = [spposright(iR,1)+spposright(iR,3)+hspace ypos+vspace+whsp(2)/2 whrbg];
                rbsavetrcpos(iR,:) = [spposright(iR,1)+spposright(iR,3)+hspace ypos+vspace whrb(1)*1.25 whrb(2)];
                rbshowallpos(iR,:) = [rbsavetrcpos(iR,1) rbsavetrcpos(iR,2)+whrb(2) whrb(1)*1.25 whrb(2)];
                ypos = ypos+vspace+whsp(2);
            end
            spposleft = flipud(spposleft); spposmid = flipud(spposmid); spposright = flipud(spposright); bgpos = flipud(bgpos); rbsavetrcpos = flipud(rbsavetrcpos); rbshowallpos = flipud(rbshowallpos);
            figure('Position', positiontrcfig, 'Name', 'Traces 1');
            savebut = uicontrol('parent', gcf, 'style', 'pushbutton', 'position', savebutpos,'string', 'Save selected','fontsize', fontsize, 'callback', {@cb_savetrcbut, figidx});
            closebut =  uicontrol('parent', gcf, 'style', 'pushbutton', 'position', closebutpos,'string', 'Close (no saving)','fontsize', fontsize, 'callback', {@cb_closetrcbut});
            iFig = 1; iSP = 1;
            for iRoi = 1:numrois
                sp1 = subplot('Position', spposleft(iSP,:),'Parent', gcf);
                colormap(cscm);
                imagesc(traceInfo(traceidx(iRoi)).plotmarked{1}, climraw);
                axis off;
                sp2 = subplot('Position', spposmid(iSP,:),'Parent', gcf);
                t1 = plot(traceInfo(traceidx(iRoi)).timestamp{1}, traceInfo(traceidx(iRoi)).binned_roi_av{1},'linewidth', 1, 'color', trcplotcol1); t1.Color(4) = trcalpha;
                hold on, t2 = plot(traceInfo(traceidx(iRoi)).timestamp{1}, traceInfo(traceidx(iRoi)).smoothed{1},'linewidth', 2, 'color', 'black');
                ty = yline(traceInfo(traceidx(iRoi)).events.threshold, 'color', [0 0 0], 'linewidth', threshlw); ty.Color(4) = threshla;
                ty = yline(traceInfo(traceidx(iRoi)).events.peakthreshold, 'color', [0 0 0], 'linewidth', threshlw/2); ty.Color(4) = threshla;
                if plotPeaks, evidx = traceInfo(traceidx(iRoi)).timestamp{1}(traceInfo(traceidx(iRoi)).events.peaks);
                else, evidx = traceInfo(traceidx(iRoi)).timestamp{1}(traceInfo(traceidx(iRoi)).events.crossings);
                end
                for iE = 1:numel(evidx), hold on; xline(evidx(iE), 'Linewidth',2, 'Color', evplotcol); end
                annotxt = {sprintf('Average | SD: %s  |  %s', num2str(traceInfo(traceidx(iRoi)).events.average,3), num2str(traceInfo(traceidx(iRoi)).events.sd,3)),...
                    sprintf('Threshold: %s', num2str(traceInfo(traceidx(iRoi)).events.threshold,3)),...
                    sprintf('Threshold peak: %s', num2str(traceInfo(traceidx(iRoi)).events.peakthreshold,3)),...
                    sprintf('Event type: %s', traceInfo(traceidx(iRoi)).events.eventtype),...
                    sprintf('Av. Inter-event-interval: %s s', num2str(traceInfo(traceidx(iRoi)).events.aviei,3)), ...
                    sprintf('CV IEI: %s', num2str(traceInfo(traceidx(iRoi)).events.cviei,3)),...
                    sprintf('Av. Amplitude: %s', num2str(traceInfo(traceidx(iRoi)).events.avamp,3)),...
                    sprintf('Av. Eventrate: %s Hz', num2str(traceInfo(traceidx(iRoi)).events.eventrate,3))};
                hold off;
                ylabel(trcylabel); set(gca, 'fontsize', fontsize);
                xlabel(trcxlabel); xlim(trcxrange);
                sp3 =  subplot('Position', spposright(iSP,:),'Parent', gcf); axis off;
                evinfo = text(0.02,0.5, annotxt, 'fontsize', fontsize, 'backgroundcolor', [1 1 1], 'edgecolor', [0 0 0]);
                bg = uibuttongroup('parent', gcf, 'visible', 'off', 'position', bgpos(iSP,:), 'SelectionChangedFcn', {@cb_traceswitchrb, sp2, traceidx(iRoi), trcxrange, evinfo});
                if roiInfo(roiidx(iRoi)).mode == 1, showalltrcrb =  uicontrol('parent', gcf,'style', 'radiobutton',  'unit', 'normal','position', rbshowallpos(iSP,:),'string', '  Show full trace', 'Value', 0 ,'fontsize', fontsize, 'callback', {@cb_wholetrcrb, sp2, traceidx(iRoi), trcxrange, evinfo}); end
                savetrcrb = uicontrol('parent', gcf,'style', 'radiobutton', 'unit', 'normal','position', rbsavetrcpos(iSP,:),'string', '  Save', 'Value', 1 ,'fontsize', fontsize, 'callback', {@cb_savetrcrb, traceidx(iRoi)});
                botrb = uicontrol('parent', bg,'style', 'radiobutton', 'position', rbbotpos,'string', '  DeltaF/F', 'Value', 0 ,'fontsize', fontsize, 'handlevisibility', 'off');
                toprb = uicontrol('parent', bg, 'style', 'radiobutton', 'position', rbtoppos,'string', '  Average', 'Value', 1, 'fontsize', fontsize, 'handlevisibility', 'off');
                bg.Visible = 'on';
                if mod(iRoi,sprow) == 0 && iRoi < numrois
                    iFig = iFig+1; iSP = 1;
                    figure('Position', positiontrcfig, 'Name', sprintf('Traces %i',iFig));
                    savebut = uicontrol('parent', gcf, 'style', 'pushbutton', 'position', savebutpos,'string', 'Save selected','fontsize', fontsize, 'callback', {@cb_savetrcbut, figidx});
                    closebut =  uicontrol('parent', gcf, 'style', 'pushbutton', 'position', closebutpos,'string', 'Close (no saving)','fontsize', fontsize, 'callback', {@cb_closetrcbut});
                else
                    iSP = iSP+1;
                end
            end
        end
    end

    function [traceidx] = create_trc(figidx, roiidx, iRoi, traceidx, pr, wsz)
        if isempty(traceInfo(1).figID), tmpidx = 1; else, tmpidx = size(traceInfo,2)+1; end
        traceInfo(tmpidx).figID = figInfo(figidx).IDs;
        traceInfo(tmpidx).roiID = roiInfo(roiidx(iRoi)).ID;
        traceInfo(tmpidx).fig_params{1,1} = pr;
        traceInfo(tmpidx).fig_params{2,1} = wsz;
        traceInfo(tmpidx).fig_params{3,1} = win_ms;
        traceInfo(tmpidx).fig_params{4,1} = threshold;
        traceInfo(tmpidx).save = true;
        traceInfo(tmpidx).currmode = '  Average';
        traceInfo(tmpidx).showtot = 0;
        % Call Analysis functions
        averaged = average_linescan(composite2D, wsz);
        markedroi = uint8(averaged(pr(1):pr(2),:));
        scmarkedroi = mark_rois(markedroi,pr, roiidx(iRoi));
%         scmarkedroi = scale_data(markedroi, climraw);
        [tmpvals,tmpdFoF, yrange, allvals, alldFoF] = calc_roi_av_trace(roiidx(iRoi), averaged, pr);
        if isempty(yrange)
            traceInfo(tmpidx).binned_roi_av = {NaN};
            traceInfo(tmpidx).dFoF_roi_av = {NaN};
            traceInfo(tmpidx).timestamp = {NaN};
            traceInfo(tmpidx).plotmarked = {NaN};
            msgbox(strcat('ROI # ', num2str(iRoi), ' lies outside the plotted range!'));
        else
            traceidx = [traceidx tmpidx];
            traceInfo(tmpidx).binned_roi_av = {tmpvals};
            traceInfo(tmpidx).dFoF_roi_av = {tmpdFoF};
            traceInfo(tmpidx).plotmarked = {scmarkedroi};
            traceInfo(tmpidx).timestamp = {((pr(1)+yrange(1)-1)*ftime:ftime:(pr(1)+yrange(2)-1)*ftime)'};
            % Extract trace parameters
            [events, smoothed] = get_trc_params(allvals, [pr(1)+yrange(1)-1 pr(1)+yrange(2)-1], [], []);
            traceInfo(tmpidx).events = events;
            traceInfo(tmpidx).smoothed = {smoothed};
        end
        if roiInfo(roiidx(iRoi)).mode == 1
            traceInfo(tmpidx).tot_binned_roi_av = {allvals};
            traceInfo(tmpidx).tot_dFoF_roi_av = {alldFoF};
            traceInfo(tmpidx).tot_timestamp = {(1*ftime:ftime:imhw(1)*ftime)'};
            % Extract trace parameters
            [events, smoothed] = get_trc_params(allvals,[],[],[]);
            traceInfo(tmpidx).tot_events = events;
            traceInfo(tmpidx).tot_smoothed = {smoothed};
        end
    end

    function [roivals,roidFoF, yrange, avvals, alldFoF] = calc_roi_av_trace(roiidx, vals, pr)
        roimask = roiInfo(roiidx).mask;
        roipr = roiInfo(roiidx).plotrange;
        % Calculate whole trace average
        collapsey = sum(roimask,1);
        xmin = find(collapsey,1,'first'); xmax = find(collapsey,1,'last');
        avvals = mean(vals(:,xmin:xmax),2);
        av = mean(avvals);
        
        skip = 0;
        if any(roipr ~= pr)
            if pr(1) <= roipr(1) && pr(2) <= roipr(1), yrange = 1;
            elseif pr(1) <= roipr(1) && pr(2) <= roipr(2), ycut = 1:pr(2)-roipr(1)+1; ypaste = roipr(1)-pr(1)+1:pr(2)-pr(1)+1;
            elseif pr(1) <= roipr(1) && pr(2) >= roipr(2), ycut = 1:roipr(2)-roipr(1)+1; ypaste = roipr(1)-pr(1)+1:roipr(2)-pr(1)+1;
            elseif pr(1) >= roipr(1) && pr(2) <= roipr(2), ycut = pr(1)-roipr(1)+1:pr(2)-roipr(1)+1; ypaste = pr(1)-pr(1)+1:pr(2)-pr(1)+1;
            elseif pr(1) >= roipr(1) && pr(2) >= roipr(2), ycut = pr(1)-roipr(1)+1:roipr(2)-roipr(1)+1; ypaste = pr(1)-pr(1)+1:roipr(2)-pr(1)+1;
            end
            newmask = false(diff(pr)+1, size(vals,2));
            newmask(ypaste,:) = roimask(ycut,:);
        else
            newmask = roimask;
        end
        
        if ~skip
            % Calculate ROI average
            collapsex = sum(newmask,2);
            ymin = find(collapsex,1,'first'); ymax = find(collapsex,1,'last');
            % Av values
            roivals = avvals(pr(1)+ymin-1:pr(1)+ymax-1,1);
            % dFoF
            roidFoF = (roivals-av)/av;
            alldFoF = (avvals-av)/av;
            yrange = [ymin ymax];
        else
            roivals = []; roidFoF = []; alldFoF = []; yrange = [];
        end
    end

    function [] = update_trc(sp2, trace, xrange, mode, showall, evinfo)
        subplot(sp2);
        if strcmp(mode, '  Average') && ~showall
            t1 = plot(traceInfo(trace).timestamp{1}, traceInfo(trace).binned_roi_av{1},'linewidth', 1, 'color', trcplotcol1); t1.Color(4) = trcalpha;
            hold on, t2 = plot(traceInfo(trace).timestamp{1}, traceInfo(trace).smoothed{1},'linewidth', 2, 'color', 'black');
            ty = yline(traceInfo(trace).events.threshold, 'color', [0 0 0],'linewidth', threshlw); ty.Color(4) = threshla;
            ty = yline(traceInfo(trace).events.peakthreshold, 'color', [0 0 0], 'linewidth', threshlw/2); ty.Color(4) = threshla;
            events = traceInfo(trace).events;
            if plotPeaks, evidx = traceInfo(trace).timestamp{1}(events.peaks);
            else, evidx = traceInfo(trace).timestamp{1}(events.crossings);
            end
            for iE = 1:numel(evidx), hold on; xline(evidx(iE),'Linewidth',2, 'Color', evplotcol); end
            hold off;
            ylabel(trcylabel); set(gca, 'fontsize', fontsize);
            xlim(xrange);
        elseif strcmp(mode, '  Average') && showall
            t1 = plot(traceInfo(trace).tot_timestamp{1}, traceInfo(trace).tot_binned_roi_av{1},'linewidth', 1, 'color', trcplotcol1); t1.Color(4) = trcalpha;
            hold on, t2 = plot(traceInfo(trace).tot_timestamp{1}, traceInfo(trace).tot_smoothed{1},'linewidth', 2, 'color', 'black');
            ty = yline(traceInfo(trace).tot_events.threshold, 'color', [0 0 0],'linewidth', threshlw); ty.Color(4) = threshla;
            ty = yline(traceInfo(trace).tot_events.peakthreshold, 'color', [0 0 0], 'linewidth', threshlw/2); ty.Color(4) = threshla;
            events = traceInfo(trace).tot_events;
            if plotPeaks, evidx = traceInfo(trace).tot_timestamp{1}(events.peaks);
            else, evidx = traceInfo(trace).tot_timestamp{1}(events.crossings);
            end
            for iE = 1:numel(evidx), hold on; xline(evidx(iE),'Linewidth',2, 'Color', evplotcol); end
            hold off;
            ylabel(trcylabel); set(gca, 'fontsize', fontsize);
            xlim([1*ftime imhw(1)*ftime]);            
        elseif strcmp(mode, '  DeltaF/F') && ~showall
            t1 = plot(traceInfo(trace).timestamp{1}, traceInfo(trace).dFoF_roi_av{1},'linewidth', 1, 'color', trcplotcol2);
            ylabel(trcylabeldFoF); set(gca, 'fontsize', fontsize);
            xlim(xrange);
            events = traceInfo(trace).events;
        elseif strcmp(mode, '  DeltaF/F') && showall
            plot(traceInfo(trace).tot_timestamp{1}, traceInfo(trace).tot_dFoF_roi_av{1},'linewidth', 1, 'color', trcplotcol2);
            ylabel(trcylabeldFoF); set(gca, 'fontsize', fontsize);
            xlim([1*ftime imhw(1)*ftime]);
            events = traceInfo(trace).tot_events;
        end
        xlabel(trcxlabel);
        annotxt = {sprintf('Average | SD: %s  |  %s', num2str(events.average,3), num2str(events.sd,3)),...
                    sprintf('Threshold: %s', num2str(events.threshold,3)),...
                    sprintf('Threshold peak: %s', num2str(events.peakthreshold,3)),...
                    sprintf('Event type: %s', events.eventtype),...
                    sprintf('Av. Inter-event-interval: %s s', num2str(events.aviei,3)), ...
                    sprintf('CV IEI: %s', num2str(events.cviei,3)),...
                    sprintf('Av. Amplitude: %s', num2str(events.avamp,3)),...
                    sprintf('Av. Eventrate: %s Hz', num2str(events.eventrate,3))};
        set(evinfo, 'string', annotxt);
    end

    function [smoothed] = smooth_data(data, win)
        hw = round(0.5*win/ftime);
        iPos = hw+1;
        smoothed = data;
        while iPos < numel(data)- hw+1
            smoothed(iPos) = mean(data(iPos-hw:iPos+hw));
            iPos = iPos+1;
        end
        smoothed(1:hw) = mean(data(1:hw));
        smoothed(end-hw+1:end) = mean(data(end-hw+1:end));
    end

    function [events, data] = get_trc_params(alldata, range, av, sd)
        if isempty(range), range = 1:numel(alldata); else, range = range(1):range(2); end
        if isempty(av) || isempty(sd)
            % Smooth data
            smthdata = smooth_data(alldata, win_ms);
            data = smthdata(range);
            % Threshold crossing: smoothed data and 1st derivative
            av = mean(smthdata);
            sd = std(smthdata);
        else
            data = smooth_data(alldata(range), win_ms);
        end
        trcthreshold = av+threshold*sd;
        trcpkthreshold = av + peakthreshold*sd;
        if sd/av < 0.05
            aboveT = false(size(data));
            crossing = aboveT; peaks = aboveT;
            crossidx = NaN; peakidx = NaN; peakheight = NaN; ieis = NaN; eventrate = NaN; cviei = NaN; evkey = 'NaN';
        else   
            aboveT = data >= trcthreshold;
            crossing = false(size(aboveT));
            iPos = 2;
            while iPos < numel(data)-2
                if aboveT(iPos-1) == 0 && all(aboveT(iPos:iPos+2)) == 1,crossing(iPos) = true; end
                iPos = iPos+1;
            end
            crossidx = find(crossing == 1);

            % Find peak after crossing
            peaks = false(size(data));
            peakidx = zeros(size(crossidx));
            peakheight = zeros(size(crossidx));
            delidx = [];
            riseThresh = ceil(risetime/ftime);
            for iCross = 1:numel(crossidx)
                iPos = crossidx(iCross);
                test = false(size(aboveT));
                while iPos <= numel(aboveT) && aboveT(iPos) == 1
                    test(iPos) = true;
                    iPos = iPos+1;
                end
                [valP,idxP] = max(data(test));
                peakidx(iCross) = crossidx(iCross)+idxP-1;
                % Calculate delta F based on ROI average
                peakheight(iCross) = valP - av;
                % Delete crossings where peak is too small
                if valP <= trcpkthreshold || (peakidx(iCross)-crossidx(iCross)) > riseThresh
                    delidx = [delidx iCross];
                end
            end
            crossing(crossidx(delidx)) = false;
            crossidx(delidx) = []; peakidx(delidx) = []; peakheight(delidx) = [];
            peaks(peakidx) = true;

            % Frequency & Inter-peak-dist
            if plotPeaks, ieis = diff(peakidx).*ftime; evkey = 'peak';
            else, ieis = diff(crossidx).*ftime; evkey = 'crossing';
            end
            cviei = std(ieis)/mean(ieis);
            eventrate = numel(crossidx)/(numel(data)*ftime); % [Hz]
        end
        events=struct();
        events.average = av;
        events.sd = sd;
        events.threshold = trcthreshold;
        events.peakthreshold = trcpkthreshold;
        events.eventtype = evkey;
        events.eventrate = eventrate;
        events.aviei = mean(ieis);
        events.avamp = mean(peakheight);
        events.cviei = cviei;
        events.suprathreshold = aboveT;
        events.crossings = crossing;
        events.peaks = peaks;
        events.crossidx = crossidx;
        events.peakidx = peakidx;
        events.amps = peakheight;
        events.ieis = ieis;
    end

    function save_files(sp, fi, mode)
        diapos = [positionfig(1)+positionfig(3)/2 positionfig(2)+positionfig(4)/2 200 100];
        diafig = figure('Position',diapos,'Name', 'Saving','toolbar', 'none', 'menu', 'none');
        set(gca, 'Color', 'none'); axis off;
        diatxtinfo = {'Files are saved', 'this may take a while...'};
        diatxt = uicontrol('parent', diafig, 'style', 'text', 'position', [10 5 180 80], 'string', diatxtinfo, 'fontsize', fontsize+1);
        pause(0.5);
        
        % Get parameters
        pr = figInfo(fi).plotrange;
        wsz = figInfo(fi).avwinsize;
        
        if ~figInfo(fi).saved
            answer = questdlg('Do you want to save the full 2D matrix (raw, average, dFoF) as .csv ?', 'Saving', 'Yes','No','Yes');
            if strcmp(answer, 'Yes')
                % Save ROI-unrelated stuff
                averaged = average_linescan(composite2D, wsz);
                dFoF = deltaFovF_linescan(averaged);
                % .csv
                set(diatxt,'string',[diatxtinfo,{'','...create .csv files.'}]);
                pause(0.2);
                writematrix(composite2D, strcat(sp,'raw.csv'));
                writematrix(averaged, strcat(sp,'AV_bin_',num2str(wsz), '.csv'));
                writematrix(dFoF, strcat(sp,'dFoF_bin_',num2str(wsz), '.csv'));
                % .tifs
                set(diatxt,'string',[diatxtinfo,{'','...create .tif files.'}]);
                pause(0.2);
                vals = uint16(composite2D(pr(1):pr(2),:));
                imwrite(vals, strcat(sp,'raw_range_', num2str(pr(1)),'-',num2str(pr(2)),'.tif'), 'tif');
                vals = uint16(averaged(pr(1):pr(2),:));
                imwrite(vals, strcat(sp,'AV_','bin_',num2str(wsz),'_range_', num2str(pr(1)),'-',num2str(pr(2)),'.tif'), 'tif');
                vals = uint16(dFoF(pr(1):pr(2),:));
                imwrite(vals, strcat(sp,'dFovF_','bin_',num2str(wsz),'_range_', num2str(pr(1)),'-',num2str(pr(2)),'.tif'), 'tif');
                figInfo(fi).saved = true;
            end
        else
            answer = 'No';
        end
        
        % Determine ROIs to save
        if strcmp(mode, 'all')
            saverois = 1:size(roiInfo,2);
            if isempty(roiInfo(end).position), saverois = saverois(1:end-1); end
        else
            saverois = find([traceInfo(:).save] == 1);
        end
        
        % Check if data already exist
        traceidx = zeros(size(saverois));
        trcexists = false(size(saverois));
        isline = false(size(saverois));
        for iRoi = 1:numel(saverois)
            figid = figInfo(fi).IDs;
            existidx1 = [traceInfo(:).figID] == figid;
            existidx2 = [traceInfo(:).roiID] == roiInfo(saverois(iRoi)).ID;
            if any(existidx1 & existidx2)
                existidx = find(existidx1 & existidx2);
                iEx = 1;
                while iEx <= numel(existidx)
                    if all(traceInfo(existidx(iEx)).fig_params{1,1} == pr) &&...
                            traceInfo(existidx(iEx)).fig_params{2,1} == wsz &&...
                            traceInfo(existidx(iEx)).fig_params{3,1} == win_ms &&...
                            traceInfo(existidx(iEx)).fig_params{4,1} == threshold % check plotrange & binning
                        trcexists(iRoi) = true;
                        traceidx(iRoi) = existidx(iEx);
                        break;
                    end
                    iEx = iEx+1;
                end
            end
            if roiInfo(iRoi).mode == 1, isline(iRoi) = true; end
        end
        
        if any(~trcexists) && strcmp(answer,'No')
            averaged = average_linescan(composite2D, wsz);
        end
        
        for iRoi = 1:numel(saverois)
            roiid = roiInfo(saverois(iRoi)).ID;
            set(diatxt,'string',[diatxtinfo,{'','...save ROI # ', num2str(roiid)}]);
            pause(0.2);
            if ~trcexists(iRoi)
                if isline(iRoi), tmppr = [1 imhw(1)]; else, tmppr = pr; end
                val = uint8(averaged(tmppr(1):tmppr(2),:));
                markedroi = mark_rois(val,tmppr, saverois(iRoi));
                scmarkedroi = scale_data(markedroi, climraw);
                [avvals,dFoFvals,yrange,allvals,~] = calc_roi_av_trace(saverois(iRoi), averaged, tmppr);
                [evinfo, smoothed] = get_trc_params(allvals,[tmppr(1)+yrange(1)-1 tmppr(1)+yrange(2)-1], [],[]);
                crossings = evinfo.crossings;
                supraT = evinfo.suprathreshold;
                peaks = evinfo.peaks;
                timestamp = ((tmppr(1)+yrange(1)-1)*ftime:ftime:(tmppr(1)+yrange(2)-1)*ftime)';
            else
                if isline(iRoi)
                    avvals = traceInfo(traceidx(iRoi)).tot_binned_roi_av{1};
                    dFoFvals = traceInfo(traceidx(iRoi)).tot_dFoF_roi_av{1};
                    smoothed = traceInfo(traceidx(iRoi)).tot_smoothed{1};
                    crossings = traceInfo(traceidx(iRoi)).tot_events.crossings;
                    supraT = traceInfo(traceidx(iRoi)).tot_events.suprathreshold;
                    peaks = traceInfo(traceidx(iRoi)).tot_events.peaks;
                    timestamp = traceInfo(traceidx(iRoi)).tot_timestamp{1};
                    evinfo =  traceInfo(traceidx(iRoi)).tot_events;
                else
                    avvals =  traceInfo(traceidx(iRoi)).binned_roi_av{1};
                    dFoFvals = traceInfo(traceidx(iRoi)).dFoF_roi_av{1};
                    smoothed = traceInfo(traceidx(iRoi)).smoothed{1};
                    crossings = traceInfo(traceidx(iRoi)).events.crossings;
                    supraT = traceInfo(traceidx(iRoi)).events.suprathreshold;
                    peaks = traceInfo(traceidx(iRoi)).events.peaks;
                    timestamp = traceInfo(traceidx(iRoi)).timestamp{1};
                    evinfo =  traceInfo(traceidx(iRoi)).events;
                end
                scmarkedroi = traceInfo(traceidx(iRoi)).plotmarked{1};
            end
            writeMtrx = [timestamp avvals smoothed dFoFvals supraT crossings peaks];
            imwrite(scmarkedroi, strcat(sp,'AV_ROI_',num2str(roiid),'_binning_', num2str(wsz), '.tif'), 'tif');
            writematrix(writeMtrx, strcat(sp,'ROI_',num2str(roiid),'_values_binning_', num2str(wsz), '_smth_', num2str(win_ms), '_thresh_', num2str(threshold),'_threshpk_', num2str(peakthreshold),'.csv'));
            % Save event info as table
            evTable = table();
            if isempty(evinfo.crossidx)
                evTable.crossing_idx = NaN;
                evTable.peaks_idx = NaN;
                evTable.peak_amps = NaN;
                evTable.ieis = NaN;
            else
                evTable.crossing_idx = evinfo.crossidx;
                evTable.peaks_idx = evinfo.peakidx;
                evTable.peak_amps = evinfo.amps;
                evTable.ieis = [evinfo.ieis; NaN];
            end
            lt = size(evinfo.crossidx,1);
            if lt > 1, filler = repelem({''},lt-1)'; else, filler = []; end
            evTable.threshold = [evinfo.threshold; filler];
            evTable.event_type = [evinfo.eventtype; filler];
            evTable.av_eventrate_Hz = [evinfo.eventrate; filler];
            evTable.av_intereventinterval_s = [evinfo.aviei; filler];
            evTable.av_peak_amp = [evinfo.avamp;filler];
            evTable.cv_iei = [evinfo.cviei; filler];
            writetable(evTable, strcat(sp,'ROI_',num2str(roiid),'_eventinfo_binning_', num2str(wsz), '_smth_', num2str(win_ms), '_thresh_', num2str(threshold),'_threshpk_', num2str(peakthreshold),'.csv'));
        end
        close(diafig);
        okbox = msgbox('Files saved', '', 'modal');
    end

%% Local callback
    function cb_metabut(~,~)
        if isempty(imInfo)
            msgbox('No metadata has been extracted!');
        else
            % Save metadata to .txt file
            disp('Metadata is written to .txt file...');
            metapath = strcat(spath,'\', fname,'_metadata.txt');
            fid = fopen(metapath,'w');
            fprintf(fid,'%s\n',imInfo.metadata{:});
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
        scmaplist = get(hObj,'String');
        selected = get(hObj,'Value');
        scmap = scmaplist{selected};
    end

    function cb_winszdd(hObj,~)
        winszlist = get(hObj,'String');
        selected = get(hObj,'Value');
        winsz = str2double(winszlist{selected});
    end

    function cb_sclimin1(hObj,~)
        climui(1) = str2double(get(hObj,'String'));
    end

    function cb_sclimin2(hObj,~)
        climui(2) = str2double(get(hObj,'String'));
    end

    function cb_rangein1(hObj,~)
        plotrange(1) = str2double(get(hObj,'String'));
    end

    function cb_rangein2(hObj,~)
        plotrange(2) = str2double(get(hObj,'String'));
    end

    function cb_roilb(hObj,~)
        roilist = get(hObj,'String');
        selected = get(hObj,'Value');
        for iF = 1:size(roiInfo,2),roiInfo(iF).selected = false; end
        for iS = 1:numel(selected)
            try roiInfo(strcmp({roiInfo(:).name}, roilist{selected(iS)})).selected = true; end
        end
    end

    function cb_deltrcbut(~,~)
        del = find([roiInfo(:).selected] == true);
        nroi = size(roiInfo,2);
        if isempty(roiInfo(end).position), nroi = nroi-1; end
        if ~isempty(traceInfo(1).figID)
            deltrc = zeros(numel(del),1);
            for iE = 1:numel(del)
                deltrc(iE) = find([traceInfo(:).roiID] == roiInfo(del(iE)).ID);
            end
            traceInfo = traceInfo(~deltrc);
        end
        if numel(del) == nroi
            roiInfo = roiInfo(1); roiInfo(1).mask = []; roiInfo(1).position = []; roiInfo(1).name = []; roiInfo(1).ID = []; roiInfo(1).selected = []; roiInfo(1).saved = []; roiInfo(1).mode = []; roiInfo(1).plotrange = [];
        else
            roiInfo(del) = [];
        end
        figures = get(groot,'Children');
        tmpfig = figures(strncmp({figures(:).Name}, 'Line',4));
        roilb = findobj(tmpfig, 'type', 'uicontrol', 'style', 'listbox');
        if isempty(roiInfo(1).name), set(roilb, 'Value', 1); set(roilb, 'string', '');
        else, set(roilb, 'Value',1); set(roilb, 'string', {roiInfo(:).name});
        end
        update_fig();
        [~,closeidx] = find(strncmp({figures(:).Name},'ROI',3)==1);
        if ~isempty(closeidx), for iRoi = 1:numel(closeidx), close(figures(closeidx(iRoi))); end; end
    end

    function cb_showtrcbut(~,~)
        if sum([roiInfo(:).selected]) == 0
            msgbox('Please select a ROI first!');
        end
        traces_overview();
    end

    function cb_subplot(~,~,tmpmtrx, figID, cscui)
        create_rois(tmpmtrx, figID, cscui);
    end

    function cb_drawroi(hObj, ~, tmpax, pr)
        try
            if isempty(roiInfo(1).ID)
                roiInfo(1).ID = 1;
                roiInfo(1).name = 0; roiInfo(1).mask = 0; roiInfo(1).position = 0; roiInfo(1).selected = true; roiInfo(1).saved = 0; roiInfo(1).mode = 0; roiInfo(1).plotrange = 0;
            end
            roicnt = size(roiInfo,2);
            if strncmp(hObj.String, 'Line',4)
                roiInfo(roicnt).mask = false(imhw);
                roi1 = images.roi.Line(tmpax, 'linewidth', 2, 'InteractionsAllowed','none');
                draw(roi1);
                roi2 = images.roi.Line(tmpax, 'linewidth', 2, 'InteractionsAllowed','none');
                draw(roi2);
                x = sort([round(mean(roi1.Position(:,1))) round(mean(roi2.Position(:,1)))], 'ascend');
                roiInfo(roicnt).mask(:,x(1):x(2)) = true;
                roiInfo(roicnt).position = roi1.Position(1,:);
                roiInfo(roicnt).mode = 1;
                roiInfo(roicnt).plotrange = [1 imhw(1)];
            elseif strncmp(hObj.String, 'Rect',4)
                roi = images.roi.Rectangle(tmpax, 'FaceAlpha', 0.2, 'InteractionsAllowed','none');
                draw(roi);
                roiInfo(roicnt).mask = createMask(roi);
                roiInfo(roicnt).position = roi.Position(1:2);
                roiInfo(roicnt).mode = 2;
                roiInfo(roicnt).plotrange = pr;
            end
        end
    end

    function cb_addroi(~, ~,fp, tmpmtrx, pr, cscui)
        roicnt = size(roiInfo,2);
        if isempty(roicnt), roicnt = 1; end
        if isempty(roiInfo(roicnt).position)
            msgbox('Please draw a ROI first!');
        else
            roiInfo(roicnt).name = sprintf('ROI #%i',roiInfo(roicnt).ID);
            roiInfo(roicnt).saved = 0;
            if roicnt == 1, roiInfo(roicnt).selected = true;
            else, roiInfo(roicnt).selected = false;
            end
            figures = get(groot,'Children');
            tmpfig = figures(strncmp({figures(:).Name}, 'Line',4));
            roilb = findobj(tmpfig, 'type', 'uicontrol', 'style', 'listbox');
            set(roilb, 'string', {roiInfo.name});
            tmpfig = figures(1);
            roistxt = findobj(tmpfig, 'type', 'uicontrol', 'style', 'text');
            set(roistxt, 'String', sprintf('ROIs Total # %i', size(roiInfo,2)));
            tmpmtrx = mark_rois(tmpmtrx, pr, 'all');
            subplot('Position', positionroiselect, 'parent', tmpfig);
            colormap(figInfo(fp).cscmap);
            imagesc(tmpmtrx, cscui);
            axis off;
            for iC = 1:size(roiInfo,2)
                if ~isempty(roiInfo(iC).position)
                    tmppos = [whfig.*positionroiselect(1:2)+...
                        roiInfo(iC).position.*(whfig.*positionroiselect(3:4)./fliplr(size(tmpmtrx))) 15 15];
                    tmppos(2) = whfig(2)-tmppos(2);
                    uicontrol('style','text','parent',tmpfig,'position', tmppos,'string',num2str(iC), 'fontsize', fontsize)
                end
            end
            roiInfo(roicnt+1).ID = max([roiInfo(:).ID])+1;
        end
    end

    function cb_delroi(~, ~,tmpmtrx, fp, cscui)
        roiInfo = roiInfo(1);
        roiInfo(1).mask = []; roiInfo(1).position = []; roiInfo(1).name = []; roiInfo(1).ID = []; roiInfo(1).selected = []; roiInfo(1).mode = []; roiInfo(1).plotrange = [];
        figures = get(groot,'Children');
        tmpfig = figures(strncmp({figures(:).Name}, 'Line',4));
        roilb = findobj(tmpfig, 'type', 'uicontrol', 'style', 'listbox');
        set(roilb, 'string', '');
        set(roilb, 'value', 0);
        tmpfig = figures(1);
        roistxt = findobj(tmpfig, 'type', 'uicontrol', 'style', 'text');
        delete(roistxt(1:end-1));
        set(roistxt(end), 'String', '');
        subplot('Position', positionroiselect, 'parent', tmpfig);
        colormap(figInfo(fp).cscmap);
        imagesc(tmpmtrx, cscui);
        axis off;
    end

    function cb_closeroi(~,~)
        figures = get(groot,'Children');
        close(figures(1));
        update_fig();
    end

    function cb_traceswitchrb(~,evdat, sp2, trace, xrange,evinfo)
        traceInfo(trace).currmode = get(evdat.NewValue, 'String');
        update_trc(sp2, trace, xrange, traceInfo(trace).currmode, traceInfo(trace).showtot, evinfo);
    end

    function cb_evsmthin(hObj,~)
        win_ms = str2double(get(hObj,'String'));
    end

    function cb_threshin(hObj,~)
        threshold = str2double(get(hObj,'String'));
    end

    function cb_peakthreshin(hObj,~)
        peakthreshold = str2double(get(hObj,'String'));
    end

    function cb_saveallbut(~,~)
        figures = get(groot,'Children');
        if numel(figures) == 1
            msgbox('You have to generate a figure first!');
        else
            findfig = strncmp({figures(:).Name}, 'Fig',3);
            currfig = figures(find(findfig,1, 'first'));
            figid = currfig.UserData;
            [~,figidx] = find([figInfo(:).IDs] == figid);
            % Saving path and name
            [spath, fname, ~] = fileparts(filepointer);
            savepath = uigetdir(spath);
            if isa(savepath, 'char')
                savepointer = strcat(savepath,'\', fname,'_');
                % Prepare for saving
                save_files(savepointer, figidx, 'all');
            end
        end
    end

    function cb_continuebut(~,~)
        for iE = 1: size(traceInfo,2), traceInfo(iE).frametime = ftime; end
        disp(strcat('Save struct files...  @ ', spath, fname));
        save(strcat(spath, fname, '_figure_info', '.mat'), 'figInfo');
        save(strcat(spath, fname, '_roi_info', '.mat'), 'roiInfo');
        save(strcat(spath, fname, '_traces_info', '.mat'), 'traceInfo');
        uiresume;
    end

    function cb_savetrcrb(hObj, ~, trace)
        traceInfo(trace).save = hObj.Value;
    end

    function cb_wholetrcrb(hObj, ~, sp2, trace, xrange, evinfo)
        traceInfo(trace).showtot = logical(hObj.Value);
        update_trc(sp2, trace, xrange, traceInfo(trace).currmode, traceInfo(trace).showtot, evinfo);
    end

    function cb_savetrcbut(~,~,figidx)
        if isempty(savepointer)
            % Saving path and name
            [spath, fname, ~] = fileparts(filepointer);
            savepath = uigetdir(spath);
            if isa(savepath, 'char'), savepointer = strcat(savepath,'\', fname,'_'); end
        end
        if ~isempty(savepointer), close gcf; save_files(savepointer, figidx, 'selected'); end
    end

    function cb_closetrcbut(~,~)
        for iE = 1: size(traceInfo,2), traceInfo(iE).save = false; end
        close gcf;
    end

end

