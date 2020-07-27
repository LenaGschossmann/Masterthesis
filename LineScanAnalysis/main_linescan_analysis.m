%% ** Script for line scan analysis **
function main_linescan_analysis()

%% Select file(s)
[filenames, pathnames] = uigetfile({'*.nd2';'*.tif'}, 'Multiselect', 'on');

if isa(filenames,'cell') || isa(filenames,'char')
    %% Initialize variables
    if ~isa(filenames,'cell'), filenames = {filenames}; pathnames = {pathnames}; end
    fullfilenames = fullfile(pathnames, filenames);
    numfiles = size(filenames,2);
    iFile = 1;
    scrsz = get(0, 'Screensize'); scrsz = scrsz(3:4);
    whfig = [800 600];
    positionfig = [100 scrsz(2)-whfig(2)-100 whfig];
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
    trcxlabel = 'time [ms]';
    savepointer = '';

    while iFile <= numfiles
        % Load file
        if numfiles > 1,filepointer = fullfilenames{iFile};
        else, filepointer = fullfilenames{1};
        end
        [spath, fname, ~] = fileparts(filepointer);
        [tmpstack, imInfo] = load_bf_file(filepointer, true);

        % New parameters
        immeta =imInfo.metadata;
        timstampidx = strncmp(immeta, 'timestamp',9);
        timemeta = immeta(find(timstampidx, 2, 'first'));
        t1 = regexp(timemeta{2},'= ','split'); t1 = str2double(t1{2});
        t2 = regexp(timemeta{1},'= ','split'); t2 = str2double(t2{2});
        ftime = abs(t2-t1);
        figcounter = 0;
        roiInfo = struct(); roiInfo(1).name = []; roiInfo(1).mask = []; roiInfo(1).position = []; roiInfo(1).ID = []; roiInfo(1).selected = []; roiInfo(1).saved = []; roiInfo(1).mode = []; roiInfo(1).plotrange = [];
        figInfo = struct('IDs', [], 'name', [], 'plotrange',[], 'cscmap',[], 'csclimits', [], 'avwinsize',[], 'saved',[]);
        traceInfo = struct('figID',[], 'fig_params',[],'roiID',[], 'binned_roi_av',[],'dFoF_roi_av',[], 'timestamp',[], 'save',[]);

        % Collapse 3rd dimension
        dims = size(tmpstack);
        composite2D = zeros(dims(1)*dims(3), dims(2));
        catrange = 1:dims(1);
        for frame = 1:dims(3)
            composite2D(catrange,:) = tmpstack(:,:,frame);
            catrange = catrange+dims(1);
        end
        imhw = size(composite2D);

        % GUI & plotting
        close all;
        ini_ctrl_box();

        uiwait;
        iFile = iFile+1;
        close all;
    end
end

%% Main analysis functions
    function [binned] = average_linescan(input, wsz)
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
        whctrl = [400 380];
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
        roilbpos = [hctr+hspace(2) deltrcbutpos(2)+whbutbg(2)+vspace(1) whdd(1) whdd(2)*5];
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
        deltrcbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', deltrcbutpos,'string', 'Delete ROIs','fontsize', fontsize, 'callback', {@cb_deltrcbut});
        showtrcbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', showtrcbutpos,'string', 'Show traces','fontsize', fontsize, 'callback', {@cb_showtrcbut});
        saveallbut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', saveallbutpos,'string', 'Save traces','fontsize', fontsize, 'callback', {@cb_saveallbut});
        continuebut = uicontrol('parent', ctrlfig, 'style', 'pushbutton', 'position', continuebutpos,'string', 'Next file','fontsize', fontsize, 'callback', {@cb_continuebut});
        
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
        raw = composite2D(plotrange(1):plotrange(2),:);
        averaged = average_linescan(raw, winsz);
        dFoF = deltaFovF_linescan(averaged);
        
        % Plot        
        set(0,'currentfig',currfig);
        figID = get(currfig,'UserData');
        colormap(scmap); 
        % Raw (scaled to 255)
        scraw = scale_data(raw, climraw);
        if ~isempty(roiInfo(1).ID), scrawroi = mark_rois(scraw,plotrange,'all'); else, scrawroi=scraw; end
        sp1 = subplot('Position', positionspul);
        im1 = imagesc(scrawroi, climraw);
        axis off;
        spname = 'Raw';
        title(spname);
        colorbar;
        set(im1, 'ButtonDownFcn', {@cb_subplot, scraw, figID, false});
        
        % Averaged (scaled to 255)
        scav = scale_data(averaged, climraw);
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
            '', strcat('Time per frame: ', num2str(round(ftime)), ' ms')};
        positionspbr = positionspbr.*[whfig whfig]; positionspbr(2) = round(positionspbr(2)*0.75);
        infotxt = uicontrol('style','text', 'parent', currfig, 'position', positionspbr, 'String', currinfo,'fontsize', fontsize);
    end
    
    function tmpmtrx = mark_rois(tmpmtrx, pr, rois)
        if strcmp(rois,'all')
            roiidx = 1:size(roiInfo,2);
            if isempty(roiInfo(end).position), roiidx= roiidx(1:end-1); end
        else
            roiidx = rois;
        end
        tmpmarkval = max(tmpmtrx,[],'all')*0.3;
        for iRoi = 1:numel(roiidx)
            % Check if ROI still lies in plotted range
            roipr = roiInfo(roiidx(iRoi)).plotrange;
            roimask = roiInfo(roiidx(iRoi)).mask;
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
        polyroipos = [rectroipos(1)+rectroipos(3)+hspace(1) vspace(1) whbut];
        roistxtpos = [polyroipos(1)+polyroipos(3)+hspace(3)*2 vspace(1) whtxt];
        addroipos = [roistxtpos(1)+roistxtpos(3)+hspace(1) vspace(1) whbut];
        delroipos = [addroipos(1)+addroipos(3)+hspace(1) vspace(1) whbut];
        closeroipos = [delroipos(1)+delroipos(3)+hspace(1) vspace(1) whbut];
        
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
        lineroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', lineroipos,'string', 'Line ROI','fontsize', fontsize, 'callback', {@cb_drawroi,roiax, tmpmtrx, pr});
        rectroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', rectroipos,'string', 'Rectangle ROI','fontsize', fontsize, 'callback', {@cb_drawroi,roiax, tmpmtrx, pr});
        polyroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', polyroipos,'string', 'Polyline ROI','fontsize', fontsize, 'callback', {@cb_drawroi,roiax, tmpmtrx, pr});
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
        figid = currfig.UserData;
        [~,figidx] = find([figInfo(:).IDs] == figid);
        pr = figInfo(figidx).plotrange;
        trcxrange = [pr(1)*ftime pr(2)*ftime];
        wsz = figInfo(figidx).avwinsize;
        cscm = figInfo(figidx).cscmap;
        cscl = figInfo(figidx).csclimits;
        
        % ROI Parameters
        [~,roiidx] = find([roiInfo(:).selected] == 1);
        numrois = numel(roiidx);
        
        % Call Analysis functions
        averaged = average_linescan(composite2D(pr(1):pr(2),:), wsz);
        dFoF = deltaFovF_linescan(averaged);
        
        scav = scale_data(averaged, climraw);
        dFoF(dFoF<cscl(1)) = cscl(1); dFoF(dFoF>cscl(2)) = cscl(2); 
        dFoF = scale_data(dFoF, cscl);
        
        traceidx = zeros(numrois,1);
        for iRoi = 1:numrois
            trcexists = false;
            existidx1 = [traceInfo(:).figID] == figid;
            existidx2 = [traceInfo(:).roiID] == roiInfo(roiidx(iRoi)).ID;
            if any(existidx1 & existidx2)
                existidx = find(existidx1 & existidx2);
                iEx = 1;
                while iEx <= numel(existidx)
                    if all(traceInfo(existidx(iEx)).fig_params{1,1} == pr) &&...
                            traceInfo(existidx(iEx)).fig_params{2,1} == wsz % check plotrange & binning
                        trcexists = true;
                        traceidx(iRoi) = existidx(iEx);
                        break;
                    end
                    iEx = iEx+1;
                end
            end
            if ~trcexists
                if isempty(traceInfo(1).figID), tmpidx = 1; else, tmpidx = size(traceInfo,2)+1; end
                traceidx(iRoi) = tmpidx;
                traceInfo(tmpidx).figID = figInfo(figidx).IDs;
                traceInfo(tmpidx).roiID = roiInfo(roiidx(iRoi)).ID;
                traceInfo(tmpidx).fig_params{1,1} = pr;
                traceInfo(tmpidx).fig_params{2,1} = wsz;
                traceInfo(tmpidx).save = true;
                [tmpvals,tmpdFoF, yrange] = calc_roi_av_trace(roiidx(iRoi), averaged, pr, 'ROI');
                traceInfo(tmpidx).binned_roi_av = {tmpvals};
                traceInfo(tmpidx).dFoF_roi_av = {tmpdFoF};
                traceInfo(tmpidx).timestamp = {(pr(1)+yrange(1)-1)*ftime:ftime:(pr(1)+yrange(2)-1)*ftime};
            end
        end
        
        % Prepare display
        spcol = 8; sprow = 3;
        hspace = 20;
        vspace = 30;
        whbut = [100 25];
        whsp = [floor((whfig(1)-5*hspace)/spcol) floor((whfig(2)-(sprow-1+4)*vspace-whbut(2))/sprow)];
        whrbg = [whsp(1) round(whsp(2)/2)];
        whrb = [round(0.8*whrbg(1)) round(whrbg(2)/3)];
        savebutpos = [round(whfig(1)/2)-round(hspace/2)-whbut(1) round(vspace/2) whbut];
        closebutpos = [round(whfig(1)/2)+round(hspace/2) round(vspace/2) whbut];
        rbbotpos = [0.15*whrbg(1) round(whrb(2)*0.5) whrb];
        rbtoppos = [0.15*whrbg(1) round(whrb(2)*1.5) whrb];
        whsp = whsp./whfig; whrbg = whrbg./whfig; hspace = hspace/whfig(1); vspace = vspace/whfig(2); whbut=whbut./whfig;
        spposleft = zeros(sprow,4);
        spposright = zeros(sprow,4);
        bgpos = zeros(sprow,4);
        rbsavetrcpos = zeros(sprow,4);
        ypos = vspace*2+whbut(2);
        for iR = 1:sprow
            spposleft(iR,:) = [hspace ypos+vspace whsp(1)*2 whsp(2)];
            spposright(iR,:) = [spposleft(iR,1)+spposleft(iR,3)+hspace*2 ypos+vspace whsp(1)*5 whsp(2)];
            bgpos(iR,:) = [spposright(iR,1)+spposright(iR,3)+hspace ypos+vspace+whsp(2)/2 whrbg];
            rbsavetrcpos(iR,:) = [(spposright(iR,1)+spposright(iR,3)+hspace+0.15*whrbg(1))*whfig(1) (ypos+vspace+whsp(2)/6)*whfig(2) whrb];
            ypos = ypos+vspace+whsp(2);
        end
        spposleft = flipud(spposleft); spposright = flipud(spposright); bgpos = flipud(bgpos); rbsavetrcpos = round(flipud(rbsavetrcpos));
        figure('Position', positionfig, 'Name', 'Traces 1');
        savebut = uicontrol('parent', gcf, 'style', 'pushbutton', 'position', savebutpos,'string', 'Save selected','fontsize', fontsize, 'callback', {@cb_savetrcbut, figidx});
        closebut =  uicontrol('parent', gcf, 'style', 'pushbutton', 'position', closebutpos,'string', 'Close (no saving)','fontsize', fontsize, 'callback', {@cb_closetrcbut});
        iFig = 1; iSP = 1;
        for iRoi = 1:numrois
            tracesp1 = subplot('Position', spposleft(iSP,:),'Parent', gcf);
            dFoFroi = mark_rois(dFoF,pr, roiidx(iRoi));
            imagesc(dFoFroi, cscl);
            axis off;
            tracesp2 = subplot('Position', spposright(iSP,:),'Parent', gcf);
            plot(traceInfo(traceidx(iRoi)).timestamp{1}, traceInfo(traceidx(iRoi)).binned_roi_av{1}, 'linewidth', 1, 'color', 'black');
            ylabel(trcylabel); set(gca, 'fontsize', fontsize);
            xlabel(trcxlabel); xlim(trcxrange);
            bg = uibuttongroup('parent', gcf, 'visible', 'off', 'position', bgpos(iSP,:), 'SelectionChangedFcn', {@cb_traceswitchrb, traceidx, iRoi, tracesp2, trcxrange});
            savetrcrb =uicontrol('parent', gcf,'style', 'radiobutton', 'position', rbsavetrcpos(iSP,:),'string', '  Save', 'Value', 1 ,'fontsize', fontsize, 'callback', {@cb_savetrcrb, traceidx, iRoi});
            botrb = uicontrol('parent', bg,'style', 'radiobutton', 'position', rbbotpos,'string', '  DeltaF/F', 'Value', 0 ,'fontsize', fontsize, 'handlevisibility', 'off');
            toprb = uicontrol('parent', bg, 'style', 'radiobutton', 'position', rbtoppos,'string', '  Average', 'Value', 1, 'fontsize', fontsize, 'handlevisibility', 'off');
            bg.Visible = 'on';
            if mod(iRoi,sprow) == 0 && iRoi < numrois
                iFig = iFig+1; iSP = 1;
                figure('Position', positionfig, 'Name', sprintf('Traces %i',iFig));
                savebut = uicontrol('parent', gcf, 'style', 'pushbutton', 'position', savebutpos,'string', 'Save selected','fontsize', fontsize, 'callback', {@cb_savetrcbut});
                closebut =  uicontrol('parent', gcf, 'style', 'pushbutton', 'position', closebutpos,'string', 'Close (no saving)','fontsize', fontsize, 'callback', {@cb_closetrcbut});
            else
                iSP = iSP+1;
            end
        end
    end

    function [tmpvals,tmpdFoF, yrange] = calc_roi_av_trace(roiidx, vals, pr, mode)
        roimask = roiInfo(roiidx).mask;
        if strcmp(mode,'ROI')
            roipr = roiInfo(roiidx).plotrange;
            if any(roipr ~= pr)
                if pr(1) <= roipr(1) && pr(2) <= roipr(1), newmask = false(size(roimask));
                elseif pr(1) <= roipr(1)&& pr(2) <= roipr(2), ycut = 1:pr(2)-roipr(1)+1; ypaste = roipr(1)-pr(1)+1:pr(2)-pr(1)+1;
                elseif pr(1) <= roipr(1) && pr(2) >= roipr(2), ycut = 1:roipr(2)-roipr(1)+1; ypaste = roipr(1)-pr(1)+1:roipr(2)-pr(1)+1;
                elseif pr(1) >= roipr(1) && pr(2) <= roipr(2), ycut = pr(1)-roipr(1)+1:pr(2)-roipr(1)+1; ypaste = pr(1)-pr(1)+1:pr(2)-pr(1)+1;
                elseif pr(1) >= roipr(1) && pr(2) >= roipr(2), ycut = pr(1)-roipr(1)+1:roipr(2)-roipr(1)+1; ypaste = pr(1)-pr(1)+1:roipr(2)-pr(1)+1;
                end
                newmask = false(size(vals));
                newmask(ypaste,:) = roimask(ycut,:);
            else
                newmask = roimask;
            end
            % Calculate ROI average
            collapsemask = sum(newmask,2);
            ymin = find(collapsemask,1,'first'); ymax = find(collapsemask,1,'last');
        else % Across all timepoints
            newmask = roimask;
            collapsemask = sum(newmask,1); % Find first and last x value in mask
            xmin = find(collapsemask,1,'first'); xmax = find(collapsemask,1,'last');
            ymin = 1; ymax = size(vals,1);
            newmask = false(size(vals));
            newmask(ymin:ymax,xmin:xmax) = true;
        end
            tmpvals = zeros(ymax-ymin+1,1);
            ival = 1;
            for iY = ymin:ymax
                tmpvals(ival,1) = mean(vals(iY,newmask(iY,:)));
                ival = ival+1;
            end
            % dFoF
            tmpvalav = mean(tmpvals);
            tmpdFoF = (tmpvals-tmpvalav)/tmpvalav;
            yrange = [ymin ymax];
    end

    function save_files(sp, fi, mode)
        diafig = figure('Position',[positionfig(1)+positionfig(3)/2 positionfig(2)+positionfig(4)/2 200 100],'Name', 'Saving','toolbar', 'none', 'menu', 'none');
        set(gca, 'Color', 'none'); axis off;
        diatxtinfo = {'Files are saved', 'this may take a while...'};
        diatxt = uicontrol('parent', diafig, 'style', 'text', 'position', [10 5 180 80], 'string', diatxtinfo, 'fontsize', fontsize+1);
        pause(0.5);
        
        % Get parameters
        pr = figInfo(fi).plotrange;
        wsz = figInfo(fi).avwinsize;
        
        raw = composite2D;
        averaged = average_linescan(raw, wsz);
        dFoF = deltaFovF_linescan(averaged);
        
        if ~figInfo(fi).saved
            % Save ROI-unrelated stuff
            % .csv
            set(diatxt,'string',[diatxtinfo,{'','...create .csv files.'}]);
            pause(0.2);
            writematrix(raw, strcat(sp,'raw.csv'));
            writematrix(averaged, strcat(sp,'AV_bin_',num2str(wsz), '.csv'));
            writematrix(dFoF, strcat(sp,'dFoF_bin_',num2str(wsz), '.csv'));
            % .tifs
            set(diatxt,'string',[diatxtinfo,{'','...create .tif files.'}]);
            pause(0.2);
            vals = uint16(raw(pr(1):pr(2),:));
            imwrite(vals, strcat(sp,'raw_range_', num2str(pr(1)),'-',num2str(pr(2)),'.tif'), 'tif');
            vals = uint16(averaged(pr(1):pr(2),:));
            imwrite(vals, strcat(sp,'AV_','bin_',num2str(wsz),'_range_', num2str(pr(1)),'-',num2str(pr(2)),'.tif'), 'tif');
            vals = uint16(dFoF(pr(1):pr(2),:));
            imwrite(vals, strcat(sp,'dFovF_','bin_',num2str(wsz),'_range_', num2str(pr(1)),'-',num2str(pr(2)),'.tif'), 'tif');
            figInfo(fi).saved = true;
        end
        
        % Determine ROIs to save
        if strcmp(mode, 'all')
            saverois = 1:size(roiInfo,2);
        else
            saverois = find([traceInfo(:).save] == 1);
        end
        
        for iRoi = 1:numel(saverois)
            set(diatxt,'string',[diatxtinfo,{'','...save ROI # ', num2str(iRoi)}]);
            pause(0.2);
            val = averaged(pr(1):pr(2),:);
%             valsc = scale_data(val, climraw);
            valroi = mark_rois(uint8(val),pr, saverois(iRoi));
            imwrite(valroi, strcat(sp,'AV_ROI_',num2str(iRoi),'.tif'), 'tif');
            if roiInfo(iRoi).mode == 1 % Line
                [avvals,dFoFvals,~] = calc_roi_av_trace(saverois(iRoi), averaged, [], 'whole');
                writematrix(avvals, strcat(sp,'ROI_',num2str(saverois(iRoi)),'_AV_timetotal.csv'));
                writematrix(dFoFvals, strcat(sp,'ROI_',num2str(saverois(iRoi)),'_dFoF_timetotal.csv'));
            else
                [avvals,dFoFvals,~] = calc_roi_av_trace(saverois(iRoi), val, pr, 'ROI');
                writematrix(avvals, strcat(sp,'ROI_',num2str(saverois(iRoi)),'_AV.csv'));
                writematrix(dFoFvals, strcat(sp,'ROI_',num2str(saverois(iRoi)),'_dFoF.csv'));
            end
        end
        close(diafig);
        okbox = msgbox('Files saved');
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
        del = [roiInfo(:).selected] == true;
        nroi = size(roiInfo,2);
        if isempty(roiInfo(end).position), del = del(1:end-1); nroi = nroi-1; end
        if ~isempty(traceInfo.figID)
            deltrc = zeros(numel(del),1);
            for iE = 1:numel(del)
                deltrc(iE) = [traceInfo(:).roiID] == roiInfo(del(iE)).ID;
            end
            traceInfo = traceInfo(~deltrc);
        end
        if sum(del) == nroi
            roiInfo = roiInfo(1); roiInfo(1).mask = []; roiInfo(1).position = []; roiInfo(1).name = []; roiInfo(1).ID = []; roiInfo(1).selected = []; roiInfo(1).saved = []; roiInfo(1).mode = []; roiInfo(1).plotrange = [];
        else
            roiInfo = roiInfo(~del);
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
        if sum(roiInfo(:).selected) == 0
            msgbox('Please select a ROI first!');
        end
        traces_overview();
    end

    function cb_subplot(~,~,tmpmtrx, figID, cscui)
        create_rois(uint8(tmpmtrx), figID, cscui);
        
    end

    function cb_drawroi(hObj, ~, tmpax, tmpmtrx, pr)
        if isempty(roiInfo(1).ID)
            roiInfo(1).ID = 1;
            roiInfo(1).name = 0; roiInfo(1).mask = 0; roiInfo(1).position = 0; roiInfo(1).selected = true; roiInfo(1).saved = 0; roiInfo(1).mode = 0; roiInfo(1).plotrange = 0;
        end
        roicnt = roiInfo(end).ID;
        if strncmp(hObj.String, 'Line',4)
            roiInfo(roicnt).mask = false(size(tmpmtrx));
            roi1 = images.roi.Line(tmpax, 'linewidth', 2, 'InteractionsAllowed','none');
            draw(roi1);
            roi2 = images.roi.Line(tmpax, 'linewidth', 2, 'InteractionsAllowed','none');
            draw(roi2);
            x = sort([round(mean(roi1.Position(:,1))) round(mean(roi2.Position(:,1)))], 'ascend');
            roiInfo(roicnt).mask(:,x(1):x(2)) = true;
            roiInfo(roicnt).position = roi1.Position(1,:);
            roiInfo(roicnt).mode = 1;
            roiInfo(roicnt).plotrange = pr;
        elseif strncmp(hObj.String, 'Rect',4)
            roi = images.roi.Rectangle(tmpax, 'FaceAlpha', 0.2, 'InteractionsAllowed','none');
            draw(roi);
            roiInfo(roicnt).mask = createMask(roi);
            roiInfo(roicnt).position = roi.Position(1:2);
            roiInfo(roicnt).mode = 2;
            roiInfo(roicnt).plotrange = pr;
        else
            roi = images.roi.Freehand(tmpax, 'closed', true, 'multiclick', true, 'InteractionsAllowed','none','linewidth', 1);
            draw(roi);
            roiInfo(roicnt).mask = createMask(roi);
            roiInfo(roicnt).position = roi.Position(1,:);
            roiInfo(roicnt).mode = 3;
            roiInfo(roicnt).plotrange = pr;
        end
    end

    function cb_addroi(~, ~,fp, tmpmtrx, pr, cscui)
        roicnt = roiInfo(end).ID;
        if isempty(roiInfo(roicnt).position)
            msgbox('Please select a ROI!');
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
            roiInfo(roicnt+1).ID = roicnt+1;
        end
    end

    function cb_delroi(~, ~,tmpmtrx, fp, cscui)
        roiInfo = roiInfo(1);
        roiInfo(1).mask = []; roiInfo(1).position = []; roiInfo(1).name = []; roiInfo(1).ID = []; roiInfo(1).selected = []; roiInfo(1).mode = []; roiInfo(1).plotrange = [];
        figures = get(groot,'Children');
        tmpfig = figures(strncmp({figures(:).Name}, 'Line',4));
        roilb = findobj(tmpfig, 'type', 'uicontrol', 'style', 'listbox');
        set(roilb, 'string', roiInfo.name);
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

    function cb_traceswitchrb(~,evdat, traceidx, roi, sp, xrange)
        seltrc = get(evdat.NewValue, 'String');
        subplot(sp);
        if strcmp(seltrc, '  DeltaF/F')
            plot(traceInfo(traceidx(roi)).timestamp{1}, traceInfo(traceidx(roi)).dFoF_roi_av{1},'linewidth', 1, 'color', 'blue');
            ylabel(trcylabeldFoF); set(gca, 'fontsize', fontsize);
        else
            plot(traceInfo(traceidx(roi)).timestamp{1}, traceInfo(traceidx(roi)).binned_roi_av{1},'linewidth', 1, 'color', 'black');
            ylabel(trcylabel); set(gca, 'fontsize', fontsize);
        end
        xlabel(trcxlabel); xlim(xrange);
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
            savepointer = strcat(savepath,'\', fname,'_');
            % Prepare for saving
            save_files(savepointer, figidx, 'all');
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

    function cb_savetrcrb(hObj, ~, traceidx, iRoi)
        traceInfo(traceidx(iRoi)).save = hObj.Value;
    end

    function cb_savetrcbut(~,~,figidx)
        if isempty(savepointer)
            % Saving path and name
            [spath, fname, ~] = fileparts(filepointer);
            savepath = uigetdir(spath);
            savepointer = strcat(savepath,'\', fname,'_');
        end
        close gcf;
        save_files(savepointer, figidx, 'selected');
    end

    function cb_closetrcbut(~,~)
        for iE = 1: size(traceInfo,2), traceInfo(iE).save = false; end
        close gcf;
    end

end

