function create_rois(tmpmtrx, figID, cscui)

% Declare globally shared variables
global WHFIG POSITIONFIG POSITIONROISELECT FONTSIZE CLIMRAW figINFO roiINFO IMHW

fp = find([figINFO(:).IDs] == figID);
if isempty(roiINFO(1).ID), roitot = 0; else, roitot = size(roiINFO,2)-1; end
pr = figINFO(fp).plotrange;

% Display parameteres
vspace = [20 20]; % [small bigg bigger]
hspace = [10 40 25]; % [ctr-aligned rest]
whbut = [80 20];
whtxt = [100 15];

% Positions
lineroipos = [hspace(2) vspace(1) whbut];
rectroipos = [lineroipos(1)+lineroipos(3)+hspace(1) vspace(1) whbut];
closeroipos = [WHFIG(1)-whbut(1)-hspace(1)*4 vspace(1) whbut];
delroipos = [closeroipos(1)-hspace(1)-whbut(1) vspace(1) whbut];
addroipos = [delroipos(1)-hspace(1)-whbut(1) vspace(1) whbut];
roistxtpos = [addroipos(1)-hspace(1)-whtxt(1) vspace(1) whtxt];

% Main window
roifig = figure('Position', POSITIONFIG, 'Name', 'ROI selection');
set(gca, 'Color', 'none'); axis off;
roiax = subplot('Position', POSITIONROISELECT, 'parent', roifig);
colormap(figINFO(fp).cscmap);
if cscui, cscui = figINFO(fp).csclimits; else, cscui = CLIMRAW; end
if ~isempty(roiINFO(1).ID), tmpmtrxroi = mark_rois(tmpmtrx, pr, 'all'); else, tmpmtrxroi=tmpmtrx; end
imagesc(tmpmtrxroi, cscui);
axis off;

% Add UI controls
lineroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', lineroipos,'string', 'Line ROI','FONTSIZE', FONTSIZE, 'callback', {@cb_drawroi,roiax, pr});
rectroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', rectroipos,'string', 'Rectangle ROI','FONTSIZE', FONTSIZE, 'callback', {@cb_drawroi,roiax, pr});
addroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', addroipos,'string', 'Add ROI','FONTSIZE', FONTSIZE, 'callback', {@cb_addroi, fp, tmpmtrx, pr, cscui});
delroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', delroipos,'string', 'Delete all','FONTSIZE', FONTSIZE, 'callback', {@cb_delroi,tmpmtrx, fp, cscui});
closeroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', closeroipos,'string', 'Close','FONTSIZE', FONTSIZE, 'callback', {@cb_closeroi});
roistxt  = uicontrol('parent', roifig, 'style', 'text','position', roistxtpos,'string', sprintf('ROIs Total # %i', roitot),'FONTSIZE', FONTSIZE);


%% Local Callbacks
    function cb_drawroi(hObj, ~, tmpax, pr)
        try
            if isempty(roiINFO(1).ID)
                roiINFO(1).ID = 1;
                roiINFO(1).name = 0; roiINFO(1).mask = 0; roiINFO(1).position = 0; roiINFO(1).selected = true; roiINFO(1).saved = 0; roiINFO(1).mode = 0; roiINFO(1).plotrange = 0;
            end
            roicnt = size(roiINFO,2);
            if strncmp(hObj.String, 'Line',4)
                roiINFO(roicnt).mask = false(IMHW);
                roi1 = images.roi.Line(tmpax, 'linewidth', 2, 'InteractionsAllowed','none');
                draw(roi1);
                roi2 = images.roi.Line(tmpax, 'linewidth', 2, 'InteractionsAllowed','none');
                draw(roi2);
                x = sort([round(mean(roi1.Position(:,1))) round(mean(roi2.Position(:,1)))], 'ascend');
                roiINFO(roicnt).mask(:,x(1):x(2)) = true;
                roiINFO(roicnt).position = roi1.Position(1,:);
                roiINFO(roicnt).mode = 1;
                roiINFO(roicnt).plotrange = [1 IMHW(1)];
            elseif strncmp(hObj.String, 'Rect',4)
                roi = images.roi.Rectangle(tmpax, 'FaceAlpha', 0.2, 'InteractionsAllowed','none');
                draw(roi);
                roiINFO(roicnt).mask = createMask(roi);
                roiINFO(roicnt).position = roi.Position(1:2);
                roiINFO(roicnt).mode = 2;
                roiINFO(roicnt).plotrange = pr;
            end
        end
    end

    function cb_addroi(~, ~,fp, tmpmtrx, pr, cscui)
        roicnt = size(roiINFO,2);
        if isempty(roicnt), roicnt = 1; end
        if isempty(roiINFO(roicnt).position)
            msgbox('Please draw a ROI first!');
        else
            roiINFO(roicnt).name = sprintf('ROI #%i',roiINFO(roicnt).ID);
            roiINFO(roicnt).saved = 0;
            if roicnt == 1, roiINFO(roicnt).selected = true;
            else, roiINFO(roicnt).selected = false;
            end
            figures = get(groot,'Children');
            tmpfig = figures(strncmp({figures(:).Name}, 'Line',4));
            roilb = findobj(tmpfig, 'type', 'uicontrol', 'style', 'listbox');
            set(roilb, 'string', {roiINFO.name});
            tmpfig = figures(1);
            roistxt = findobj(tmpfig, 'type', 'uicontrol', 'style', 'text');
            set(roistxt, 'String', sprintf('ROIs Total # %i', size(roiINFO,2)));
            tmpmtrx = mark_rois(tmpmtrx, pr, 'all');
            subplot('Position', POSITIONROISELECT, 'parent', tmpfig);
            colormap(figINFO(fp).cscmap);
            imagesc(tmpmtrx, cscui);
            axis off;
            for iC = 1:size(roiINFO,2)
                if ~isempty(roiINFO(iC).position)
                    tmppos = [WHFIG.*POSITIONROISELECT(1:2)+...
                        roiINFO(iC).position.*(WHFIG.*POSITIONROISELECT(3:4)./fliplr(size(tmpmtrx))) 15 15];
                    tmppos(2) = WHFIG(2)-tmppos(2);
                    uicontrol('style','text','parent',tmpfig,'position', tmppos,'string',num2str(iC), 'fontsize', FONTSIZE)
                end
            end
            roiINFO(roicnt+1).ID = max([roiINFO(:).ID])+1;
        end
    end

    function cb_delroi(~, ~,tmpmtrx, fp, cscui)
        roiINFO = roiINFO(1);
        roiINFO(1).mask = []; roiINFO(1).position = []; roiINFO(1).name = []; roiINFO(1).ID = []; roiINFO(1).selected = []; roiINFO(1).mode = []; roiINFO(1).plotrange = [];
        figures = get(groot,'Children');
        tmpfig = figures(strncmp({figures(:).Name}, 'Line',4));
        roilb = findobj(tmpfig, 'type', 'uicontrol', 'style', 'listbox');
        set(roilb, 'string', ''); set(roilb, 'value', 1);
        tmpfig = figures(1);
        roistxt = findobj(tmpfig, 'type', 'uicontrol', 'style', 'text');
        delete(roistxt(1:end-1));
        set(roistxt(end), 'String', '');
        subplot('Position', POSITIONROISELECT, 'parent', tmpfig);
        colormap(figINFO(fp).cscmap);
        imagesc(tmpmtrx, cscui);
        axis off;
    end

    function cb_closeroi(~,~)
        figures = get(groot,'Children');
        close(figures(1));
        update_fig();
    end

end