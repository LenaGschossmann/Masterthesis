function create_rois(tmpmtrx, csclims)
% Opens a figure with the linescan to let the user draw ROIs (line ROI:
% allows user to select an ROI of certain width that covers the total
% selected range)
% Input:
% tmpmtrx: (2D matrix) Linescan image on which ROI is drawn
% figID: ID of calling figure to obtain plotting parameters
% csclims: ([min max]) Colorscale limits for plotting

% Declare globally shared variables
global WHFIG POSITIONFIG POSITIONROISELECT FONTSIZE roiINFO...
    ROICNTER IMHW SCMAP

% Initialize variables and unpack
newroi = [];
addstate = false;

% Display parameteres
vspace = [20 20]; % [small bigg bigger]
hspace = [10 40 25]; % [ctr-aligned rest]
whbut = [80 20];
whtxt = [100 15];

% Positions
lineroipos = [hspace(2) vspace(1) whbut];
closeroipos = [WHFIG(1)-whbut(1)-hspace(1)*4 vspace(1) whbut];
delroipos = [closeroipos(1)-hspace(1)-whbut(1) vspace(1) whbut];
addroipos = [delroipos(1)-hspace(1)-whbut(1) vspace(1) whbut];
roistxtpos = [addroipos(1)-hspace(1)-whtxt(1) vspace(1) whtxt];

% Main window
roifig = figure('Position', POSITIONFIG, 'Name', 'ROI selection');
set(gca, 'Color', 'none'); axis off;
roiax = subplot('Position', POSITIONROISELECT, 'parent', roifig);
colormap(SCMAP);
if ROICNTER > 0, tmpmtrxroi = mark_rois(tmpmtrx, 'all'); else, tmpmtrxroi=tmpmtrx; end
imagesc(tmpmtrxroi, csclims);
axis off;

% Add UI controls
lineroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', lineroipos,'string', 'Line ROI','FONTSIZE', FONTSIZE, 'callback', {@cb_drawroi,roiax});
addroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', addroipos,'string', 'Add ROI','FONTSIZE', FONTSIZE, 'callback', {@cb_addroi, tmpmtrx});
delroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', delroipos,'string', 'Delete last','FONTSIZE', FONTSIZE, 'callback', {@cb_delroi,tmpmtrx});
closeroi = uicontrol('parent', roifig, 'style', 'pushbutton', 'position', closeroipos,'string', 'Close','FONTSIZE', FONTSIZE, 'callback', {@cb_closeroi});
roistxt  = uicontrol('parent', roifig, 'style', 'text','position', roistxtpos,'string', sprintf('ROIs Total # %i', ROICNTER),'FONTSIZE', FONTSIZE);

%% Local Callbacks
    function cb_drawroi(~, ~, tmpax)
        addstate = false;
        newroi = ROICNTER+1;
        roiINFO(newroi).mask = false(IMHW);
        roi1 = images.roi.Line(tmpax, 'linewidth', 2, 'InteractionsAllowed','none');
        draw(roi1);
        roi2 = images.roi.Line(tmpax, 'linewidth', 2, 'InteractionsAllowed','none');
        draw(roi2);
        x = sort([round(mean(roi1.Position(:,1))) round(mean(roi2.Position(:,1)))], 'ascend');
        roiINFO(newroi).mask(:,x(1):x(2)) = true;
        roiINFO(newroi).position = [x(1) 1 diff(x)+1 IMHW(1)];
        roiINFO(newroi).plotrange = [1 IMHW(1)];
    end

    function cb_addroi(~, ~, tmpmtrx)
        if isempty(newroi)
            msgbox('Please draw a ROI first!');
        else
            roiINFO(newroi).name = sprintf('ROI #%i',newroi);
            roiINFO(newroi).saved = 0;
            if newroi == 1, roiINFO(newroi).selected = true; else, roiINFO(newroi).selected = false; end
            roiINFO(newroi).ID = newroi;
            ROICNTER = ROICNTER+1;
            addstate = true;
            newroi = [];
            figures = get(groot,'Children');
            tmpfig = figures(strncmp({figures(:).Name}, 'Line',4));
            roilb = findobj(tmpfig, 'type', 'uicontrol', 'style', 'listbox');
            set(roilb, 'string', {roiINFO.name});
            tmpfig = figures(1);
            roistxt = findobj(tmpfig, 'type', 'uicontrol', 'style', 'text');
            set(roistxt, 'String', sprintf('ROIs Total # %i', ROICNTER));
            tmpmtrx = mark_rois(tmpmtrx, 'all');
            subplot('Position', POSITIONROISELECT, 'parent', tmpfig);
            colormap(SCMAP);
            imagesc(tmpmtrx, csclims);
            axis off;
        end
    end

    function cb_delroi(~, ~,tmpmtrx)
        if addstate
            if ROICNTER > 1,roiINFO = roiINFO(1:ROICNTER-1);
            else, roiINFO = struct('name', [], 'mask', [], 'position', [], 'ID', [], 'selected', [], 'saved', [], 'plotrange', []);
            end
            ROICNTER = ROICNTER-1;
            figures = get(groot,'Children');
            tmpfig = figures(strncmp({figures(:).Name}, 'Line',4));
            roilb = findobj(tmpfig, 'type', 'uicontrol', 'style', 'listbox');
            set(roilb, 'string', ''); set(roilb, 'value', 1);
            tmpfig = figures(1);
            roistxt = findobj(tmpfig, 'type', 'uicontrol', 'style', 'text');
            set(roistxt, 'String',sprintf('ROIs Total # %i', ROICNTER));
            subplot('Position', POSITIONROISELECT, 'parent', tmpfig);
            colormap(SCMAP);
            imagesc(tmpmtrx, csclims);
            axis off;
        end
        subplot(roiax);
        hold off;
        colormap(SCMAP);
        if ROICNTER > 0, tmpmtrxroi = mark_rois(tmpmtrx, 'all'); else, tmpmtrxroi=tmpmtrx; end
        imagesc(tmpmtrxroi, csclims);
        axis off;
    end

    function cb_closeroi(~,~)
        figures = get(groot,'Children');
        close(figures(1));
        update_fig();
    end

end