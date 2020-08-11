function update_fig()

% Declare globally shared variables
global WHFIG POSITIONFIG FONTSIZE SCMAP WINSZ PLOTRANGE CLIMRAW CLIMUI figINFO roiINFO FIGCOUNTER IMHW COMPOSITE2D FTIME

figures = get(groot,'Children');
if numel(figures) == 1
    FIGCOUNTER = FIGCOUNTER+1;
    currfig = figure('Position', POSITIONFIG);
    set(gcf, 'UserData', FIGCOUNTER);
    set(gcf, 'Name', sprintf('Figure #%i', FIGCOUNTER));
    currcnt = FIGCOUNTER;
    figINFO(currcnt).IDs = FIGCOUNTER;
    figINFO(currcnt).name = sprintf('Figure #%i', FIGCOUNTER);
    figINFO(currcnt).saved = false;
else
    findfig = strncmp({figures(:).Name}, 'Fig',3);
    currfig = figures(find(findfig,1, 'first'));
    currcnt = get(currfig, 'UserData');
end
figINFO(currcnt).plotrange = PLOTRANGE;
figINFO(currcnt).avwinsize = WINSZ;
figINFO(currcnt).cscmap = SCMAP;
figINFO(currcnt).csclimits = CLIMUI;

% Parameters
spcol = 2; sprow = 2;
hspace = 30;
vspace = 60;
whsp = [floor((WHFIG(1)-3*hspace)/spcol) floor((WHFIG(2)-3*vspace)/sprow)];
whsp = whsp./WHFIG; hspace = hspace/WHFIG(1); vspace = vspace/WHFIG(2);
positionspbl = [hspace vspace whsp];
positionspbr = [positionspbl(1)+positionspbl(3)+hspace vspace whsp];
positionspul = [hspace positionspbl(2)+positionspbl(4)+vspace whsp];
positionspur = [positionspul(1)+positionspul(3)+hspace positionspbl(2)+positionspbl(4)+vspace whsp];

% Call Analysis functions
if PLOTRANGE(2)>IMHW(1), PLOTRANGE=[1 IMHW(1)]; end % Make sure the plot range doesnt exceed sampling timepoints
rawmtrx = COMPOSITE2D(PLOTRANGE(1):PLOTRANGE(2),:);
averaged = average_linescan(rawmtrx, WINSZ);
dFoF = deltaFovF_linescan(averaged);

% Plot
set(0,'currentfig',currfig);
figID = get(currfig,'UserData');
colormap(SCMAP);
% Raw (scaled to 255)
%         scraw = scale_data(raw, CLIMRAW);
scraw = uint8(rawmtrx);
if ~isempty(roiINFO(1).ID), scrawroi = mark_rois(scraw,PLOTRANGE,'all'); else, scrawroi=scraw; end
sp1 = subplot('Position', positionspul);
im1 = imagesc(scrawroi, CLIMRAW);
axis off;
spname = 'Raw';
title(spname);
colorbar;
set(im1, 'ButtonDownFcn', {@cb_subplot, scraw, figID, false});

% Averaged (scaled to 255)
%         scav = scale_data(averaged, CLIMRAW);
scav = uint8(averaged);
if ~isempty(roiINFO(1).ID), scavroi = mark_rois(scav,PLOTRANGE,'all'); else, scavroi=scav; end
sp2 = subplot('Position', positionspur, 'parent', currfig);
im2 = imagesc(scavroi, CLIMRAW);
axis off;
spname = 'Averaged';
title(spname);
colorbar;
set(im2, 'ButtonDownFcn', {@cb_subplot, scav, figID, false});

% DeltaFovF (scaled to user-defined scale)
dFoF(dFoF<CLIMUI(1)) = CLIMUI(1); dFoF(dFoF>CLIMUI(2)) = CLIMUI(2);
dFoF = scale_data(dFoF, CLIMUI);
if ~isempty(roiINFO(1).ID), dFoFroi = mark_rois(dFoF, PLOTRANGE,'all'); else, dFoFroi=dFoF; end
sp3 = subplot('Position', positionspbl, 'parent', currfig);
im3 = imagesc(dFoFroi, CLIMUI);
axis off;
spname = 'Delta F over F';
title(spname);
colorbar;
set(im3, 'ButtonDownFcn', {@cb_subplot,dFoF, figID, true});

% Info
currinfo = {'Information:', '', sprintf('Plotrange: %i - %i', figINFO(figID).plotrange),...
    '', sprintf('Line averaging window size: %i',figINFO(figID).avwinsize),...
    '', sprintf('Colormap: %s', figINFO(figID).cscmap),...
    '', sprintf('Color scale limits: %i - %i', figINFO(figID).csclimits),...
    '', sprintf('Time per frame: %s s', num2str(round(FTIME,3)))};
positionspbr = positionspbr.*[WHFIG WHFIG]; positionspbr(2) = round(positionspbr(2)*0.75);
infotxt = uicontrol('style','text', 'parent', currfig, 'position', positionspbr, 'String', currinfo,'FONTSIZE', FONTSIZE);


%% Local callbacks
    function cb_subplot(~,~,tmpmtrx, figID, cscui)
        create_rois(tmpmtrx, figID, cscui);
    end


end