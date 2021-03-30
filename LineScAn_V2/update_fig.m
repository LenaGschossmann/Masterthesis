function update_fig()
% Plots and updates overview of raw, averaged, and dFoF image

% Declare globally shared variables
global WHFIG FONTSIZE SCMAP WINSZ PLOTRANGE CLIMRAW CLIMUI figINFO...
    ROICNTER IMHW COMPOSITE2D FTIME dFWIN

deltawinsz = round(dFWIN/FTIME);

figures = get(groot,'Children');
findfig = strncmp({figures(:).Name}, 'Fig',3);
currfig = figures(find(findfig,1, 'first'));
currcnt = get(currfig, 'UserData');
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
rawmtrx = COMPOSITE2D(PLOTRANGE(1):PLOTRANGE(2),:);
averaged = average_linescan(rawmtrx, WINSZ);
lin_averaged = averaged';
[~, dFoF] = rollBase_dFoF(lin_averaged,deltawinsz,size(lin_averaged,2), 'roll');
dFoF = dFoF';

% Plot data
scraw = scale_data(rawmtrx, CLIMRAW);
scav = scale_data(averaged, CLIMRAW);
dFoF(dFoF<CLIMUI(1)) = CLIMUI(1); dFoF(dFoF>CLIMUI(2)) = CLIMUI(2);
dFoF = scale_data(dFoF, CLIMUI);

% Check if any ROIs already created and add to images
if ROICNTER > 0
    scrawroi = mark_rois(scraw,PLOTRANGE,'all');
    scavroi = mark_rois(scav,PLOTRANGE,'all');
    dFoFroi = mark_rois(dFoF, PLOTRANGE,'all');
else
    scrawroi=scraw;
    scavroi=scav;
    dFoFroi=dFoF;
end

%% Plot
set(0,'currentfig',currfig);
figID = get(currfig,'UserData');
colormap(SCMAP);
% Raw image
subplot('Position', positionspul, 'parent', currfig);
im1 = imagesc(scrawroi, CLIMRAW);
axis off;
spname = 'Raw';
title(spname);
colorbar;
set(im1, 'ButtonDownFcn', {@cb_subplot, scraw, figID, CLIMRAW});
% Averaged image
subplot('Position', positionspur, 'parent', currfig);
im2 = imagesc(scavroi, CLIMRAW);
axis off;
spname = 'Averaged';
title(spname);
colorbar;
set(im2, 'ButtonDownFcn', {@cb_subplot, scav, figID, CLIMRAW});
% DeltaFovF
subplot('Position', positionspbl, 'parent', currfig);
im3 = imagesc(dFoFroi, CLIMUI);
axis off;
spname = 'Delta F over F';
title(spname);
colorbar;
set(im3, 'ButtonDownFcn', {@cb_subplot,dFoF, figID, CLIMUI});

%% Display Info
currinfo = {'Information:', '', sprintf('Plotrange: %i - %i', figINFO(figID).plotrange),...
    '', sprintf('Image size (post binning): W: %i | H: %i',IMHW(2), IMHW(1)),...
    '', sprintf('Line averaging window size: %i',figINFO(figID).avwinsize),...
    '', sprintf('Colormap: %s', figINFO(figID).cscmap),...
    '', sprintf('Color scale limits: %i - %i', figINFO(figID).csclimits),...
    '', sprintf('Time per frame: %s s', num2str(round(FTIME,3)))};
positionspbr = positionspbr.*[WHFIG WHFIG]; positionspbr(2) = round(positionspbr(2)*0.75);
infotxt = uicontrol('style','text', 'parent', currfig, 'position', positionspbr, 'String', currinfo,'FONTSIZE', FONTSIZE);


%% Local callbacks
    function cb_subplot(~,~,tmpmtrx, figID, csclims)
        create_rois(tmpmtrx, figID, csclims);
    end

end