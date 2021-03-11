function [val_mtrx, im_cell] = define_bg_threshold(aiplist, selfiles, bg_perc)
%% Function to manually define background threshold value

im_cell = cell(numel(selfiles),3);
val_mtrx = NaN(numel(selfiles),3);

for iSel = 1:numel(selfiles)
    avim = aiplist{selfiles(iSel)};
    val_mtrx(iSel,2:3) = [min(avim,[],'all') max(avim,[],'all')]; 
    [c_abs,e] = histcounts(avim); % Count frequency of intensity values
    c_rel = c_abs./sum(c_abs); % Relative frequency
    c_cum = cumsum(c_rel); % Cummulative frequency
    [~,idx] = min(abs(c_cum - bg_perc));
    val_mtrx(iSel,1) = e(idx);
    im_cell{iSel,1} = avim;
    im_cell{iSel,2} = (avim./val_mtrx(iSel,3)).*2^16;
    im_cell{iSel,3} = true(size(avim,1),size(avim,2));
    im_cell{iSel,3}(avim > val_mtrx(iSel,1)) = false;
end

cf = 1;

%% GUI
scrsz = get(0,'ScreenSize');
ctrl_pos_2 = [round(scrsz(3)/3) round(scrsz(4)/5) 410 50];
sl_pos = [5 5 250 30];
but1_pos = [260 0 40 30];
but2_pos = [305 0 40 30];
but3_pos = [350 0 40 30];
bg_fig = figure('toolbar', 'none', 'menu', 'none');
ctrl_fig = figure('Position', ctrl_pos_2,'toolbar', 'none', 'menu', 'none');
sl = uicontrol('Style','slider','Min',val_mtrx(cf,2),'Max',val_mtrx(cf,3),...
    'SliderStep',[1 1]./100,'Value',val_mtrx(cf,1), 'Position', sl_pos, 'Callback', {@sl_cb, bg_fig});
but1 = uicontrol('parent', ctrl_fig, 'style','pushbutton', 'Position', but1_pos, 'String', '<-', 'Callback', {@but12_cb});
but2 = uicontrol('parent', ctrl_fig, 'style','pushbutton', 'Position', but2_pos, 'String', '->', 'Callback', {@but12_cb});
but3 = uicontrol('parent', ctrl_fig, 'style','pushbutton', 'Position', but3_pos, 'String', 'OK', 'Callback', {@but3_cb});
if numel(selfiles) == 1, but1.Enable = 'off'; but2.Enable = 'off'; end

figure(bg_fig);
colormap('gray');
imagesc(uint16(im_cell{cf,2})); axis('off');
hold on, imagesc(uint16(im_cell{cf,3}.*2^16), 'AlphaData',.4);
uiwait;

%% Local function
function sl_cb(hOb,~, bg_fig)
val_mtrx(cf,1) = round(hOb.Value);
im_cell{cf,3}(:,:) = true;
im_cell{cf,3}(im_cell{cf,1} > val_mtrx(cf,1)) = false;
figure(bg_fig);
clf;
colormap('gray');
imagesc(uint16(im_cell{cf,2}));axis('off');
hold on, imagesc(uint16(im_cell{cf,3}.*2^16), 'AlphaData',.4);
end

function but12_cb(hOb,~)
    if strcmp(hOb.String, '<-')
        if cf > 1, cf = cf-1; end
    else
        if cf < numel(selfiles), cf = cf+1; end
    end
    sl.Min = val_mtrx(cf,2);
    sl.Max = val_mtrx(cf,3);
    sl.SliderStep = [1 1]./val_mtrx(cf,3);
    sl.Value = val_mtrx(cf,1);
    figure(bg_fig);
    clf;
    colormap('gray');
    imagesc(uint16(im_cell{cf,2}));axis('off');
    hold on, imagesc(uint16(im_cell{cf,3}.*2^16), 'AlphaData',.4);
end

function but3_cb(~,~)
uiresume;
figures = get(groot,'Children');
findfig = strncmp({figures(:).Name}, 'ROI',3);
close(figures(~findfig));
end


end