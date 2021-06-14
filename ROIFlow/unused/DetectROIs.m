%% Script for automatic ROI detection in iGluSnFR fluorescence data %%

%% Add paths
wd = pwd(); p = split(wd,'\'); p = p(1:end-1); pp = []; for ii=1:numel(p), pp = fullfile(pp,p{ii}); end
pp=fullfile(pp,'OpenBFiles'); addpath(genpath(pp)), clear('p','pp');

%% Parameters
[params] = get_predefined_params('perc_fac', 1.05, 'connc_px_thresh', 20);
global bg_thresh;

t_all = tic;

%% Import timeseries and ask for frame time
[files, paths] = uigetfile({'*.tif*'; '*.lsm'; '*.nd2'}, [],'Multiselect', 'off');
if ~isa(files, 'cell'), files = {files}; end
if ~isa(paths, 'cell'), paths = {paths}; end
fullfilename = strcat(paths{1},files{1});

if strcmp(fullfilename(end-3:end),'.lsm')
    [inputim,~] = load_bf_file(fullfilename,false);
    imsize = [size(inputim,2) size(inputim,1)];
    num_px = imsize(1)*imsize(2);
    n_frames = size(inputim,3);
    destpath = strcat(paths{1}, 'dF_', files{1}(1:end-3),"tif");
else
    Info=imfinfo(fullfilename);
    imsize=[Info(1).Width,Info(1).Height];
    num_px = imsize(1)*imsize(2);
    n_frames = length(Info);
    destpath = strcat(paths{1}, 'dF_', files{1});
    
    % Read timeseries
    inputim = zeros(imsize(2), imsize(1),n_frames);
    for frame=1:n_frames
        inputim(:,:,frame) = single(imread(fullfilename,frame));
    end
end
imsize = size(inputim);
n_px = imsize(1)*imsize(2);
ft_ms = inputdlg('Enter ft_s in ms:');
ft_s = str2num(ft_ms{1})/1000;

%% Subtract dark count or background

%% Bleaching correction

%% Classify as background
aip = mean(inputim,3);
graylims = [min(aip,[],'all') max(aip,[],'all')];
uint8_aip = (aip./graylims(2)).*255;

[c_abs,e] = histcounts(aip); % Count frequency of intensity values
c_rel = c_abs./sum(c_abs); % Relative frequency
c_cum = cumsum(c_rel); % Cummulative frequency
[~,idx] = min(abs(c_cum - params.bg_perc));
bg_thresh = e(idx);

background = true(size(aip,1),size(aip,2));
background(aip > bg_thresh) = false;

%% GUI
bg_fig = figure();
ctrl_fig = figure('Position', [750 150 300 40]);
sl = uicontrol('Style','slider','Min',graylims(1),'Max',graylims(2),...
    'SliderStep',[1 1]./graylims(2),'Value',bg_thresh, 'Position', [0 0 250 30], 'Callback', {@sl_cb, background, aip, uint8_aip, bg_fig});
pb = uicontrol('parent', ctrl_fig, 'style','pushbutton', 'Position', [255 0 40 30], 'String', 'OK', 'Callback', {@pb_cb});
figure(bg_fig);
colormap('gray');
imagesc(uint8(uint8_aip)); axis('off');
hold on, imagesc(uint8(background.*255), 'AlphaData',.4);
uiwait;

background(:,:) = true;
background(aip > bg_thresh) = false;

%% Calculate dF/F with rolling baseline
disp('Calculate dF over F stack...');
winsize = round(params.fbase_winsize_s/ft_s);
lintraces = reshape(inputim,[n_px n_frames]);
avtraces = lintraces;
avsteps = NaN(n_frames,2);
avsteps(1:winsize+1,1) = 1; avsteps(1:winsize+1,2) = winsize;
avsteps(winsize+2:n_frames,1) = 2:n_frames-winsize;
avsteps(winsize+2:n_frames,2) = winsize+1:n_frames-1;

for iAv = 1:n_frames
    avtraces(:,iAv) = mean(lintraces(:,avsteps(iAv,1):avsteps(iAv,2)),2);
end

dFtraces = lintraces-avtraces;
px_av_val = mean(lintraces,2);
dFoFtraces = dFtraces./px_av_val;

dFoFim = reshape(dFoFtraces,imsize);
clear('lintraces','avtraces','dFtraces', 'dFoFtraces', 'avsteps', 'px_av_val');
% for idx = 1:size(dFoFim,3)
%     imwrite(uint16(dFoFim(:,:,idx)), destpath, 'tiff', 'WriteMode', 'append');
% end

%% Detect responding px
disp('Detect responding pixel...');
% Responding px
threshold = prctile(dFoFim,params.perc_thresh,3);
response = dFoFim > params.perc_fac.*threshold;
% Eliminate background
mask1 = repmat(background, [1 1 imsize(3)]);
response(mask1) = false;

% for idx = 1:size(response,3)
%     imwrite(response(:,:,idx).*255, fullfile(paths{1},'responding_px.tif'), 'tiff', 'WriteMode', 'append');
% end

%% ******** Assign px to ROI ********
%% Component detection in 2D
disp('Find connected components...');
all_connc_cell=cell(n_frames,1);
tic;
parfor t =1:n_frames
    bw=response(:,:,t);
    connc = bwconncomp(bw);
    all_connc_cell{t,1} = connc.PixelIdxList;
end
toc

tic;
num_cc_px = 0;
for t =1:n_frames
    conn_px = cellfun(@numel, all_connc_cell{t,1}, 'UniformOutput', false);
    conn_px = cellfun(@(x) x, conn_px);
    connc_cell{t,1} = all_connc_cell{t,1}(conn_px > params.connc_px_thresh);
    num_cc_px = num_cc_px+sum(conn_px(conn_px > params.connc_px_thresh));
end
keep = ~cellfun(@isempty, connc_cell);
connc_cell = connc_cell(keep,:);
connc_frame_id = find(keep==1);
toc

tic;
roi_idx_mtrx = NaN(num_cc_px,4);
cnter = 0; cccter=1;
for t = 1:numel(connc_frame_id)
    ncc = numel(connc_cell{t,1});
    for cc = 1:ncc
        [r,c] = ind2sub([imsize(1) imsize(2)], connc_cell{t,1}{cc});
        n = numel(r);
        roi_idx_mtrx(cnter+1:cnter+n,:) = [repelem(cccter,n)' r c repelem(connc_frame_id(t),n)'];
        cnter=cnter+n; cccter=cccter+1;
    end
end
toc
tot_num_cc = cccter-1;

if tot_num_cc ~= 0
    %% Fill gaps
    disp('Fill gaps in connected regions...');
    tic;
    roi_idx_mtrx_filled = fill_roi_gaps(imsize, roi_idx_mtrx,[2 4], tot_num_cc,params.fill_thresh);
    toc
    
    % for idx = 1:size(response,3)
    %     imwrite(response_filled(:,:,idx).*255, fullfile(paths{1},'response_filled.tif'), 'tiff', 'WriteMode', 'append');
    % end
    
    %% Test overlap between ROIs and merge - in a loop until also combined, filled ROIs dont overlap
    disp('Check for overlapping ROIs...');
    % Symmetric matrix with relativ overlap
    overlap_cell = cell(1);
    roi_idx_cell = cell(1,2);
    roi_idx_cell{1,1} = roi_idx_mtrx_filled;
    roi_idx_cell{1,2} = roi_idx_mtrx_filled(:,1);
    n_rois = tot_num_cc;
    combine_more = true;
    rndcnt = 1;
    
    while combine_more
        disp('*New round*');
        tmp_idx_mtrx = roi_idx_cell{end,1};
        roi_ids = unique(tmp_idx_mtrx(:,1));
        tmp_num_cc = n_rois(end);
        overlap_mtrx = NaN(tmp_num_cc, tmp_num_cc);
        
        disp('Calculate overlap...');tic;
        for r = 1:tmp_num_cc
            tmpid1 = roi_ids(r);
            tmpcc1 = tmp_idx_mtrx(tmp_idx_mtrx(:,1) == tmpid1,2:3);
            overlap_mtrx(r,r) = 1;
            for c = r+1:tmp_num_cc
                tmpid2 = roi_ids(c);
                tmpcc2 = tmp_idx_mtrx(tmp_idx_mtrx(:,1) == tmpid2,2:3);
                cnt_overlap = 0;
                for ipx = 1:size(tmpcc2,1)
                    if any(all((tmpcc1 - tmpcc2(ipx,:))==0,2)), cnt_overlap = cnt_overlap+1; end
                end
                px_tot = min([size(tmpcc2,1) size(tmpcc1,1)]);
                overlap_mtrx(r,c) = cnt_overlap/px_tot;
            end
        end
        toc
        
        %% Based on overlap, combine ROIs
        overlap_comb = logical(overlap_mtrx > params.overlap_thresh);
        if any(overlap_comb & logical(overlap_mtrx < 1), 'all')
            disp('Combine ROIs...');tic;
            r = 1;
            while r <= size(overlap_comb,1)
                addrows = find(overlap_comb(r,:)==1);
                if ~isempty(addrows) && ~all(sum(overlap_comb(addrows(2:end),:),2) == 0)
                    addrows = addrows(2:end);
                    for addr = 1:numel(addrows), overlap_comb(r,:) = overlap_comb(r,:) | overlap_comb(addrows(addr),:); end
                    overlap_comb(addrows,:) = false;
                else
                    r = r+1;
                end
            end
            keep = ~all(~overlap_comb,2);
            new_comb = overlap_comb(keep,:);
            num_new_cc = size(new_comb,1);
            
            new_cc_id = NaN(size(tmp_idx_mtrx,1),1);
            ori_pter = roi_idx_cell{1,2};
            new_ori_cc_id = NaN(numel(ori_pter),1);
            for cc = 1:num_new_cc
                ids = find(new_comb(cc,:) == 1);
                change1 = tmp_idx_mtrx(:,1) == ids; change1 = logical(sum(change1,2));
                new_cc_id(change1) = cc;
                change2 = ori_pter == ids; change2 = logical(sum(change2,2));
                new_ori_cc_id(change2) = cc;
            end
            
            %% Discard duplicate px
            tmp_idx_mtrx_combined = [];
            for cc = 1:num_new_cc
                tmp_px = tmp_idx_mtrx(new_cc_id == cc,:);
                [sorted,sortidx] = sortrows(tmp_px, [2 3]);
                sorted = sorted(:,2:3);
                doubles = find(sum(abs(diff(sorted,1)),2) == 0);
                tmp_px(doubles,:) = [];
                tmp_idx_mtrx_combined = [tmp_idx_mtrx_combined; repelem(cc,size(tmp_px,1))', tmp_px];
            end
            
            %% Fill gaps second time
            tic;
            tmp_idx_mtrx_comb_filled = fill_roi_gaps(imsize, tmp_idx_mtrx_combined, [3 4], num_new_cc,params.fill_thresh);
            toc
            
            overlap_cell{rndcnt,1} = overlap_mtrx;
            overlap_cell{rndcnt,2} = new_comb;
            rndcnt = rndcnt+1;
            roi_idx_cell{rndcnt,1} = tmp_idx_mtrx_comb_filled;
            roi_idx_cell{rndcnt-1,2} = new_cc_id;
            roi_idx_cell{1,2} = new_ori_cc_id;
            n_rois(rndcnt) = num_new_cc;
        else
            %% Delete left overlapping parts
            for cc = 1:tmp_num_cc
                tmp_px = tmp_idx_mtrx(tmp_idx_mtrx(:,1) == cc,:);
                del=[];
                for ipx = 1:size(tmp_px,1)
                    del_test = sum(abs(tmp_idx_mtrx(:,2:3)-tmp_px(ipx,2:3)),2);
                    if any(logical(del_test==0) & logical(tmp_idx_mtrx(:,1)~=cc))
                        del = [del; find(del_test==0 & tmp_idx_mtrx(:,1)~=cc)];
                    end
                end
                tmp_idx_mtrx(del,:) = [];
            end
            roi_idx_cell{rndcnt,1} = tmp_idx_mtrx;
            roi_idx_cell{rndcnt,2} = tmp_idx_mtrx(:,1);
            
            roi_ids1= unique(roi_idx_cell{1,2});
            roi_ids2 = unique(tmp_idx_mtrx(:,1));
            pter = roi_idx_cell{rndcnt,1}(:,1);
            for iID = 1:numel(roi_ids1)
                if all(~(roi_ids1(iID) == roi_ids2))
                    roi_idx_cell{rndcnt,1}(pter == roi_ids1(iID)) = [];
                    roi_idx_cell{1,2}(pter == roi_ids1(iID),:) = [];
                    pter(pter == roi_ids1(iID),:) = [];
                end
            end
            combine_more = false;
        end
    end
    
    
    %% Select combinatino round
    rnd = size(roi_idx_cell,1);
    
    %% Create colormap
    c1=[1 1 0]; %G
    c2=[0 1 1]; %Y
    c3=[1 0 0]; %R
    % According to time:
    % n1=round(imsize(3)/10); n2=round(imsize(3)/10);
    % According to combination:
    n1=n_rois(rnd); n2=n1;
    cmap=[linspace(c1(1),c2(1),n1);linspace(c1(2),c2(2),n1);linspace(c1(3),c2(3),n1)];
    cmap(:,end+1:end+n2)=[linspace(c2(1),c3(1),n2);linspace(c2(2),c3(2),n2);linspace(c2(3),c3(3),n2)];
    
    %% Extract ROI boundaries - combined ROIs
    disp('Find ROI boundaries...');tic;
    num_cc = n_rois(rnd);
    bounds_cell_comb = cell(num_cc,2);
    cnt=1; cc=1;
    idx_mtrx = roi_idx_cell{rnd,1};
    while cc <= num_cc
        cc_px = idx_mtrx(idx_mtrx(:,1) == cc,2:3);
        bw = zeros(imsize(1),imsize(2),1);
        for ipx = 1:size(cc_px,1), bw(cc_px(ipx,1), cc_px(ipx,2)) = 1; end
        [B,~] = bwboundaries(bw,8,'noholes');  % Find boundary
        bounds_cell_comb{cnt,1} = B;
        % Color code according to z
        %     colidx = ceil(mean(cc_px(:,3))/10);
        %     bounds_cell_comb{cnt,2} = cmap(:,colidx);
        bounds_cell_comb{cnt,2} = cmap(:,cc);
        cnt = cnt+1;
        cc = cc+1;
    end
    bounds_cell_comb = bounds_cell_comb(1:cnt-1,:);
    
    %% Extract ROI boundaries - all ROIs
    cmap = get_colormap([1 0 0],[1 1 0],[0 1 1],n_rois(1));
    num_cc = n_rois(1);
    bounds_cell = cell(num_cc,3);
    cnt=1; cc=1;
    idx_mtrx = roi_idx_cell{1,1};
    cc_pter = roi_idx_cell{1,2};
    while cc <= num_cc
        cc_px = idx_mtrx(idx_mtrx(:,1) == cc,2:4);
        bw = zeros(imsize(1),imsize(2),1);
        for ipx = 1:size(cc_px,1), bw(cc_px(ipx,1), cc_px(ipx,2)) = 1; end
        [B,~] = bwboundaries(bw,8,'noholes');  % Find boundary
        bounds_cell{cnt,1} = B;
        %     colidx = ceil(mean(cc_px(:,3))/10); % Color code according to z
        colidx = unique(cc_pter(idx_mtrx(:,1)==cc)); % Color code according to component it belongs to
        bounds_cell{cnt,2} = cmap(:,colidx);
        bounds_cell{cnt,3} = colidx;
        cnt = cnt+1;
        cc = cc+1;
    end
    bounds_cell = bounds_cell(1:cnt-1,:);
    toc
else
    fprintf('**************\nNo ROIs found\n**************\n');
end

t_proc = toc(t_all);
fprintf('Total processing time: %i min\n', round(t_proc/60));

%% Local function
function sl_cb(hOb,~, background, aip, uint8_aip, bg_fig)
global bg_thresh
bg_thresh = round(hOb.Value);
background(:,:) = true;
background(aip > bg_thresh) = false;
figure(bg_fig);
clf;
colormap('gray');
imagesc(uint8(uint8_aip));axis('off');
hold on, imagesc(uint8(background.*255), 'AlphaData',.4);
end

function pb_cb(~,~)
uiresume;
figures = get(groot,'Children');
findfig = strncmp({figures(:).Name}, 'ROI',3);
close(figures(~findfig));
end
