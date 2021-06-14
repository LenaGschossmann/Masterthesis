function run_roi_detection(infostruct, params)
%% Function for automatic detection of ROI defined by intensity
% Detect ROIs and save as struct to load and process further

%% Add paths
wd = pwd(); p = split(wd,'\'); p = p(1:end-1); pp = []; for ii=1:numel(p), pp = fullfile(pp,p{ii}); end
addpath(genpath(pp)), clear('p','pp');

nfiles = size(infostruct,2);
% bleach_corr_win_s = 5;
t_proc = 600;
% dark_count = 0.528155;
rectype = questdlg('Select type of recording','Recording type','Confocal', '2PM', 'Confocal');

%% Set parameters dependend on recording type if not given as input
if strncmp(rectype,'C',1)
    if isempty(params)        
        params = get_predefined_params('connc_px_thresh_low', 15,'connc_px_thresh_high', 5,...
            'responding_px_base_low', 0.9, 'responding_px_base_high', 1.8,...
            'responding_px_factor_low', 2, 'responding_px_factor_high', 1, ...
            'responding_px_slope', [5 20],...
            'px_exclusion', 0.25, 'corr_threshold', 0.25, 'fill_thresh',0.25);
    end
else
    if isempty(params)
        params = get_predefined_params('connc_px_thresh_low', 16,'connc_px_thresh_high', 6,...
            'responding_px_base_low', 1.3, 'responding_px_base_high', 2.9,...
            'responding_px_factor_low', 6, 'responding_px_factor_high', 1, ...
            'px_exclusion', 0.25, 'corr_threshold', 0.25, 'fill_thresh',0.25);
        % Slope 20-15
        %         'responding_px_base_low', 1.6, 'responding_px_base_high', 2.5,...
        %             'responding_px_factor_low', 2, 'responding_px_factor_low', 3, ...
    end
end

%% Main loop
t_all = tic;
for iF = 1:nfiles
    roidata = struct();
    save_empty = true;
    
    fprintf('\n***************\nProcessing file %i of %i | Expected time remaining: %s min\n***************\n', iF,nfiles,num2str(round(nfiles*t_proc/60,2)));
    
    %% Unpack information
    filepointer = infostruct(iF).path;
    [filepath,filename, ~] = fileparts(filepointer);
    savepath = fullfile(filepath,'ROIFlow_Analysis');
    if ~isfolder(savepath), mkdir(savepath); end
    ft_s = infostruct(iF).frametime_s;
    %     bg_thresh = infostruct(iF).bg_threshold;
    bg_corr_val = infostruct(iF).bg_2pm_corr;
    mocorr = infostruct(iF).mocorr;
    shift = params.dFoF_baseline_shift;
    %     sd_factor = params.responding_px_factor;
    %     bleachcorr = infostruct(iF).bleachcorr;
    xy_scale = infostruct(iF).xy_scale;
    winsize = round(params.dFoF_baseline_s/ft_s);
    if winsize < params.dFoF_baseline_frames, winsize = params.dFoF_baseline_frames; end
    if winsize > 20, winsize = 20; end
    if strncmp(rectype,'C',1)
        pmt = NaN;
        offset = NaN;
    else
        pmt = infostruct(iF).pmt;
        offset = infostruct(iF).offset;
    end
    
    %% Import file (& Motion Correction)
    if mocorr
        disp('Run motion correction...');
        inputim = run_mocorr(filepointer, true, true);
        %         inputim = inputim{1};
    else
        [inputim] = read_timeseries(filepointer, 'load');
    end
    inputim = double(inputim);
    
    imsize = size(inputim);
    n_px = imsize(1)*imsize(2);
    n_frames = imsize(3);
    lintraces = reshape(inputim,[n_px n_frames]);
    
    %% Detect background px
    sd_alltrc = std(lintraces,[],2);
    sd_range = [prctile(sd_alltrc,5,'all') prctile(sd_alltrc,95,'all')];
    sd_thresh = sd_range(1) + diff(sd_range) * params.px_exclusion;
    exclude_px = sd_alltrc < sd_thresh;
    bg_thresh = mean(lintraces(exclude_px,:), 'all');
    background = exclude_px;
    
    %     avtr = mean(lintraces,2);
    %     f=figure();
    %     subplot(1,2,1), imshow(reshape(uint8(exclude_px.*255),[imsize(1) imsize(2)]));
    %     subplot(1,2,2), imshow(reshape(uint16(avtr.*5),[imsize(1) imsize(2)]));
    %     tmpp = strcat(filepath, '\resp_px_crit1_', infostruct(iF).name, '.png');
    %     saveas(f, tmpp);
    %     close(f);
    
    %% Subtract dark count or background
    if strncmp(rectype, '2',1)
        subtract_value = - bg_corr_val; % Negative sign necessary, if difference had been calculated as 0-offset - rec-specific offset
    else
        subtract_value = prctile(lintraces(exclude_px,:),5, 'all');
    end
    
    %     bgvalues = lintraces(background==1);
    %     [c_abs,e] = histcounts(bgvalues); % Count frequency of intensity values
    %     c_rel = c_abs./sum(c_abs); % Relative frequency
    %     c_cum = cumsum(c_rel); % Cummulative frequency
    %     [~,idx] = min(abs(c_cum - params.subtract_perc));
    %     subtract_value = e(idx);
    
    corr_traces = lintraces - subtract_value;
    
    
    %% Calculate area information
    fov_area = n_px*xy_scale*xy_scale;
    n_px_bg = sum(background);
    bg_area = n_px_bg*xy_scale*xy_scale;
    dend_area = fov_area-bg_area;
    
    %% Calculate dFoF stack
    disp('Calculate dF/F stack...');
    [~, FoFtraces, dFoFtraces] = rollBase_dFoF(corr_traces, winsize,shift, 'px');
    
    %     [~, ~, dFoFtraces] = rollBase_dFoF(corr_traces, winsize,shift, 'roll');
    tmpim = reshape(FoFtraces,[imsize(1) imsize(2) n_frames]);
    tmpp = strcat(filepath, '\FoF_', infostruct(iF).name, '.tif');
    for idx = 1:n_frames
        imwrite(uint8(tmpim(:,:,idx).*30), tmpp, 'tiff', 'WriteMode', 'append');
    end
    
    %% Detect responding px
    disp('Detect responding pixel...');
    std_noise = NaN(n_px,1);
    for ipx = 1:n_px
        trc_neg = FoFtraces(ipx,:); trc_neg = trc_neg(trc_neg < 1);
        trc_synt_symm = [trc_neg 1+(1-trc_neg)];
        std_noise(ipx,1) = std(trc_synt_symm);
    end
    
    % * Adjustible threshold *
%     test = FoFtraces-1;
%     test = reshape(test, [n_px n_frames]);
%     test(exclude_px,:) = 0;
%     testfof = reshape(test,[imsize(1) imsize(2) n_frames]);
%     testsd = reshape(std_noise,[imsize(1) imsize(2)]);
%     testfof = (testfof)./testsd;
%     pos = testfof(testfof>0);
%     cnttot = numel(pos);
%     thetas = 0:0.1:4;
%     results = NaN(numel(thetas),1);
%     for ii=1:numel(thetas)
%         tmptheta = thetas(ii);
%         cnt=sum(pos>=tmptheta);
%         results(ii) = 100*cnt/cnttot;
%     end
%     
%     slope = 0.15/(thetas(find(results<=params.responding_px_slope(2),1,'first'))-thetas(find(results<=params.responding_px_slope(1),1,'first')));
%     sd_factor_low = params.responding_px_base_low - params.responding_px_factor_low * slope;
%     sd_factor_high = params.responding_px_base_high - params.responding_px_factor_high * slope;
    % **
    sd_factor_low = 1.2;
    sd_factor_high = 1.5;
    
    resp_thresh = 1+ std_noise .* sd_factor_low;
    response_low = FoFtraces > resp_thresh;
    response_low(exclude_px,:) = 0;
    response_low = reshape(response_low,[imsize(1) imsize(2) n_frames]);
    
    resp_thresh2 = 1+ std_noise .* sd_factor_high;
    response_high = FoFtraces > resp_thresh2;
    response_high(exclude_px,:) = 0;
    response_high = reshape(response_high,[imsize(1) imsize(2) n_frames]);
    
    test = zeros(size(response_low));
    
    %     tmpp = strcat(filepath, '\response_low_', infostruct(iF).name, '.tif');
    %     for idx = 1:n_frames
    %         imwrite(uint8(response_low(:,:,idx).*255), tmpp, 'tiff', 'WriteMode', 'append');
    %     end
    %
    %     tmpp = strcat(filepath, '\response_high_', infostruct(iF).name, '.tif');
    %     for idx = 1:n_frames
    %         imwrite(uint8(response_high(:,:,idx).*255), tmpp, 'tiff', 'WriteMode', 'append');
    %     end
    
    %% ******** Assign px to ROI ********
    %% Component detection in 2D
    disp('Find connected components...');
    % Detect connected pixel (connectivity: 8) in each frame and save as linear index
    % list (2D array indices!)
    all_connc_cell_low=cell(n_frames,1);
    all_connc_cell_high=cell(n_frames,1);
    tic;
    for t = 1:winsize
        all_connc_cell_low{t,1} = [];
        all_connc_cell_high{t,1} = [];
    end
    parfor t = winsize+1:n_frames
        bw=response_low(:,:,t);
        connc = bwconncomp(bw,4);
        all_connc_cell_low{t,1} = connc.PixelIdxList;
        
        bw=response_high(:,:,t);
        connc2 = bwconncomp(bw,4);
        all_connc_cell_high{t,1} = connc2.PixelIdxList;
    end
    toc
    
    % Discard all components with less than minimum number of connected
    % pixel (preset threshold)
    disp('Discard components with too few px...');
    tic;
    num_cc_px = 0;
    connc_cell_low = cell(n_frames,1);
    connc_cell_high = cell(n_frames,1);
    for t = winsize+1:n_frames
        conn_px_low = cellfun(@numel, all_connc_cell_low{t,1}, 'UniformOutput', false);
        conn_px_low = cellfun(@(x) x, conn_px_low);
        connc_cell_low{t,1} = all_connc_cell_low{t,1}(conn_px_low > params.connc_px_thresh_low);
        %         num_cc_px = num_cc_px+sum(conn_px_low(conn_px_low > params.connc_px_thresh_low));
        
        conn_px_high = cellfun(@numel, all_connc_cell_high{t,1}, 'UniformOutput', false);
        conn_px_high = cellfun(@(x) x, conn_px_high);
        connc_cell_high{t,1} = all_connc_cell_high{t,1}(conn_px_high > params.connc_px_thresh_high);
    end
    keep = ~cellfun(@isempty, connc_cell_low);
    connc_cell_low = connc_cell_low(keep,:);
    connc_frame_id = find(keep==1); % Indices of frames with connected components
    
    connc_cell_high = connc_cell_high(keep,:);
    toc
    
    % Discard connected components if not detected by high and low bar
    % detection run
    tic;
    del = [];
    linidx = cell(size(connc_cell_low,1),1);
    for idx = 1:size(connc_cell_low,1)
        for idx2 = 1:numel(connc_cell_high{idx,1})
            linidx{idx} = [linidx{idx}; connc_cell_high{idx,1}{idx2}];
        end
    end
    
    for idx = 1:size(connc_cell_low,1)
        keep = false(numel(connc_cell_low{idx,1}),1);
        tmpn = numel(connc_cell_low{idx,1});
        for idx2 = 1:tmpn % loop through components on this frame
            %             [r,c] = ind2sub([imsize(1) imsize(2)], connc_cell_low{idx,1}{idx2});
            for idx3 = 1:numel(connc_cell_low{idx,1}{idx2})
                if any(connc_cell_low{idx,1}{idx2}(idx3) == linidx{idx}), keep(idx2) = true; end
            end
            %                 if sum(response_high(connc_cell_low{idx,1}{idx2},t)) >= round(tmpn/2), keep(idx2) = true; end
        end
        if all(~keep), del = [del idx];
        else
            connc_cell_low{idx,1} = connc_cell_low{idx,1}(keep);
            num_cc_px = num_cc_px+sum(cellfun(@numel, connc_cell_low{idx,1}));
        end
    end
    connc_cell_low(del) = [];
    connc_frame_id(del) = []; % Indices of frames with connected components
    toc
    
    % Create matrix with information of all px within connected components
    disp('Create master matrix...');
    tic;
    roi_idx_mtrx = NaN(num_cc_px,4);
    % Master matrix with px information: col1: conn.component id, col2: row idx,
    % col3: column idx, col4: frame of conn.component
    cnter = 0; cccter=1;
    for t = 1:numel(connc_frame_id)
        ncc = numel(connc_cell_low{t,1}); % number connected components in this frame
        for cc = 1:ncc
            [r,c] = ind2sub([imsize(1) imsize(2)], connc_cell_low{t,1}{cc});
            n = numel(r); % Number of px in component
            roi_idx_mtrx(cnter+1:cnter+n,:) = [repelem(cccter,n)' r c repelem(connc_frame_id(t),n)'];
            
            for ii = 1:size(r,1), test(r(ii),c(ii),connc_frame_id(t)) = 1; end
            
            cnter=cnter+n; cccter=cccter+1;
        end
    end
    toc
    tot_num_cc = cccter-1;
    
    %     tmpp = strcat(filepath, '\bin_components_single_', infostruct(iF).name, '.tif');
    %     for idx = 1:n_frames
    %         imwrite(uint8(test(:,:,idx).*255), tmpp, 'tiff', 'WriteMode', 'append');
    %     end
    
    if tot_num_cc ~= 0
        roi_idx_mtrx_filled = roi_idx_mtrx;
        
        %% Combine ROIs
        disp('Combine ROIs...');
        % Symmetric matrix with relativ overlap
        roi_idx_cell = cell(1,2);
        roi_idx_cell{1,1} = roi_idx_mtrx_filled; % Master matrix of all px and their assigned components - round 1
        roi_idx_cell{1,2} = roi_idx_mtrx_filled(:,1); % Component assignment of each px
        n_rois = numel(unique(roi_idx_mtrx_filled(:,1)));
        
        tmp_idx_mtrx = roi_idx_cell{1,1};
        roi_ids = unique(tmp_idx_mtrx(:,1));
        tmp_num_cc = n_rois(end);
        expand8_yx = [-1 -1; -1 0; -1 1; 0 -1; 0 0; 0 1; 1 -1; 1 0; 1 1];
        
        % Find 2nd order components
        test_im = zeros(imsize(1), imsize(2));
        for iCC = 1:tmp_num_cc
            cc_coords = tmp_idx_mtrx(tmp_idx_mtrx(:,1) == roi_ids(iCC),2:3);
            for iPx = 1:size(cc_coords,1)
                test_im(cc_coords(iPx,1), cc_coords(iPx,2)) = 1;
            end
        end
        conn_comp = bwconncomp(test_im,4);
        
        % Find components contributing to 2nd order component
        new_mapping = NaN(size(tmp_idx_mtrx,1),2);
        new_mapping(:,1) = tmp_idx_mtrx(:,1);
        for iCC = 1:tmp_num_cc
            tmppter = find(tmp_idx_mtrx(:,1) == roi_ids(iCC));
            tmplinidx = sub2ind([imsize(1) imsize(2)], tmp_idx_mtrx(tmppter(1),2), tmp_idx_mtrx(tmppter(1),3));
            for iCCC = 1:size(conn_comp.PixelIdxList,2)
                if any(conn_comp.PixelIdxList{iCCC} == tmplinidx)
                    new_mapping(tmppter,2) = iCCC;
                end
            end
        end
        
        % Check correlation of 1st order components
        for iCCC = 1:size(conn_comp.PixelIdxList,2)
            tmp_cc = unique(new_mapping(new_mapping(:,2) == iCCC,1));
            tmp_n_cc = numel(tmp_cc);
            tmp_corr_mtrx = zeros(tmp_n_cc, tmp_n_cc);
            for i = 1:tmp_n_cc, tmp_corr_mtrx(i,i) = 1; end
            
            %             if tmp_n_cc>1
            %                 figure(); spcter=0;
            %                 for r = 1:tmp_n_cc-1, for c = r+1:tmp_n_cc, spcter = spcter+1; end, end
            %                 spc = 2; spr = ceil(spcter/spc); spcter2=1;
            %             end
            for r = 1:tmp_n_cc-1
                tmpid1 = tmp_cc(r);
                tmpcc1 = tmp_idx_mtrx(tmp_idx_mtrx(:,1) == tmpid1,2:3); % y x
                
                for c = r+1:tmp_n_cc
                    tmpid2 = tmp_cc(c);
                    tmpcc2 = tmp_idx_mtrx(tmp_idx_mtrx(:,1) == tmpid2,2:3);
                    ovlap_cc2 = false(size(tmpcc2,1),1);
                    for ipx = 1:size(tmpcc2,1)
                        if any(all((tmpcc1 - tmpcc2(ipx,:))==0,2))
                            ovlap_cc2(ipx) = true;
                        end
                    end
                    cnt_overlap = sum(ovlap_cc2);
                    
                    test_cc1 = tmpcc1; test_cc2 = tmpcc2;
                    if cnt_overlap > 0
                        % Expand overlapping region
                        ovlap_info = tmpcc2(ovlap_cc2,:);
                        exp_px = [];
                        
                        for ipx = 1:size(ovlap_info,1)
                            exp_px = [exp_px; ovlap_info(ipx,:) + expand8_yx];
                        end
                        
                        for ipx = 1:size(exp_px,1)
                            del = find(all((exp_px - exp_px(ipx,:))==0,2) == 1);
                            if numel(del) > 1, exp_px(del(2:end),:) = NaN; end
                        end
                        exp_px(isnan(exp_px(:,1)),:) = [];
                        
                        for ipx  = 1:size(exp_px,1)
                            del = find(all((test_cc1 - exp_px(ipx,:))==0,2) == 1);
                            test_cc1(del,:) = [];
                            del = find(all((test_cc2 - exp_px(ipx,:))==0,2) == 1);
                            test_cc2(del,:) = [];
                        end
                    end
                    
                    if ~isempty(test_cc1) && ~isempty(test_cc2)
                        lindices1 = sub2ind([imsize(1) imsize(2)], test_cc1(:,1), test_cc1(:,2));
                        avtrc1 = mean(FoFtraces(lindices1,:),1);
                        lindices2 = sub2ind([imsize(1) imsize(2)], test_cc2(:,1), test_cc2(:,2));
                        avtrc2 = mean(FoFtraces(lindices2,:),1);
                        av1 = mean(avtrc1); av2 = mean(avtrc2);
                        tmp_corr_mtrx(r,c) = sum((av1-avtrc1).*(av2-avtrc2)) / sqrt(sum((av1-avtrc1).^2)*sum((av2-avtrc2).^2));
                        %                         subplot(spr,spc,spcter2), plot(avtrc1), hold on, plot(avtrc2,'r'), title(strcat('r= ',num2str(tmp_corr_mtrx(r,c)), ',px b/r:', num2str(size(lindices1,1)),',',num2str(size(lindices2,1))))
                        %                         spcter2 = spcter2+1;
                    else
                        tmp_corr_mtrx(r,c) = 1;
                    end
                end
                
            end
            tmp_corr_mtrx = tmp_corr_mtrx + rot90(fliplr(tmp_corr_mtrx));
            cnter = max(new_mapping(:,2))+1;
            for r = 1:tmp_n_cc
                if all(tmp_corr_mtrx(r,:) <= params.corr_threshold | tmp_corr_mtrx(r,:) == 2)
                    new_mapping(new_mapping(:,1) == tmp_cc(r),2) = cnter;
                    cnter = cnter+1;
                end
            end
        end
        
        % Combine
        rndcnt = 2;
        tmp_idx_mtrx(:,1) = new_mapping(:,2);
        roi_idx_cell{rndcnt,1} = tmp_idx_mtrx;
        roi_idx_cell{rndcnt,2} = new_mapping(:,2);
        roi_idx_cell{1,2} = new_mapping(:,2);
        n_rois(rndcnt) = numel(unique(tmp_idx_mtrx(:,1)));
        
        
        %% Delete left overlapping parts
        ccs = unique(tmp_idx_mtrx(:,1));
        for cc1 = 1:numel(ccs)-1
            tmpidx1 = tmp_idx_mtrx(:,1) == ccs(cc1);
            tmpcc1 = tmp_idx_mtrx(tmpidx1,2:3);
            if ~isempty(tmpcc1)
                for cc2 = cc1+1:numel(ccs)
                    tmpidx2 = tmp_idx_mtrx(:,1) == ccs(cc2);
                    tmpcc2 = tmp_idx_mtrx(tmpidx2,2:3);
                    
                    if ~isempty(tmpcc2)
                        ovlap_cc2 = false(size(tmpcc2,1),1);
                        ovlap_cc1 = false(size(tmpcc1,1),1);
                        for ipx = 1:size(tmpcc2,1)
                            if any(all((tmpcc1 - tmpcc2(ipx,:))==0,2))
                                ovlap_cc2(ipx) = true;
                                ovlap_cc1(all((tmpcc1 - tmpcc2(ipx,:))==0,2)) = true;
                            end
                        end
                        
                        if sum(ovlap_cc1) > 0
                            % Expand overlapping region
                            ovlap_info = tmpcc2(ovlap_cc2,:);
                            exp_px = [];
                            test_cc1 = tmpcc1; test_cc2 = tmpcc2;
                            
                            for ipx = 1:size(ovlap_info,1)
                                exp_px = [exp_px; ovlap_info(ipx,:) + expand8_yx];
                            end
                            
                            for ipx = 1:size(exp_px,1)
                                del = find(all((exp_px - exp_px(ipx,:))==0,2) == 1);
                                if numel(del) > 1, exp_px(del(2:end),:) = NaN; end
                            end
                            exp_px(isnan(exp_px(:,1)),:) = [];
                            
                            for ipx  = 1:size(exp_px,1)
                                del = find(all((test_cc1 - exp_px(ipx,:))==0,2) == 1);
                                test_cc1(del,:) = [];
                                del = find(all((test_cc2 - exp_px(ipx,:))==0,2) == 1);
                                test_cc2(del,:) = [];
                            end
                            
                            if ~isempty(test_cc1) && ~isempty(test_cc2)
                                lindices1 = sub2ind([imsize(1) imsize(2)], test_cc1(:,1), test_cc1(:,2));
                                avtrc1 = mean(dFoFtraces(lindices1,:),1);
                                lindices2 = sub2ind([imsize(1) imsize(2)], test_cc2(:,1), test_cc2(:,2));
                                avtrc2 = mean(dFoFtraces(lindices2,:),1);
                                lindices_ovlap = sub2ind([imsize(1) imsize(2)], ovlap_info(:,1), ovlap_info(:,2));
                                avtrc_ovlap = mean(dFoFtraces(lindices_ovlap,:),1);
                                av1 = mean(avtrc1); av2 = mean(avtrc2); av_ovlap = mean(avtrc_ovlap);
                                r1 = sum((av1-avtrc1).*(av_ovlap-avtrc_ovlap)) / sqrt(sum((av1-avtrc1).^2)*sum((av_ovlap-avtrc_ovlap).^2));
                                r2 = sum((av2-avtrc2).*(av_ovlap-avtrc_ovlap)) / sqrt(sum((av2-avtrc2).^2)*sum((av_ovlap-avtrc_ovlap).^2));
                                
                                if r1 > r2
                                    tmpidx2_2 = find(tmpidx2==1); tmp_idx_mtrx(tmpidx2_2(ovlap_cc2),:) = NaN;
                                else, tmpidx1_2 = find(tmpidx1==1); tmp_idx_mtrx(tmpidx1_2(ovlap_cc1),:) = NaN;
                                end
                            elseif isempty(test_cc1) && isempty(test_cc2)
                                tmp_idx_mtrx(tmpidx1,:) = NaN;
                                tmp_idx_mtrx(tmpidx2,:) = NaN;
                            elseif isempty(test_cc1), tmp_idx_mtrx(tmpidx1,:) = NaN;
                            elseif isempty(test_cc2), tmp_idx_mtrx(tmpidx2,:) = NaN;
                            end
                        end
                    end
                end
            end
            tmp_idx_mtrx(isnan(tmp_idx_mtrx(:,1)),:) = [];
        end
        
        % Delete all components that are too small
        ccs = unique(tmp_idx_mtrx(:,1));
        del = false(size(tmp_idx_mtrx,1),1);
        for cc = 1:numel(ccs)
            if sum(tmp_idx_mtrx(:,1) == ccs(cc)) <= params.connc_px_thresh_low
                del(tmp_idx_mtrx(:,1) == ccs(cc),:) = true;
            end
        end
        tmp_idx_mtrx(del,:) = [];
        
        % Fill gaps for the last time
        if ~isempty(tmp_idx_mtrx)
            tic;
            tmp_idx_mtrx = fill_roi_gaps(imsize, tmp_idx_mtrx, [2 3],params.fill_thresh);
            toc
        end
        
        % Delete all components that are too small
        ccs = unique(tmp_idx_mtrx(:,1));
        del = false(size(tmp_idx_mtrx,1),1);
        for cc = 1:numel(ccs)
            if sum(tmp_idx_mtrx(:,1) == ccs(cc)) <= params.connc_px_thresh_low
                del(tmp_idx_mtrx(:,1) == ccs(cc),:) = true;
            end
        end
        tmp_idx_mtrx(del,:) = [];
        
        if ~isempty(tmp_idx_mtrx)
            roi_idx_cell{rndcnt,1} = tmp_idx_mtrx;
            roi_idx_cell{rndcnt,2} = tmp_idx_mtrx(:,1);
            
            roi_ids1= unique(roi_idx_cell{1,2});
            roi_ids2 = unique(roi_idx_cell{end,2});
            pter = roi_idx_cell{1,2};
            for iID = 1:numel(roi_ids1)
                if all(~(roi_ids1(iID) == roi_ids2))
                    roi_idx_cell{1,2}(pter == roi_ids1(iID)) = NaN;
                end
            end
            % Sort ROIs new
            new_roi_ids = NaN(size(tmp_idx_mtrx,1),1);
            new_ori_ids =  NaN(size(roi_idx_cell{1,2},1),1);
            for iNum = 1:numel(roi_ids2)
                new_roi_ids(tmp_idx_mtrx(:,1) == roi_ids2(iNum)) = iNum;
                new_ori_ids(roi_idx_cell{1,2}(:,1) == roi_ids2(iNum)) = iNum;
            end
            
            roi_idx_cell{end,1}(:,1) = new_roi_ids;
            roi_idx_cell{1,2} = new_ori_ids;
            
            n_rois(rndcnt) = numel(unique(roi_idx_cell{rndcnt,2}));
            
            
            
            %% Extract ROI area and centroids
            disp('Extract ROI area and centroids...');
            tic;
            num_cc = n_rois(end);
            idx_mtrx = roi_idx_cell{end,1};
            roi_area_um = NaN(num_cc,1);
            roi_xy_center_um = NaN(num_cc,2);
            for cc = 1:num_cc
                cc_px = idx_mtrx(idx_mtrx(:,1) == cc,2:3);
                roi_area_um(cc,1) = size(cc_px,1)*xy_scale*xy_scale;
                roi_xy_center_um(cc,1) = mean(cc_px(:,2))*xy_scale;
                roi_xy_center_um(cc,2) = mean(cc_px(:,1))*xy_scale;
            end
            toc
            
            %% Sort ROIs by distance to origin
            dist_info = sqrt(roi_xy_center_um(:,1).^2 + roi_xy_center_um(:,2).^2);
            [~, sorted_id] = sort(dist_info, 'ascend');
            % Apply sorting
            roi_area_um = roi_area_um(sorted_id,:);
            roi_xy_center_um = roi_xy_center_um(sorted_id,:);
            
            sorted_roi_info = roi_idx_cell{end,1};
            new_id = NaN(size(sorted_roi_info,1),1);
            sorted_ori_id = NaN(size(roi_idx_cell{1,2},1),1);
            for cc = 1:num_cc
                new_id(sorted_roi_info(:,1) == sorted_id(cc),1) = cc;
                sorted_ori_id(roi_idx_cell{1,2}==sorted_id(cc)) = cc;
            end
            sorted_roi_info(:,1) = new_id; clear('new_id');
            [~, sorted] = sort(sorted_roi_info(:,1),'ascend');
            sorted_roi_info = sorted_roi_info(sorted,:);
            
            %% Extract ROI boundaries - combined ROIs
            disp('Find ROI boundaries...');tic;
            num_cc = n_rois(end);
            bounds_cell_comb = cell(num_cc,2);
            cnt=1; cc=1;
            idx_mtrx = sorted_roi_info;
            while cc <= num_cc
                cc_px = idx_mtrx(idx_mtrx(:,1) == cc,2:3);
                bw = zeros(imsize(1),imsize(2),1);
                for ipx = 1:size(cc_px,1), bw(cc_px(ipx,1), cc_px(ipx,2)) = 1; end
                [B,~] = bwboundaries(bw,8,'noholes');  % Find boundary
                bounds_cell_comb{cnt,1} = B;
                % Color as indicator of number of individual components
                bounds_cell_comb{cnt,2} = numel(unique(roi_idx_cell{1,1}(sorted_ori_id == cc,1)));
                %             bounds_cell_comb{cnt,2} = cc;
                cnt = cnt+1;
                cc = cc+1;
            end
            bounds_cell_comb = bounds_cell_comb(1:cnt-1,:);
            toc
            
            %% Extract traces
            % dF/F calculated using rolling average for F_b and F_0
            % F/F calculated as corrected raw intensity normalized by average
            % across timepoints
            disp('Extract average ROI traces...');
            tic;
            num_cc = n_rois(end);
            roi_av_traces = NaN(num_cc,n_frames);
            roi_av_FoF_traces = NaN(num_cc,n_frames);
            %         roi_av_dFoF_traces = NaN(num_cc,n_frames);
            mean_traces = NaN(num_cc,n_frames);
            idx_mtrx = sorted_roi_info;
            for cc = 1:num_cc
                cc_px = idx_mtrx(idx_mtrx(:,1) == cc,2:3);
                px_lin = sub2ind([imsize(1) imsize(2)],cc_px(:,1),cc_px(:,2));
                roi_av_traces(cc,:) = mean(corr_traces(px_lin,:),1);
                %             roi_av_dFoF_traces(cc,:) = mean(dFoFtraces(px_lin,:),1);
                av = mean(roi_av_traces(cc,:));
                roi_av_FoF_traces(cc,:) = roi_av_traces(cc,:)./ av;
                mean_traces(cc,:) = mean(corr_traces(px_lin,:),1);
            end
            [~, ~,roi_av_dFoF_traces] = rollBase_dFoF(mean_traces, winsize,shift, 'roll');
            toc
            
            %% AIP & MIP Projection
            corr_aip = reshape(mean(corr_traces,2),[imsize(1) imsize(2)]);
            tmpmin = min(corr_aip,[],'all');
            uint16_aip = corr_aip-tmpmin; tmpmax = max(uint16_aip,[],'all');
            uint16_aip = (uint16_aip./tmpmax).*(2^8);
            corr_mip = reshape(max(corr_traces,[],2),[imsize(1) imsize(2)]);
            tmpmin = min(corr_mip,[],'all');
            uint16_mip = corr_mip-tmpmin; tmpmax = max(uint16_mip,[],'all');
            uint16_mip = (uint16_mip./tmpmax).*(2^8);
            clear('tmpmin','tmpmax');
            
            %         % Save overview of all ROIs using AIP
            scrsz = get(0,'ScreenSize');
            red=[1 0 1]; green=[0 1 1]; blue=[0 .7 1];
            cmap_n = max(cellfun(@max, bounds_cell_comb(:,2)));
            cmap = get_colormap(red,green,blue,cmap_n);
            %         disp('Save overview...');
            %         ovwfig_aip = figure('Position', [0 0 scrsz(3:4)], 'PaperSize', scrsz(3:4));
            %         imagesc(uint16_aip);colormap('gray');axis('off');daspect([1 1 1]);
            %         for cc = 1:num_cc
            %             B = bounds_cell_comb{cc,1};
            %             hold on
            %             for k = 1:length(B)
            %                 boundary = B{k};
            %                 plot(boundary(:,2), boundary(:,1), 'Color',cmap(:,bounds_cell_comb{cc,2})', 'LineWidth', 1)
            %             end
            %         end
            %         outputname = strcat(savepath,'\ROIs_AIP_',filename,'.png');
            %         if isa(outputname,'cell'), outputname= outputname{1}; end
            %         saveas(ovwfig_aip, outputname);
            %         pause(0.3);
            %         close gcf;
            
            % Save overview of all ROIs using MIP
            ovwfig_mip = figure('Position', [0 0 scrsz(3:4)], 'PaperSize', scrsz(3:4));
            imagesc(uint16_mip);colormap('gray');axis('off');daspect([1 1 1]);
            for cc = 1:num_cc
                B = bounds_cell_comb{cc,1};
                hold on
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:,2), boundary(:,1), 'Color',cmap(:,bounds_cell_comb{cc,2})', 'LineWidth', 1)
                end
            end
            outputname = strcat(savepath,'\ROIs_MIP_',filename,'.png');
            if isa(outputname,'cell'), outputname= outputname{1}; end
            saveas(ovwfig_mip, outputname);
            pause(0.3);
            close gcf;
            
            %         ccmask = false(imsize(1), imsize(2));
            %         for ipx = 1:size(sorted_roi_info,1)
            %             ccmask(sorted_roi_info(ipx,2), sorted_roi_info(ipx,3)) = true;
            %         end
            %         figure, imagesc(ccmask*255);
            
            %% Pack information and save
            disp('Save information...');
            roidata.filepath = filepointer;
            roidata.frametime_s = ft_s;
            roidata.xy_scale = xy_scale;
            roidata.baseline_frames = winsize;
            roidata.bg_threshold = bg_thresh;
            roidata.subtract_value = subtract_value;
            roidata.offset = offset;
            roidata.pmt = pmt;
            roidata.params = params;
            roidata.mocorr_run = mocorr;
            roidata.aip = corr_aip;
            roidata.mip = corr_mip;
            roidata.background_mask = background;
            roidata.response_threshold = resp_thresh;
            roidata.dendritic_area_um = dend_area;
            roidata.background_area_um = bg_area;
            roidata.fov_area_um = fov_area;
            roidata.roi_bounds = bounds_cell_comb;
            roidata.roi_idx = roi_idx_cell{end,1};
            roidata.original_roi_idx = roi_idx_cell{1,1};
            roidata.n_rois = n_rois(end);
            roidata.n_original_rois = n_rois(1);
            roidata.dFoF_traces = roi_av_dFoF_traces;
            roidata.FoF_traces = roi_av_FoF_traces;
            roidata.traces = roi_av_traces;
            roidata.roi_centroids_um = roi_xy_center_um;
            roidata.roi_area_um = roi_area_um;
            
            outputname = strcat(savepath,'\ROIdata_',filename,'.mat');
            if isa(outputname,'cell'), outputname= outputname{1}; end
            save(outputname,'roidata', '-v7.3');
            clear('roidata', 'corr_mip','corr_aip','filepointer','mocorr','bg_thresh',...
                'ft_s','resp_thresh','background','bounds_cell_comb','roi_idx_cell');
            save_empty = false;
        end
    end
    if save_empty
        fprintf('**************\nNo ROIs found\n**************\n');
        
        corr_aip = reshape(mean(corr_traces,2),[imsize(1) imsize(2)]);
        corr_mip = reshape(max(corr_traces,[],2),[imsize(1) imsize(2)]);
        
        roidata.filepath = filepointer;
        roidata.frametime_s = ft_s;
        roidata.xy_scale = xy_scale;
        roidata.baseline_frames = winsize;
        roidata.bg_threshold = bg_thresh;
        roidata.subtract_value = subtract_value;
        roidata.params = params;
        roidata.mocorr_run = mocorr;
        roidata.offset = offset;
        roidata.pmt = pmt;
        roidata.aip = corr_aip;
        roidata.mip = corr_mip;
        roidata.background_mask = background;
        roidata.response_threshold = resp_thresh;
        roidata.dendritic_area_um = dend_area;
        roidata.background_area_um = bg_area;
        roidata.fov_area = fov_area;
        roidata.n_rois = 0;
        roidata.roi_bounds = [];
        roidata.roi_idx = [];
        roidata.original_roi_idx = [];
        roidata.n_original_rois = [];
        roidata.dFoF_traces = [];
        roidata.FoF_traces = [];
        roidata.traces = [];
        roidata.roi_centroids_um = [];
        roidata.roi_area_um = [];
        outputname = strcat(savepath,'\0_ROIDATA_',filename,'.mat');
        if isa(outputname,'cell'), outputname= outputname{1}; end
        save(outputname,'roidata');
    end
    
        tmpp = strcat(filepath, '\bin_components_', infostruct(iF).name, '.tif');
        for idx = 1:n_frames
            imwrite(uint8(test(:,:,idx).*255), tmpp, 'tiff', 'WriteMode', 'append');
        end
    
    clear('all_connc_cell_low', 'all_connc_cell_high', 'background', 'bg_area', 'bg_thresh',...
        'bleachcorr', 'ccs', 'ccter', 'cnter',...
        'conn_px_low', 'conn_px_high', 'connc_cell_low', 'connc_cell_high', 'connc_frame_id',...
        'corr_aip', 'corr_mip', 'corr_traces','rolldFoFtraces',...
        'dend_area', 'dFoFtraces', 'fov_area', 'ft_s', 'ii', 'imsize', 'inputim',...
        'keep', 'lintraces', 'mpcprr', 'n_frames', 'n_px', 'n_px_bg', 'num_cc_px', 'resp_thresh',...
        'response_low', 'response_high', 'roi_idx_mtrx', 'roidata', 'subtract_value', 't',...
        'tot_num_cc', 'xy_scale');
end
t_proc = toc(t_all);
fprintf('Total processing time: %s min\n', num2str(round(t_proc/60,2)));
msgbox('ROI detection finished');

end