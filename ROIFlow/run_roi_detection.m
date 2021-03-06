function run_roi_detection(infostruct, params)
%% Function for automatic detection of ROI defined by intensity
% Detect ROIs and save as struct to load and process further

%% Add paths
wd = pwd(); p = split(wd,'\'); p = p(1:end-1); pp = []; for ii=1:numel(p), pp = fullfile(pp,p{ii}); end
addpath(genpath(pp)), clear('p','pp');

nfiles = size(infostruct,2);
t_proc = 600;
rectype = questdlg('Select type of recording','Recording type','Confocal', '2PM', 'Confocal');
connectivity = 4;

%% Set parameters dependend on recording type if not given as input
% corr_threshold: range between 0.3 - 0.5 recommended if smoothed data
% are used, and 0.1 - 0.25 if unsmoothed
if strncmp(rectype,'C',1)
    if isempty(params)
        params = get_predefined_params('connc_px_thresh', 16,...
            'responding_px_factor', 1.15, 'px_exclusion', 0.25,...
            'lut_xlow', 1, 'lut_log_minval', 5, 'lut_log_maxval', 6.5, 'lut_log_growthrate', 0.2, 'lut_log_midpoint', 25,...
            'corr_threshold', 0.15, 'fill_thresh',0.25);
    end   
else
    if isempty(params)
        params = get_predefined_params('connc_px_thresh', 18,...
            'responding_px_factor', 1.15, 'px_exclusion', 0.25,...
            'lut_xlow', 1.4, 'lut_log_minval', 7, 'lut_log_maxval', 9, 'lut_log_growthrate', 0.3, 'lut_log_midpoint', 25,...
            'corr_threshold', 0.15, 'fill_thresh',0.25); % 7-10,0.2
    end
end

%% Create LUT
xval = 1:50;
cc_px_factor_lut = lut_create(xval, params);
xlut = params.lut_xlow+xval./100;
xluta = [1:0.1:xlut(1)-0.1];
xlutb = [xlut(end)+1:0.1:2.5];
cc_px_factor_lut = [repelem(cc_px_factor_lut(1),numel(xluta)) cc_px_factor_lut repelem(cc_px_factor_lut(end),numel(xlutb))];
cc_px_factor_lut = [xluta xlut xlutb; cc_px_factor_lut];
prctiles_rec_all = NaN(nfiles,1);

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
    bg_corr_val = infostruct(iF).corr_2pm;
    mocorr = infostruct(iF).mocorr;
    shift = params.dFoF_baseline_shift;
    xy_scale = infostruct(iF).xy_scale;
    winsize = round(params.dFoF_baseline_s/ft_s);
    if winsize < params.dFoF_baseline_frames, winsize = params.dFoF_baseline_frames; end
    if winsize > 20, winsize = 20; end
    if strncmp(rectype,'C',1)
        pmt_gain = NaN;
        offset = NaN;
    else
        pmt_gain = infostruct(iF).pmt_gain;
        offset = infostruct(iF).offset;
    end
    sd_factor = params.responding_px_factor;
    min_corr = params.corr_threshold;
    min_conn_px = params.connc_px_thresh;
    corr_smth_win = 3;
    
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
    if strncmp(rectype,'C',1), disp('Correct for background contamination');
    else, disp('Correct for dark count contamination'); 
    end
    if strncmp(rectype, '2',1)
        subtract_value = - bg_corr_val; % Negative sign necessary, if difference had been calculated as 0-offset - rec-specific offset
    else
        subtract_value = prctile(lintraces(exclude_px,:),5, 'all');
    end
    corr_traces = lintraces - subtract_value;
    
    
    %% Calculate area information
    fov_area = n_px*xy_scale*xy_scale;
    n_px_bg = sum(background);
    bg_area = n_px_bg*xy_scale*xy_scale;
    dend_area = fov_area-bg_area;
    
    %% Calculate dFoF stack
    disp('Calculate F/F values...');
    [~, FoFtraces, ~] = rollBase_dFoF(corr_traces, winsize,shift, 'roll');
    smth_traces = FoFtraces; %response_lin; % smooth_data(FoFtraces, ones(1,corr_smth_win));  %FoFtraces;
        
    %     tmpim = reshape(FoFtraces,[imsize(1) imsize(2) n_frames]);
    %     tmpp = strcat(filepath, '\FoF_', infostruct(iF).name, '.tif');
    %     for idx = 1:n_frames
    %         imwrite(uint8(tmpim(:,:,idx).*30), tmpp, 'tiff', 'WriteMode', 'append');
    %     end
    
    %% Detect responding px
    disp('Detect responding time points...');
    std_noise = NaN(n_px,1);
    for ipx = 1:n_px
        trc_neg = FoFtraces(ipx,:); trc_neg = trc_neg(trc_neg < 1);
        trc_synt_symm = [trc_neg 1+(1-trc_neg)];
        std_noise(ipx,1) = std(trc_synt_symm);
    end
    
    resp_thresh = 1+ std_noise .* sd_factor;
    response_lin = FoFtraces > resp_thresh;
    response_lin(exclude_px,:) = 0;
    response = reshape(response_lin,[imsize(1) imsize(2) n_frames]);
    
    %     test = zeros(size(response));
    %         tmpp = strcat(filepath, '\response_low_', infostruct(iF).name, '.tif');
    %         for idx = 1:n_frames
    %             imwrite(uint8(response(:,:,idx).*255), tmpp, 'tiff', 'WriteMode', 'append');
    %         end
    
    %% Detect in-frame components
    % Detect connected pixel (connectivity: 4) in each frame and save as linear index
    % list (2D array indices!)
    disp('Find in-frame components...');
    all_connc_cell_low=cell(n_frames,1);
    for t = 1:winsize, all_connc_cell_low{t,1} = []; end
    parfor t = winsize+1:n_frames
        bw=response(:,:,t);
        connc = bwconncomp(bw,connectivity);
        all_connc_cell_low{t,1} = connc.PixelIdxList;
    end
    
    %% Discard all components with less than minimum number of connected
    % pixel (preset threshold)
    disp('Filter components by size...');
    n_tot_cc_px = 0;
    connc_cell_low = cell(n_frames,1);
    for t = winsize+1:n_frames
        conn_px_low = cellfun(@numel, all_connc_cell_low{t,1}, 'UniformOutput', false);
        conn_px_low = cellfun(@(x) x, conn_px_low);
        connc_cell_low{t,1} = all_connc_cell_low{t,1}(conn_px_low > min_conn_px);
    end
    keep = ~cellfun(@isempty, connc_cell_low);
    connc_cell_low = connc_cell_low(keep,:);
    connc_frame_id = find(keep==1); % Indices of frames with connected components
    
    %% Derive adaptive factor based on SDnoise of component average traces
    % and discard components if its average response in frame t too small
    disp('Define recording-specific factor to filter in-frame components...');
    cc_n_tot = 0;
    for ii = 1:size(connc_cell_low,1), cc_n_tot=cc_n_tot+numel(connc_cell_low{ii,1}); end
    cc_trc = NaN(cc_n_tot,n_frames);
    cc_noise = NaN(cc_n_tot,1);
    cnt=0;
    for t = 1:size(connc_cell_low,1)
        for iCC = 1:numel(connc_cell_low{t,1}) % loop through components on this frame
            cnt = cnt + 1;
            tmp_px = connc_cell_low{t,1}{iCC};
            cc_trc(cnt,:) = mean(corr_traces(tmp_px,:),1);
        end
    end
    
    [~, cc_fof, ~] = rollBase_dFoF(cc_trc, winsize,shift, 'roll');
    for ii = 1:cc_n_tot
        cc_neg = cc_fof(ii,:); cc_neg = cc_neg(cc_neg < 1);
        cc_synt_symm = [cc_neg 1+(1-cc_neg)];
        cc_noise(ii,1) = std(cc_synt_symm);
    end
    
    cc_pooled = (cc_fof-1)./cc_noise;
    cc_pooled_pos = cc_pooled(cc_pooled>0);
    
    [~,midx] = min(abs(prctile(cc_pooled_pos, 80,'all') - cc_px_factor_lut(1,:)));
    cc_sd_factor_mltpl = cc_px_factor_lut(2,midx);
    cc_sd_factor = cc_sd_factor_mltpl * prctile(cc_pooled_pos, 80);
    prctiles_rec_all(iF) = prctile(cc_pooled_pos, 80);
    
    disp('Filter in-frame components based on response strength...');
    cnt = 0; del = [];
    for t = 1:size(connc_cell_low,1)
        keep = false(numel(connc_cell_low{t,1}),1);
        for iCC = 1:numel(connc_cell_low{t,1}) % loop through components on this frame
            cnt = cnt + 1;
            if cc_fof(cnt, connc_frame_id(t)) > (1+ cc_noise(cnt) * cc_sd_factor), keep(iCC) = true; end
        end
        if all(~keep), del = [del t];
        else
            connc_cell_low{t,1} = connc_cell_low{t,1}(keep);
            n_tot_cc_px = n_tot_cc_px+sum(cellfun(@numel, connc_cell_low{t,1}));
        end
    end
    connc_cell_low(del) = [];
    connc_frame_id(del) = []; % Indices of frames with connected components
    
    if cc_n_tot ~= 0
        %% Create matrix with information of all px within connected components
        px_cc_master_mtrx = NaN(n_tot_cc_px,4);
        % Master matrix with px information: col1: conn.component id, col2: row idx,
        % col3: column idx, col4: frame of conn.component
        pxcnter = 0; cc_n_tot=0;
        for t = 1:numel(connc_frame_id)
            ncc = numel(connc_cell_low{t,1}); % number connected components in this frame
            for cc = 1:ncc
                cc_n_tot = cc_n_tot+1;
                [r,c] = ind2sub([imsize(1) imsize(2)], connc_cell_low{t,1}{cc});
                n = numel(r); % Number of px in component
                px_cc_master_mtrx(pxcnter+1:pxcnter+n,:) = [repelem(cc_n_tot,n)' r c repelem(connc_frame_id(t),n)'];
                %             for ii = 1:size(r,1), test(r(ii),c(ii),connc_frame_id(t)) = 1; end
                pxcnter=pxcnter+n;
            end
        end
        
        %% Prepare for component combination
        cc_ids = unique(px_cc_master_mtrx(:,1));
        cc_n_curr = numel(cc_ids);
        expand8_yx = [-1 -1; -1 0; -1 1; 0 -1; 0 0; 0 1; 1 -1; 1 0; 1 1];
        expand4_yx = [ -1 0; 0 -1; 0 0; 0 1; 1 0];
        
        %% Detect frame-reduced components
        disp('Detect frame-reduced components...');
        flat_red_cc = zeros(imsize(1), imsize(2));
        for cc = 1:cc_n_curr
            cc_coords = px_cc_master_mtrx(px_cc_master_mtrx(:,1) == cc_ids(cc),2:3);
            for iPx = 1:size(cc_coords,1), flat_red_cc(cc_coords(iPx,1), cc_coords(iPx,2)) = 1; end
        end
        conn_comp = bwconncomp(flat_red_cc,connectivity);
        
        % Map in-frame to frame-red. components
        ccmap_px_mtrx = NaN(size(px_cc_master_mtrx,1),2); % Matrix to map in-frame to frame-red. component
        ccmap_px_mtrx(:,1) = px_cc_master_mtrx(:,1); % First column: in-frame components
        for cc1 = 1:cc_n_curr
            tmppter = find(px_cc_master_mtrx(:,1) == cc_ids(cc1));
            tmplinidx = sub2ind([imsize(1) imsize(2)], px_cc_master_mtrx(tmppter(1),2), px_cc_master_mtrx(tmppter(1),3));
            for cc2 = 1:size(conn_comp.PixelIdxList,2)
                if any(conn_comp.PixelIdxList{cc2} == tmplinidx), ccmap_px_mtrx(tmppter,2) = cc2; end
            end
        end
        
        %% Check correlation of in-frame components
        disp('Correlation between in-frame components within frame-reduced component...');
        for cc1 = 1:size(conn_comp.PixelIdxList,2) % Loop through frame-red. components
            inframe_id = unique(ccmap_px_mtrx(ccmap_px_mtrx(:,2) == cc1,1));
            tmp_n_cc = numel(inframe_id);
            adj_mtrx = false(tmp_n_cc, tmp_n_cc);
            corr_mtrx = zeros(tmp_n_cc, tmp_n_cc);
            for i = 1:tmp_n_cc, corr_mtrx(i,i) = 1; end
            
            % Loop through in-frame components of current frame-red.
            % component
            for r = 1:tmp_n_cc-1
                tmpid1 = inframe_id(r);
                tmpcc1 = px_cc_master_mtrx(px_cc_master_mtrx(:,1) == tmpid1,2:3); % y x
                
                % Find XY-adjacent/overlapping in-frame components
                % Expand each px of current component
                exp_px_1 = [];
                for ipx = 1:size(tmpcc1,1), exp_px_1 = [exp_px_1; tmpcc1(ipx,:) + expand4_yx]; end
                
                % Delete multiple occurrences in exp_px
                for ipx = 1:size(exp_px_1,1)
                    del = find(all((exp_px_1 - exp_px_1(ipx,:))==0,2) == 1);
                    if numel(del) > 1, exp_px_1(del(2:end),:) = NaN; end
                end
                exp_px_1(isnan(exp_px_1(:,1)),:) = [];
                
                for c = r+1:tmp_n_cc
                    tmpid2 = inframe_id(c);
                    tmpcc2 = px_cc_master_mtrx(px_cc_master_mtrx(:,1) == tmpid2,2:3);
                    
                    % Check if any other component shares px with the expanded
                    % region
                    adj_pter = false(size(tmpcc2,1),1);
                    for ipx = 1:size(exp_px_1,1)
                        adj_pter = adj_pter | all((tmpcc2 - exp_px_1(ipx,:))==0,2); % Mark every px of other components with 1 if adjacent/overlapping
                    end
                    % Mark in adjacent matrix
                    if any(adj_pter), adj_mtrx(r,c) = true; adj_mtrx(c, r) = true; end
                    
                    % Check if overlap exists
                    ovlap_cc2 = false(size(tmpcc2,1),1);
                    for ipx = 1:size(tmpcc2,1)
                        if any(all((tmpcc1 - tmpcc2(ipx,:))==0,2)), ovlap_cc2(ipx) = true; end
                    end
                    cnt_overlap = sum(ovlap_cc2);
                    
                    test_cc1 = tmpcc1; test_cc2 = tmpcc2;
                    if cnt_overlap > 0
                        % Expand overlapping region
                        ovlap_info = tmpcc2(ovlap_cc2,:);
                        exp_px = [];
                        for ipx = 1:size(ovlap_info,1), exp_px = [exp_px; ovlap_info(ipx,:) + expand8_yx]; end
                        % Delete multiple occurrences in exp_px
                        for ipx = 1:size(exp_px,1)
                            del = find(all((exp_px - exp_px(ipx,:))==0,2) == 1);
                            if numel(del) > 1, exp_px(del(2:end),:) = NaN; end
                        end
                        exp_px(isnan(exp_px(:,1)),:) = [];
                        % Delete overlapping px + 1-px buffer zone from
                        % components for correlation test
                        for ipx  = 1:size(exp_px,1)
                            del = find(all((test_cc1 - exp_px(ipx,:))==0,2) == 1);
                            test_cc1(del,:) = [];
                            del = find(all((test_cc2 - exp_px(ipx,:))==0,2) == 1);
                            test_cc2(del,:) = [];
                        end
                    end
                    
                    % Test correlation of the average trace of pairs of
                    % in-frame components within the current frame-red. component
                    if ~isempty(test_cc1) && ~isempty(test_cc2)
                        lindices1 = sub2ind([imsize(1) imsize(2)], test_cc1(:,1), test_cc1(:,2));
                        avtrc1 = mean(smth_traces(lindices1,:),1);
                        lindices2 = sub2ind([imsize(1) imsize(2)], test_cc2(:,1), test_cc2(:,2));
                        avtrc2 = mean(smth_traces(lindices2,:),1);
                        av1 = mean(avtrc1); av2 = mean(avtrc2);
                        corr_mtrx(r,c) = sum((av1-avtrc1).*(av2-avtrc2)) / sqrt(sum((av1-avtrc1).^2)*sum((av2-avtrc2).^2));
                        corr_mtrx(c,r) = corr_mtrx(r,c);
%                         figure, plot(avtrc1), hold on, plot(avtrc2);
                    else
                        corr_mtrx(r,c) = 1;
                        corr_mtrx(c,r) = 1;
                    end
                end
            end
            %             corr_mtrx = corr_mtrx + rot90(fliplr(corr_mtrx)); % Create fully diagonal matrix
            %             if numel(corr_mtrx) >= 16, figure(); histogram(corr_mtrx, 10); end
            
            %% Check if in-frame component has sufficiently high correlation
            % with all other in-frame components to integrate it into current frame-red.
            % component
            check_pter = [];
            if ~all(corr_mtrx > min_corr, 'all') %all(corr_mtrx(r,:) <= params.corr_threshold | corr_mtrx(r,:) == 1)
                check_pter = 1:tmp_n_cc;
            end
            
            %% Check if in-frame components which were not included in current frame-red. component
            % are adjacent and split from current frame-red. component
            disp('Combine in-frame components if sufficiently correlated...');
            if ~isempty(check_pter)
                n_check = numel(check_pter);
                check_adj_mtrx = false(n_check, n_check);
                check_corr_mtrx = zeros(n_check, n_check);
                %                 for i = 1:n_check, check_adj_mtrx(i,i) = 1; end
                %                 for i = 1:n_check, check_corr_mtrx(i,i) = 1; end
                
                % Adjacent in-frame comp.
                for r = 1:n_check-1
                    for c = r+1:n_check
                        % Correlation of in-frame comp.
                        check_corr_mtrx(r,c) = corr_mtrx(check_pter(r),check_pter(c));
                        check_corr_mtrx(c,r) = corr_mtrx(check_pter(r),check_pter(c));
                        % Check if adjacent/overlapping
                        if adj_mtrx(check_pter(r),check_pter(c))
                            check_adj_mtrx(r,c) = true;
                            check_adj_mtrx(c,r) = true;
                        end
                    end
                end
                
                % Combine based on the highest correlation
                % Operate only on components 1-n_check incl in check_pter
                % and reference them by 1-n_check
                check_adjcorr_mtrx = check_corr_mtrx;
                check_adjcorr_mtrx(check_adjcorr_mtrx <= min_corr) = 0; % set all correlations <= threshold to 0
                check_adjcorr_mtrx(~check_adj_mtrx) = 0; % set all correlations of non-adjacent pairs to 0
                res_cc_bin = true(1, n_check);
                add_cc_cnter = max(ccmap_px_mtrx(:,2));
                
                while any(res_cc_bin)
                    % First combined component consists of the adjacent
                    % components with the highest correlation
                    [maxcorr, maxidx] = max(check_adjcorr_mtrx,[],'all', 'linear');
                    if maxcorr ~= 0 % there exist adjacent components with sufficient correlation
                        [cc1, cc2] = ind2sub([n_check n_check], maxidx);
                        conn_cc = [cc1, cc2];
                        cc_endpoints = conn_cc;
                        res_cc_bin(conn_cc) = false;
                        while ~isempty(cc_endpoints)
                            % Find components adjacent to current component
                            bridge_cc = cc_endpoints(1);
                            pot_partners = find(check_adjcorr_mtrx(bridge_cc,:) > 0 & res_cc_bin);
                            [~,sortidx] = sort(check_adjcorr_mtrx(bridge_cc,check_adjcorr_mtrx(bridge_cc,:) > 0 & res_cc_bin),'descend');
                            pot_partners = pot_partners(sortidx);
                            %                             for ii = 1:cc_endpoints, pot_partners(pot_partners == cc_endpoints(ii)) = []; end
                            if ~isempty(pot_partners)
                                for ii = 1:numel(pot_partners)
                                    % If component to be tested has no adjacent component
                                    % to which it is higher correlated and if its
                                    % correlation with all other connected components is
                                    % sufficiently high, add it to complex
                                    if all(check_corr_mtrx(pot_partners(ii),conn_cc) > min_corr)
                                        conn_cc = [conn_cc pot_partners(ii)]; % add to combined components
                                        cc_endpoints = [cc_endpoints pot_partners(ii)]; % add as new endpoint
                                        res_cc_bin(pot_partners(ii)) = false; % remove from residual components
                                    end
                                end
                            end
                            cc_endpoints(cc_endpoints == bridge_cc) = [];
                        end
                        
                        % Mark as belonging to new frame-red. component
                        add_cc_cnter = add_cc_cnter + 1;
                        for ii = 1:numel(conn_cc)
%                             fprintf('%i: %i \n', add_cc_cnter, conn_cc(ii));
                            ccmap_px_mtrx(ccmap_px_mtrx(:,1) == inframe_id(check_pter(conn_cc(ii))),2) = add_cc_cnter;
                            % Remove from adjacency correlation matrix
                            check_adjcorr_mtrx(conn_cc(ii),:) = 0;
                            check_adjcorr_mtrx(:,conn_cc(ii)) = 0;
                        end
                    else
                        single_cc = find(res_cc_bin, 1, 'first');
                        add_cc_cnter = add_cc_cnter + 1;
%                         fprintf('%i: %i \n', add_cc_cnter, single_cc);
                        ccmap_px_mtrx(ccmap_px_mtrx(:,1) == inframe_id(check_pter(single_cc)),2) = add_cc_cnter;
                        res_cc_bin(single_cc) = false;
                        % Remove from adjacency correlation matrix
                        check_adjcorr_mtrx(single_cc,:) = 0;
                        check_adjcorr_mtrx(:,single_cc) = 0;
                    end
                end
            end
        end
        px_cc_master_mtrx(:,1) = ccmap_px_mtrx(:,2);
        px_cc_master_mtrx(:,end+1) = ccmap_px_mtrx(:,1);
        cc_n_curr = numel(unique(px_cc_master_mtrx(:,1)));
        
        %% Fill gaps
        disp('Fill holes in components...');
        if ~isempty(px_cc_master_mtrx)
            px_cc_master_mtrx = fill_roi_gaps(imsize, px_cc_master_mtrx, [2 3],params.fill_thresh);
        end
        
        %% Delete left overlapping parts
        disp('Delete residual overlap...');
        ccs = unique(px_cc_master_mtrx(:,1));
        for cc1 = 1:numel(ccs)-1
            tmpidx1 = px_cc_master_mtrx(:,1) == ccs(cc1);
            tmpcc1 = px_cc_master_mtrx(tmpidx1,2:3);
            if ~isempty(tmpcc1)
                for cc2 = cc1+1:numel(ccs)
                    tmpidx2 = px_cc_master_mtrx(:,1) == ccs(cc2);
                    tmpcc2 = px_cc_master_mtrx(tmpidx2,2:3);
                    
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
%                             bw1 = zeros(imsize(1),imsize(2),1);
%                             for ipx = 1:size(tmpcc1,1), bw1(tmpcc1(ipx,1), tmpcc1(ipx,2)) = 1; end
%                             bw2 = zeros(imsize(1),imsize(2),1);
%                             for ipx = 1:size(tmpcc2,1), bw2(tmpcc2(ipx,1), tmpcc2(ipx,2)) = 1; end
%                             figure, imshow(bw1);
%                             figure, imshow(bw2);
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
                                avtrc1 = mean(smth_traces(lindices1,:),1);
                                lindices2 = sub2ind([imsize(1) imsize(2)], test_cc2(:,1), test_cc2(:,2));
                                avtrc2 = mean(smth_traces(lindices2,:),1);
                                lindices_ovlap = sub2ind([imsize(1) imsize(2)], ovlap_info(:,1), ovlap_info(:,2));
                                avtrc_ovlap = mean(smth_traces(lindices_ovlap,:),1);
                                av1 = mean(avtrc1); av2 = mean(avtrc2); av_ovlap = mean(avtrc_ovlap);
                                r1 = sum((av1-avtrc1).*(av_ovlap-avtrc_ovlap)) / sqrt(sum((av1-avtrc1).^2)*sum((av_ovlap-avtrc_ovlap).^2));
                                r2 = sum((av2-avtrc2).*(av_ovlap-avtrc_ovlap)) / sqrt(sum((av2-avtrc2).^2)*sum((av_ovlap-avtrc_ovlap).^2));                           
                                if r1 > r2
                                    cc2_idxlist = find(tmpidx2==1);
                                    px_cc_master_mtrx(cc2_idxlist(ovlap_cc2),:) = NaN;
                                else
                                    cc1_idxlist = find(tmpidx1==1);
                                    px_cc_master_mtrx(cc1_idxlist(ovlap_cc1),:) = NaN;
                                end
                            elseif isempty(test_cc1) && isempty(test_cc2)
                                px_cc_master_mtrx(tmpidx1,:) = NaN;
                                px_cc_master_mtrx(tmpidx2,:) = NaN;
                            elseif isempty(test_cc1)
                                px_cc_master_mtrx(tmpidx1,:) = NaN;
                            elseif isempty(test_cc2)
                                px_cc_master_mtrx(tmpidx2,:) = NaN;
                            end
                        end
                    end
                end
            end
            px_cc_master_mtrx(isnan(px_cc_master_mtrx(:,1)),:) = [];
        end
        
        %% Delete all components that are too small
        disp('Filter components by size...');
        ccs = unique(px_cc_master_mtrx(:,1));
        del = false(size(px_cc_master_mtrx,1),1);
        cc = 1;
        while cc <= numel(ccs)
            cc_px = px_cc_master_mtrx(px_cc_master_mtrx(:,1) == ccs(cc),2:3);
            bw = zeros(imsize(1),imsize(2),1);
            for ipx = 1:size(cc_px,1), bw(cc_px(ipx,1), cc_px(ipx,2)) = 1; end
            connc = bwconncomp(bw,connectivity);
            nsub = size(connc.PixelIdxList,2);
            if nsub > 1
                [pxnum, idx] = max(cellfun(@numel, connc.PixelIdxList));
                tmpdel = true(nsub,1);
                if pxnum > min_conn_px, tmpdel(idx) = false; end
                for ccc = 1:nsub
                    if tmpdel(ccc)
                        [delpx_r,delpx_c] = ind2sub([imsize(1) imsize(2)], connc.PixelIdxList{ccc});
                        for ipx = 1:size(delpx_r,1)
                            del(all(px_cc_master_mtrx(:,2:3) - [delpx_r(ipx) delpx_c(ipx)] == 0,2)) = true;
                        end
                    end
                end
            elseif size(cc_px,1) <= min_conn_px
                del(px_cc_master_mtrx(:,1) == ccs(cc),:) = true;
            end
            cc = cc+1;
        end
        px_cc_master_mtrx(del,:) = [];
        
        if ~isempty(px_cc_master_mtrx)
            roi_ids = unique(px_cc_master_mtrx(:,1));
            %% Renumber components
            new_roi_ids = NaN(size(px_cc_master_mtrx,1),1);
            for iNum = 1:numel(roi_ids)
                new_roi_ids(px_cc_master_mtrx(:,1) == roi_ids(iNum)) = iNum;
            end
            px_cc_master_mtrx(:,1) = new_roi_ids;
            cc_n_curr = numel(unique(new_roi_ids));
            
            %% Extract component area and centroids
            disp('Extract ROI area and centroids...');
            idx_mtrx = px_cc_master_mtrx;
            roi_area_um = NaN(cc_n_curr,1);
            roi_xy_center_um = NaN(cc_n_curr,2);
            for cc = 1:cc_n_curr
                cc_px = idx_mtrx(idx_mtrx(:,1) == cc,2:3);
                roi_area_um(cc,1) = size(cc_px,1)*xy_scale*xy_scale;
                roi_xy_center_um(cc,1) = mean(cc_px(:,2))*xy_scale;
                roi_xy_center_um(cc,2) = mean(cc_px(:,1))*xy_scale;
            end
            
            %% Sort ROIs by distance to origin
            dist_info = sqrt(roi_xy_center_um(:,1).^2 + roi_xy_center_um(:,2).^2);
            [~, sorted_id] = sort(dist_info, 'ascend');
            % Apply sorting
            roi_area_um = roi_area_um(sorted_id,:);
            roi_xy_center_um = roi_xy_center_um(sorted_id,:);
            
            px_cc_sorted_mtrx = px_cc_master_mtrx;
            new_id = NaN(size(px_cc_sorted_mtrx,1),1);
            for cc = 1:cc_n_curr
                new_id(px_cc_sorted_mtrx(:,1) == sorted_id(cc),1) = cc;
            end
            px_cc_sorted_mtrx(:,1) = new_id; clear('new_id');
            [~, sorted] = sort(px_cc_sorted_mtrx(:,1),'ascend');
            px_cc_sorted_mtrx = px_cc_sorted_mtrx(sorted,:);
%             px_cc_master_mtrx = px_cc_sorted_mtrx;
            
            %% Extract ROI boundaries - combined ROIs
            disp('Find ROI boundaries...');
            bounds_cell_comb = cell(cc_n_curr,2);
            cnt=1; cc=1;
            idx_mtrx = px_cc_sorted_mtrx;
            while cc <= cc_n_curr
                cc_px = idx_mtrx(idx_mtrx(:,1) == cc,2:3);
                bw = zeros(imsize(1),imsize(2),1);
                for ipx = 1:size(cc_px,1), bw(cc_px(ipx,1), cc_px(ipx,2)) = 1; end
                [B,~] = bwboundaries(bw,8,'noholes');  % Find boundary
                if size(B,1) > 1
                    disp('A');
                end
                bounds_cell_comb{cnt,1} = B;
                % Color as indicator of number of individual components
                bounds_cell_comb{cnt,2} = numel(unique(idx_mtrx(idx_mtrx(:,1) == cc,end)));
                %             bounds_cell_comb{cnt,2} = cc;
                cnt = cnt+1;
                cc = cc+1;
            end
            bounds_cell_comb = bounds_cell_comb(1:cnt-1,:);
            
            %% Extract traces
            % dF/F calculated using rolling average for F_b and F_0
            % F/F calculated as corrected raw intensity normalized by average
            % across timepoints
            disp('Extract average ROI traces...');
            roi_av_traces = NaN(cc_n_curr,n_frames);
            roi_av_FoF_traces = NaN(cc_n_curr,n_frames);
            %         roi_av_dFoF_traces = NaN(cc_n_curr,n_frames);
            mean_traces = NaN(cc_n_curr,n_frames);
            idx_mtrx = px_cc_sorted_mtrx;
            for cc = 1:cc_n_curr
                cc_px = idx_mtrx(idx_mtrx(:,1) == cc,2:3);
                px_lin = sub2ind([imsize(1) imsize(2)],cc_px(:,1),cc_px(:,2));
                roi_av_traces(cc,:) = mean(corr_traces(px_lin,:),1);
                %             roi_av_dFoF_traces(cc,:) = mean(dFoFtraces(px_lin,:),1);
                av = mean(roi_av_traces(cc,:));
                roi_av_FoF_traces(cc,:) = roi_av_traces(cc,:)./ av;
                mean_traces(cc,:) = mean(corr_traces(px_lin,:),1);
            end
            [~, ~,roi_av_dFoF_traces] = rollBase_dFoF(mean_traces, winsize,shift, 'roll');
            
            %% AIP & MIP Projection
            corr_aip = reshape(mean(corr_traces,2),[imsize(1) imsize(2)]);
            tmpmin = min(corr_aip,[],'all');
            uint16_aip = corr_aip-tmpmin; tmpmax = max(uint16_aip,[],'all');
            uint16_aip = (uint16_aip./tmpmax).*(2^8);
            corr_mip = reshape(max(corr_traces,[],2),[imsize(1) imsize(2)]);
            tmpmin = min(corr_mip,[],'all');
            uint16_mip = corr_mip-tmpmin; tmpmax = max(uint16_mip,[],'all');
            uint8_mip = (uint16_mip./tmpmax).*(2^8);
            clear('tmpmin','tmpmax');
            
            %         % Save overview of all ROIs using AIP
            scrsz = get(0,'ScreenSize');
            red=[1 0 1]; green=[0 1 1]; blue=[0 .7 1];
            cmap_n = max(cellfun(@max, bounds_cell_comb(:,2)));
            cmap = get_colormap(red,green,blue,cmap_n);
            
            % Save overview of all ROIs using MIP
            ovwfig_mip = figure('Position', [0 0 scrsz(3:4)], 'PaperSize', scrsz(3:4));
            imagesc(uint8_mip);colormap('gray');axis('off');daspect([1 1 1]);
            for cc = 1:cc_n_curr
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
            
            %% Pack information and save
            disp('Save information...');
            roidata.filepath = filepointer;
            roidata.frametime_s = ft_s;
            roidata.xy_scale = xy_scale;
            roidata.baseline_frames = winsize;
            roidata.bg_threshold = bg_thresh;
            roidata.subtract_value = subtract_value;
            roidata.offset = offset;
            roidata.pmt_gain = pmt_gain;
            roidata.params = params;
            roidata.mocorr_run = mocorr;
            roidata.aip = corr_aip;
            roidata.mip = corr_mip;
            roidata.background_mask = background;
            roidata.response_multiplier = [sd_factor cc_sd_factor];
            roidata.dendritic_area_um = dend_area;
            roidata.background_area_um = bg_area;
            roidata.fov_area_um = fov_area;
            roidata.roi_bounds = bounds_cell_comb;
            roidata.roi_idx = unique(px_cc_master_mtrx(:,1));
            roidata.n_rois = cc_n_curr;
            roidata.dFoF_traces = roi_av_dFoF_traces;
            roidata.FoF_traces = roi_av_FoF_traces;
            roidata.traces = roi_av_traces;
            roidata.roi_centroids_um = roi_xy_center_um;
            roidata.roi_area_um = roi_area_um;
            
            outputname = strcat(savepath,'\ROIdata_',filename,'.mat');
            if isa(outputname,'cell'), outputname= outputname{1}; end
            save(outputname,'roidata', '-v7.3');
            clear('roidata', 'corr_mip','corr_aip','filepointer','mocorr','bg_thresh',...
                'ft_s','resp_thresh','background','bounds_cell_comb');
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
        roidata.pmt_gain = pmt_gain;
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
    
    %         tmpp = strcat(filepath, '\bin_components_', infostruct(iF).name, '.tif');
    %         for idx = 1:n_frames
    %             imwrite(uint8(test(:,:,idx).*255), tmpp, 'tiff', 'WriteMode', 'append');
    %         end
    
    clear('all_connc_cell_low', 'all_connc_cell_high', 'background', 'bg_area', 'bg_thresh',...
        'bleachcorr', 'ccs', 'ccter', 'cnter',...
        'conn_px_low', 'conn_px_high', 'connc_cell_low', 'connc_cell_high', 'connc_frame_id',...
        'corr_aip', 'corr_mip', 'corr_traces','rolldFoFtraces',...
        'dend_area', 'dFoFtraces', 'fov_area', 'ft_s', 'ii', 'imsize', 'inputim',...
        'keep', 'lintraces', 'mpcprr', 'n_frames', 'n_px', 'n_px_bg', 'n_tot_cc_px', 'resp_thresh',...
        'response', 'response_high', 'px_cc_master_mtrx', 'roidata', 'subtract_value', 't',...
        'tot_num_cc', 'xy_scale');
end
t_proc = toc(t_all);
fprintf('Total processing time: %s min\n', num2str(round(t_proc/60,2)));

figure(); title('Responding component filter: LUT & SDnoise 80%-Percentile values');
yyaxis left
histogram(prctiles_rec_all,cc_px_factor_lut(1,:),'FaceAlpha',.3, 'EdgeAlpha',.3);
ylabel('Counts, 80%-Percentile values');
yyaxis right
plot(cc_px_factor_lut(1,:),cc_px_factor_lut(2,:));
ylim([1 11]);
ylabel('LUT - SD multiple - component threshold');
outputname = strcat(savepath,'\Resp_component_percentiles_LUT.png');
if isa(outputname,'cell'), outputname= outputname{1}; end
saveas(gcf, outputname);
pause(0.3);
close gcf

msgbox('ROI detection finished');

end