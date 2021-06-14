function mtrx_histo = Threshold_test_run_roi_detection(infostruct, params)
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
        params = get_predefined_params('connc_px_thresh', 15,'responding_px_factor', 1.5,...
            'px_exclusion', 15, 'corr_threshold', 0.15, 'fill_thresh',0.25);
    end
else
    if isempty(params)
        params = get_predefined_params('connc_px_thresh', 15,'responding_px_factor', 1.5,...
            'px_exclusion', 15, 'corr_threshold', 0.15, 'fill_thresh',0.25);
    end
end

set_fac = [1.3 1.5 1.75 1.85];

sets_n = size(set_fac,2);

cell_pxcnt = cell(nfiles*sets_n,1);
cellcnter = 1;

edges = 0:1:50;

%% Main loop
t_all = tic;
for iF = 1:nfiles
    roidata = struct();
    
    fprintf('\n***************\nProcessing file %i of %i | Expected time remaining: %s min\n***************\n', iF,nfiles,num2str(round(nfiles*t_proc/60,2)));
    
    %% Unpack information
    filepointer = infostruct(iF).path;
    [filepath,filename, ~] = fileparts(filepointer);
    savepath = fullfile(filepath,'ROIFlow_Analysis');
    if ~isfolder(savepath), mkdir(savepath); end
    ft_s = infostruct(iF).frametime_s;
    bg_thresh = infostruct(iF).bg_threshold;
    bg_corr_val = infostruct(iF).bg_2pm_corr;
    mocorr = infostruct(iF).mocorr;
    shift = params.dFoF_baseline_shift;
%     sd_factor = params.responding_px_factor;
    %     bleachcorr = infostruct(iF).bleachcorr;
    xy_scale = infostruct(iF).xy_scale;
    winsize = round(params.dFoF_baseline_s/ft_s);
    if winsize < params.dFoF_baseline_frames, winsize = params.dFoF_baseline_frames; end
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
    
    %% Classify as background
    aip = mean(lintraces,2);
    background = true(n_px,1);
    background(aip > bg_thresh) = false;
    clear('aip');
    
    %% Calculate area information
    fov_area = n_px*xy_scale*xy_scale;
    n_px_bg = sum(background);
    bg_area = n_px_bg*xy_scale*xy_scale;
    dend_area = fov_area-bg_area;
    
    %% Subtract dark count or background
    if strncmp(rectype, 'C',1)
        bgvalues = lintraces(background==1);
        [c_abs,e] = histcounts(bgvalues); % Count frequency of intensity values
        c_rel = c_abs./sum(c_abs); % Relative frequency
        c_cum = cumsum(c_rel); % Cummulative frequency
        [~,idx] = min(abs(c_cum - params.subtract_perc));
        subtract_value = e(idx);
    else
        subtract_value = - bg_corr_val; % Negative sign necessary, if difference had been calculated as 0-offset - rec-specific offset
    end
    
    corr_traces = lintraces - subtract_value;
    
    %% Calculate dFoF stack
    disp('Calculate dF/F stack...');
    [~, FoFtraces, dFoFtraces] = rollBase_dFoF(corr_traces, winsize,shift, 'px');
    
    %% Detect responding px
    disp('Detect responding pixel...');
    % Criterion 1
    sd_alltrc = std(corr_traces,[],2);
    sd_thresh = prctile(sd_alltrc, params.px_exclusion);
    exclude_px = sd_alltrc < sd_thresh;
   
    for iSet = 1:sets_n
        fprintf('Set: %i\n',iSet);
        sd_factor = set_fac(iSet);
        % Criterion 2
        std_noise = NaN(n_px,1);
        numHiVals = NaN(n_px,1); numLoVals = NaN(n_px,1);
        for ipx = 1:n_px
            trc_neg = FoFtraces(ipx,:); trc_neg = trc_neg(trc_neg < 1);
            trc_synt_symm = [trc_neg 1+(1-trc_neg)];
            std_noise(ipx,1) = std(trc_synt_symm);
            numHiVals(ipx,1) = sum(FoFtraces(ipx,:) > 1+std_noise(ipx,1)*sd_factor);
            numLoVals(ipx,1) = sum(FoFtraces(ipx,:) < 1-std_noise(ipx,1)*sd_factor);
        end
        
        resp_thresh = 1+ std_noise .* sd_factor;
        response = FoFtraces > resp_thresh;
        response(exclude_px,:) = 0;
        response = reshape(response,[imsize(1) imsize(2) n_frames]);
        
        %% ******** Assign px to ROI ********
        %% Component detection in 2D
        disp('Find connected components...');
        % Detect connected pixel (connectivity: 8) in each frame and save as linear index
        % list (2D array indices!)
        all_connc_cell=cell(n_frames,1);
        tic;
        for t = 1:winsize
            all_connc_cell{t,1} = [];
        end
        parfor t = winsize+1:n_frames
            bw=response(:,:,t);
            connc = bwconncomp(bw,4);
            all_connc_cell{t,1} = connc.PixelIdxList;
        end
        toc
        
        % Discard all components with less than minimum number of connected
        % pixel (preset threshold)
        tic;
        num_cc_px = 0;
        connc_cell = cell(n_frames,1);
        all_n = zeros(1,numel(edges)-1);
        for t = winsize+1:n_frames
            conn_px = cellfun(@numel, all_connc_cell{t,1}, 'UniformOutput', false);
            conn_px = cellfun(@(x) x, conn_px);
            all_n = all_n + histcounts(conn_px,edges);
        end
        toc
        
        cell_pxcnt{cellcnter} = all_n;
        cellcnter = cellcnter+1;
    end
    
end

cellcnter = cellcnter-1;
mtrx_histo = NaN(cellcnter+1,numel(edges)-1);
mtrx_histo(1,:) = edges(2:end);
for ii = 1:cellcnter
    mtrx_histo(ii+1,:) = cell_pxcnt{ii};
end



end