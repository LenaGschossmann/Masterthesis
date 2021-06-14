%% ** Script to calculate relative position in molecular layer of ROIs in PC dendritic tree **

clear all;
addpath('C:\Users\lena_\Projects\code_extern\Matlab\WriteImageJRoi');

exp_dir_path = uigetdir('C:\Users\lena_\Projects\iGluSnFR\Analysis\Confocal\ROIFlow');
info_dir_path = uigetdir('C:\Users\lena_\Projects\iGluSnFR\Analysis\Confocal\ML_Position');

%% 1) Select folder of analyzed event data of recordings within one experiment
exp_dir = dir(exp_dir_path);
isnotdir = [exp_dir.isdir] == 0;
exp_dir(isnotdir) = []; exp_dir(1:2) = [];
nrecs = size(exp_dir,1);
rec_names = cell(nrecs,1); rec_paths = cell(nrecs,1);
for iRec = 1:nrecs
    rec_names{iRec} = exp_dir(iRec).name;
    rec_paths{iRec} = fullfile(exp_dir(1).folder, rec_names{iRec});
end 

% Load information
rec_cell = cell(nrecs,13);
% col1: roi names, col2: ctr coordinates, col3-6: FOV info (fov ctr x,y,width,px res), col7-8: reference,
% col9: ctr relative to ref, col10: FOV corner points, col11: rel.ML pos,
% ML thickness, col12: rel. ML pos corner points, col13: ROI parameters
del_rec = [];
for iRec = 1:nrecs
    tmpinfo = dir(strcat(rec_paths{iRec}, '\SUMMARY*'));
    try
        tmp = readtable(fullfile(tmpinfo.folder, tmpinfo.name));
        rec_cell{iRec,1} = table2cell(tmp(:,1));
        rec_cell{iRec,2} = table2array(tmp(:,16:17));
        rec_cell{iRec,13} = table2array(tmp(:,[2,5,12,14,15,18]));
    catch
        del_rec = [del_rec iRec];
    end
end
rec_names(del_rec) = [];
rec_paths(del_rec) = [];
rec_cell(del_rec,:) = [];
nrecs = size(rec_names,1);

%% 2) Select .xls with FOV position of selected recordings & coordinates of PC somata in reference
% Table for FOV position (of timeseries and references:
% .xls with column 1 containing recording ID,
% column 2 containing px width of FOV, column 3 containing px resolution in
% µm, column 4&5 containing the x and y coordinate of FOV center (- is left
% (x) of ctr and above (y) reference ctr)
% ! Names must be same as filenames

% Table for PC somata coordinates: .xls with column 1&2 containing x&y
% coordinate, column 3 containing angle (irrelevant), column 4 containing
% length; last 3 rows line ROIs marking molecular layer thickness

% ! Make sure in folder are only those 1 table for FOV position and tables
% with somata coordinates!

info_dir = dir(info_dir_path);
info_dir(1:2) = [];
del=[];
for i = 1:size(info_dir,1), if strncmp(info_dir(i).name, 'ML_pos_Ma',9), del = [del i]; end, end
info_dir(del) = []; 
nfiles = size(info_dir,1);
info_paths = cell(nfiles,1);
fov_cell = []; % contains information of FOV position, resolution, and reference pointer
ref_cell = cell(nfiles-1,8); % contains coordinates of somata and length of ML
refcnt = 1;
for iF = 1:nfiles
    if strcmp(info_dir(iF).name, 'fov_pos.xlsx')
        fov_cell = readtable(fullfile(info_dir(1).folder, info_dir(iF).name));
        fovrows = size(fov_cell,1);
    else
        ref_cell{refcnt,1} = info_dir(iF).name;
        ref_cell{refcnt,2} = readmatrix(fullfile(info_dir(1).folder, info_dir(iF).name));
        refcnt = refcnt+1;
    end
end
refcnt = refcnt-1;

%% Assign FOV info
% Split fov_cell
ref_fov = strcmp(fov_cell{:,6},'');
fov_cell_ref = table2cell(fov_cell(ref_fov,1:3));
fov_cell_rec = table2cell(fov_cell(~ref_fov,:));
del_rec = [];
for i = 1:size(fov_cell_rec,1)
    rec=[]; ref=[];
    check = regexp(rec_names, fov_cell_rec{i,1});
    for ii = 1:size(check,1), if ~isempty(check{ii}), rec = ii; end, end
    if ~isempty(rec)
        rec_cell(rec,3:7) = fov_cell_rec(i,2:6); % Add information about FOV coordinates to respective rec.
        % Add index of reference file
        check = regexp(ref_cell(:,1), fov_cell_rec{i,end});
        for ii = 1:size(check,1), if ~isempty(check{ii}), ref = ii; end, end
        if ~isempty(ref)
            rec_cell{rec,8} = ref; % Add index of reference belonging to respective rec. (in ref_cell)
        else
            fprintf('NO reference found for: %s.\n', rec_names{rec});
            del_rec = [del_rec rec];
        end
    else
        fprintf('NO fitting recording detected for: %s (row %i).\n', fov_cell_rec{i,end}, i);
    end
end
rec_names(del_rec) = [];
rec_paths(del_rec) = [];
rec_cell(del_rec,:) = [];
nrecs = size(rec_names,1);

for i = 1:size(fov_cell_ref,1)
    sel = [];
    check = regexp(ref_cell(:,1), fov_cell_ref{i,1});
    for ii = 1:size(check,1), if ~isempty(check{ii}), sel = ii; end, end
    if ~isempty(sel)
        % Split information about soma centers and ML length
        tmp = ref_cell{sel,2}; % Unpack information
        ref_cell{sel,3} = tmp(tmp(:,4) == 0,2:3); % Soma ctr coordinates
        ref_cell{sel,4} = tmp(tmp(:,4) ~= 0,5); % ML length
        ref_cell(sel,5:6) = fov_cell_ref(i,2:3); % Px width and resolution
    else
        fprintf('Not a reference: %s.\n', fov_cell_ref{i,1});
    end
end

%% 3) Load boundaries of ROIs
roidata_dir = dir(strcat(exp_dir_path,'\ROIDATA*.mat'));
roi_cell = cell(nrecs,3);
for i = 1:size(roidata_dir,1)
    idx = [];
    check = regexp(rec_names, roidata_dir(i).name(9:end-4));
    for ii = 1:size(check,1), if ~isempty(check{ii}), idx = ii; end, end
    if ~isempty(idx)
        load(fullfile(roidata_dir(i).folder,roidata_dir(i).name)); % load roidata struct
        % check which roi boundaries to import
        tmp_roi_name = rec_cell{idx,1};
        if roidata.n_rois > numel(tmp_roi_name)
            roi_idx = [];
            for ii = 1:numel(tmp_roi_name)
                tmp_idx = split(tmp_roi_name{ii},'_');
                roi_idx = [roi_idx str2num(tmp_idx{end})];
                %             roi_idx(roi_idx>roidata.n_rois) = [];
            end
        else, roi_idx = 1:roidata.n_rois;
        end
        roi_cell{idx,1} = roi_idx;
        roi_cell{idx,2} = roidata.roi_bounds(roi_idx,1); % [y x] (0,0 in top left corner)
    end
end

%% 4) Calculate Rel. position of ROIs to reference
for iRec = 1:nrecs
    % Resolution of FOV
    fov_res = rec_cell{iRec,4};
    fov_w = rec_cell{iRec,3}*fov_res;
    % Resolution of reference
    ref_res = ref_cell{rec_cell{iRec,8},6};
    ref_w = ref_cell{rec_cell{iRec,8},5}*ref_res; 
    % ROI centers
    tmp_ctr = rec_cell{iRec,2};
    % Set ctr (im um) relative to center of FOV
    tmp_ctr(:,1) = tmp_ctr(:,1) - 0.5*fov_w;
    tmp_ctr(:,2) = tmp_ctr(:,2) - 0.5*fov_w;
    % Set ctr (im um) relative to center of reference
    tmp_ctr(:,1) = 0.5*ref_w + rec_cell{iRec,5} + tmp_ctr(:,1);
    tmp_ctr(:,2) = 0.5*ref_w + rec_cell{iRec,6} + tmp_ctr(:,2);
    rec_cell{iRec,9} = tmp_ctr;
    
    % Limits (4 corner points)
    tmpfovx = 0.5*ref_w + tmp_ctr(:,1);
    tmpfovy = 0.5*ref_w + tmp_ctr(:,2);
    rec_cell{iRec,10} = [tmpfovx-fov_w/2 tmpfovy-fov_w/2;...
        tmpfovx-fov_w/2 tmpfovy+fov_w/2;...
        tmpfovx+fov_w/2 tmpfovy-fov_w/2;...
        tmpfovx+fov_w/2 tmpfovy+fov_w/2];
    
    % ROI boundaries
    tmp_bounds = roi_cell{iRec,2};
    for ii = 1:numel(tmp_bounds)
        bounds_x = tmp_bounds{ii}{1,1}(:,2);
        bounds_y = tmp_bounds{ii}{1,1}(:,1);
        % Set boundary (im um) relative to center of FOV
        bounds_x = bounds_x - 0.5*rec_cell{iRec,3}; % x
        bounds_y = bounds_y - 0.5*rec_cell{iRec,3}; % y
        % Set boundary (im um) relative to center of reference
        bounds_x = 0.5*ref_cell{rec_cell{iRec,8},5} + rec_cell{iRec,5}/ref_res + bounds_x.*(fov_res/ref_res); % x
        bounds_y = 0.5*ref_cell{rec_cell{iRec,8},5} + rec_cell{iRec,6}/ref_res + bounds_y.*(fov_res/ref_res); % y
        % Rescale to px
%         bounds_x = bounds_x;
%         bounds_y = bounds_y;
        
        tmp_bounds{ii}{1,1}(:,2) = bounds_x;
        tmp_bounds{ii}{1,1}(:,1) = bounds_y;
    end
    roi_cell{iRec,3} = tmp_bounds;
end

%% 5) Calculate relative ML layer position
% Fit line through somata & Calculate ML extension
for iRef = 1:refcnt
    tbl.x = ref_cell{iRef,3}(:,1);
    tbl.y =  ref_cell{iRef,3}(:,2);
    lnfit = fitlm(tbl.x, tbl.y, 'linear');
    coeffs = lnfit.Coefficients.Estimate;
    numcoeffs = numel(coeffs);
    c1 = coeffs(2); % X Coeff.
    c2 = 1; % Y Coeff.
    c3 = coeffs(1); % Intercept
    ref_cell{iRef,7} = [c1 c2 c3];
    ref_cell{iRef,8} = mean(ref_cell{iRef,4});
end

rec_out_idx = []; roicnt=0;
for iRec = 1:nrecs
        mlext = ref_cell{rec_cell{iRec,8},8};
        % Calculate distance of ROI centroids in ML of Cb from somata
        c1 = ref_cell{rec_cell{iRec,8},7}(1);
        c2 = ref_cell{rec_cell{iRec,8},7}(2);
        c3 = ref_cell{rec_cell{iRec,8},7}(3);
        numroictr = numel(rec_cell{iRec,1});
        roidistances = NaN(numroictr,1);
        for i = 1:numroictr
            roilocs_x = rec_cell{iRec,9}(:,1);
            roilocs_y = rec_cell{iRec,9}(:,2);
            roidistances(i) = abs(c1*roilocs_x(i) - c2*roilocs_y(i) + c3)/ sqrt(c1^2 + c2^2);
        end
        reldistances = roidistances.*100./mlext;
        rec_cell{iRec,11} = [roidistances repelem(mlext,numroictr)' reldistances];
        
        % Corner points
        limdistances = NaN(4,1);
        tmplim = rec_cell{iRec,10};
        for i = 1:4,limdistances(i) = abs(c1*tmplim(i,1) - c2*tmplim(i,2) + c3)/ sqrt(c1^2 + c2^2);end
        rellimdistances = limdistances.*100./mlext;
        rec_cell{iRec,12} = [min(limdistances) max(limdistances) min(rellimdistances) max(rellimdistances)];
        
        rec_out_idx = [rec_out_idx iRec];
        roicnt = roicnt + numroictr;
end

%% Save information
% ML position
outputpath = fullfile(info_dir_path,'ML_pos_Master_ROIs.xls'); if isa(outputpath,'cell'), outputpath=outputpath{1};end

% output_cell = cell(numel(rec_out_idx),8);
% for i = 1:numel(rec_out_idx)
%     tmpidx = rec_out_idx(i);
%     output_cell{i,1} = rec_names(tmpidx);
%     output_cell{i,2} = mean(rec_cell{tmpidx,11}(:,1));
%     output_cell{i,3} = mean(rec_cell{tmpidx,11}(:,2));
%     output_cell{i,4} = mean(rec_cell{tmpidx,11}(:,3));
%     output_cell{i,5} = rec_cell{tmpidx,12}(1);
%     output_cell{i,6} = rec_cell{tmpidx,12}(2);
%     output_cell{i,7} = rec_cell{tmpidx,12}(3);
%     output_cell{i,8} = rec_cell{tmpidx,12}(4);
% end
% 
% output_tbl = renamevars(cell2table(output_cell),{'output_cell1',...
%     'output_cell2', 'output_cell3', 'output_cell4', 'output_cell5', 'output_cell6', 'output_cell7', 'output_cell8'},...
%     {'Recording', 'planar_dist_soma_um', 'ML_thickness_um', 'rel_ML_pos', 'FOV_low_bound_um', 'FOV_up_bound_um',...
%     'FOV_rel_low_bound', 'FOV_rel_up_bound'});
% writetable(output_tbl,outputpath);

% Output information for individual ROIs
output_cell = cell(numel(roicnt),11);
cnter = 1;
for i = 1:numel(rec_out_idx)
    tmpidx = rec_out_idx(i);
    for ii = 1:size(rec_cell{tmpidx,11},1)
        output_cell{cnter,1} = rec_names(tmpidx);
        output_cell{cnter,2} = rec_cell{tmpidx,1}(ii,1);
        output_cell{cnter,3} = rec_cell{tmpidx,11}(ii,1);
        output_cell{cnter,4} = rec_cell{tmpidx,11}(ii,2);
        output_cell{cnter,5} = rec_cell{tmpidx,11}(ii,3);
        output_cell{cnter,6} = rec_cell{tmpidx,13}(ii,1);
        output_cell{cnter,7} = rec_cell{tmpidx,13}(ii,2);
        output_cell{cnter,8} = rec_cell{tmpidx,13}(ii,3);
        output_cell{cnter,9} = rec_cell{tmpidx,13}(ii,4);
        output_cell{cnter,10} = rec_cell{tmpidx,13}(ii,5);
        output_cell{cnter,11} = rec_cell{tmpidx,13}(ii,6);
        cnter = cnter+1;
    end
end

output_tbl = renamevars(cell2table(output_cell),{'output_cell1',...
    'output_cell2', 'output_cell3', 'output_cell4', 'output_cell5', 'output_cell6',...
    'output_cell7', 'output_cell8', 'output_cell9','output_cell10',  'output_cell11'},...
    {'Recording', 'ROI','planar_dist_soma_um', 'ML_thickness_um', 'rel_ML_pos',...
    'Num_Events', 'dFoF_Amp','IEI', 'EvRate', 'Area', 'Norm_Area'});
writetable(output_tbl,outputpath);

% %% Save ROIs in reference frame as .roi files
% ref_pointer = [rec_cell{:,8}];
% roi_paths = cell(refcnt,1);
% newpath = strcat(info_dir_path(1:end-4),'ROIs\ROI_set_overview\');
% if ~isfolder(newpath), mkdir(newpath); end
% 
% for iRef = 1:refcnt
%     tmp_pter = find(ref_pointer==iRef);
%     if ~isempty(tmp_pter)
%         roi_paths{iRef} = strcat(newpath, 'ROIs_',ref_cell{iRef,1}(1:end-4),'\');
%         if ~isfolder(roi_paths{iRef}), mkdir(roi_paths{iRef}); end
%         roi_names_cell = cell(100,1); ncnter = 1;
%         % Loop through recordings with common reference
%         for i = 1:numel(tmp_pter)
%             tmp_rois = roi_cell{tmp_pter(i),3};
%             tmpnm = rec_names{tmp_pter(i),1};
%             tmpnm = split(tmpnm, 'sl');
%             % Loop through ROIs
%             for ii = 1:numel(tmp_rois)
%                 roi_name = strcat('sl',tmpnm{2}, '_', num2str(roi_cell{tmp_pter(i),1}(ii)));
%                 if isa(roi_name,'cell'), roi_name=roi_name{1};end
%                 writeImageJROI(tmp_rois{ii}{1,1}, 3, 0,0,0, roi_paths{iRef}, roi_name);
%                 roi_names_cell{ncnter} = strcat(roi_paths{iRef}, roi_name,'.roi');
%                 ncnter = ncnter+1;
%             end
%         end
%         % Zip .roi files
%         roi_names_cell(ncnter:end) = [];
%         roi_zip = strcat(roi_paths{iRef}(1:end-1),'.zip'); if isa(roi_zip,'cell'), roi_zip=roi_zip{1};end
%         zip(roi_zip,roi_names_cell);
%         rmdir(roi_paths{iRef}(1:end-1), 's');
%     end
% end


% %% **** For recordings without activity: extend of FOV
% tmpinfo = NaN(1,6);
% somainfo = NaN(1,4);
% 
% % Resolution of reference
% ref_res = tmpinfo(1,6);
% ref_w = tmpinfo(1,5)*ref_res;
% 
% % Fit line through somata & Calculate ML extension
% tbl.x = somainfo(somainfo(:,3) == 0,1);
% tbl.y =  somainfo(somainfo(:,3) == 0,2);
% lnfit = fitlm(tbl.x, tbl.y, 'linear');
% coeffs = lnfit.Coefficients.Estimate;
% numcoeffs = numel(coeffs);
% c1 = coeffs(2); % X Coeff.
% c2 = 1; % Y Coeff.
% c3 = coeffs(1); % Intercept
% mlext = mean(somainfo(somainfo(:,3) ~= 0,4));
% 
% out_mtrx = NaN(size(tmpinfo,1),5);
% 
% for iF = 1:size(tmpinfo,1)
%     % Resolution of FOV
%     fov_res = tmpinfo(iF,2);
%     fov_w = tmpinfo(iF,1)*fov_res;
%     
%     % Limits (4 corner points)
%     tmpx = 0.5*ref_w +tmpinfo(iF,3);
%     tmpy = 0.5*ref_w +tmpinfo(iF,4);
%     lim = [tmpx-fov_w/2 tmpy-fov_w/2;...
%         tmpx-fov_w/2 tmpy+fov_w/2;...
%         tmpx+fov_w/2 tmpy-fov_w/2;...
%         tmpx+fov_w/2 tmpy+fov_w/2];
% 
%     % Calculate distance of ROI centroids in ML of Cb from somata
%     roidistances = NaN(4,1);
%     for i = 1:4,roidistances(i) = abs(c1*lim(i,1) - c2*lim(i,2) + c3)/ sqrt(c1^2 + c2^2);end
%     reldistances = roidistances.*100./mlext;
%     out_mtrx(iF,1) = mlext;
%     out_mtrx(iF,2) = min(roidistances);
%     out_mtrx(iF,3) = max(roidistances);
%     out_mtrx(iF,4) = min(reldistances);
%     out_mtrx(iF,5) = max(reldistances);
% end
