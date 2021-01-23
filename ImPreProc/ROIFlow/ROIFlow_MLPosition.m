%% ***** Calculate ML Position *****

% ROI Ceter information from xlsx file (from previous processing)
[roifiles, roipath] = uigetfile({'SUMMARY*.xlsx';'*.xls';'*.csv'}, 'Select files with ROI information','Multiselect', 'on');
if isa(roifiles,'char'), roifiles = {roifiles}; end
nrecfiles = numel(roifiles);

% Soma file
[somafile, somapath] = uigetfile({'*.csv';'*.xlsx';'*.xls'}, 'Select file with SOMATA information','Multiselect', 'off');
soma_info = table2array(readtable(fullfile(somapath,somafile)));
somares = soma_info(end,:);
somalocs = soma_info(1:end-1,:);

% ROI Rec FoV information
goOn = false;
while ~goOn
quest = questdlg('How do you want to enter the FoV position?', 'FoV position', 'Load List', 'Manually', 'Load List');
if strncmp(quest,'Load',4)
    [navfile, navpath] = uigetfile({'*.xlsx';'*.xls';'*.csv'}, 'Select Reference File','Multiselect', 'off');
    fov_pos_mtrx = readmatrix(fullfile(navfile,navpath));
    checkxypos = num2str(fov_pos_mtrx);
    checkfig = figure();
    title('Please check if the values are right!');
    uicontrol('style','text', 'units','normalized', 'Position', [0.05 0.05 0.45 0.9], 'String', roifiles');
    uicontrol('style','text', 'units','normalized', 'Position', [0.55 0.05 0.45 0.9], 'String', checkxypos);
    contanswer = questdlg('Continue?', 'Check', 'Yes', 'GoBack');
    if strcmp(contanswer,'Yes'),goOn = true; end
else
    fov_pos_mtrx = NaN(nrecfiles,2);
    pos_info = inputdlg(roifiles, 'Navigator Pos. (X,Y)', [1 100]);
    for iFile = 1:nrecfiles
        tmp = split(pos_info{iFile},',');
        fov_pos_mtrx(iFile,1) = str2num(tmp{1});
        fov_pos_mtrx(iFile,1) = str2num(tmp{2});
    end
    goOn=true;
end
end

for iFile = 1:nrecfiles
    roi_tbl = table2struct(readtable(fullfile(roipath,roifiles{iFile})));
    nrois = size(roi_tbl,1);
    roilocs = NaN(nrois,2);
    roilocs(:,1) = [roi_tbl(:).Ctr_X_um]';
    roilocs(:,2) = [roi_tbl(:).Ctr_Y_um]';
    
    % Calculate location of ROIs in z-stack projection
    roilocs(:,1) = 0.5* somares(1) + roilocs(:,1);
    roilocs(:,2) = 0.5* somares(2) + roilocs(:,2);
    
    % Fit line through somata
    tbl = table();
    tbl.x = somalocs(:,1);
    tbl.y =  somalocs(:,2);
    lnfit = fitlm(tbl.x, tbl.y, 'linear');
    coeffs = lnfit.Coefficients.Estimate;
    numcoeffs = numel(coeffs);
    c1 = coeffs(2); % X Coeff.
    c2 = 1; % Y Coeff.
    c3 = coeffs(1); % Intercept
    
    % Calculate distance of ROI centroids in ML of Cb from somata
    roidistances = zeros(nrois,1);
    for iR = 1:nrois
        tmpctr = roilocs(iR,1:2);
        roidistances(iR) = abs(c1*tmpctr(1) - c2*tmpctr(2) + c3)/ sqrt(c1^2 + c2^2);
    end
    
    % Save Distance Information
    [~,name,~] = fileparts(fullfile(roipath,roifiles{iFile}));
    savename = strcat(roipath,'DISTSOMA',name(8:end),'.xlsx');
    if isa(savename,'cell'), savename = savename{1};end
    
    dist_tbl = struct('ROI',[],'Distance_um',[]);
    for iRoi =1:nrois
        dist_tbl(iRoi).ROI = roi_tbl(iRoi).ROI;
        dist_tbl(iRoi).Distance_um = roidistances(iRoi);
    end
    
    writetable(struct2table(dist_tbl), savename, 'Sheet', 'DistToSoma');
    
end


