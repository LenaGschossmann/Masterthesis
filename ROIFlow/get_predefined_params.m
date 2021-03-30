function [params] = get_predefined_params(varargin)
%% Set predefined parameters for the automated analysis of spontaneous events in imaging data
% Possible input as name-value pair:
% ['fbase_winsize_s', 'bg_perc', 'perc_thresh', 'perc_range', 'perc_cutoff', 'perc_fac',
% 'connc_px_thresh', 'fill_thresh', 'ev_perc_thresh', 'safe_fac', 'critical_fac',
%'overlap_thresh', 'subtract_perc', 'peak_win']

params = struct();

if numel(varargin) > 0
    args = reshape(varargin, [2 numel(varargin)/2]);
else
    args = [];
end

%% ROI detection parameters
% Rolling window for the calculation of baseline fluorescence for dF calculation
if ~isempty(args) && any(strncmp(args(1,:),'fb', 2))
    params.fbase_winsize_s = args{2,strncmp(args(1,:),'fb', 2)};
else, params.fbase_winsize_s = 0.5;
end
% Initial value for setting the background threshold
if ~isempty(args) && any(strncmp(args(1,:),'bg',2))
    params.bg_perc = args{2,strncmp(args(1,:),'bg',2)};
else, params.bg_perc = 0.25;
end
% Percentile of background values used for background correction
if ~isempty(args) && any(strncmp(args(1,:),'su',2))
    params.subtract_perc = args{2,strncmp(args(1,:),'su',2)};
else, params.subtract_perc = 0.25;
end
% Percentile and factor by which it is multiplied for the classification of a px as responding
if ~isempty(args) && any(strcmp(args(1,:),'perc_thresh'))
    params.perc_thresh = args{2,strcmp(args(1,:),'perc_thresh')};
else, params.perc_thresh = 95;
end
if ~isempty(args) && any(strcmp(args(1,:),'perc_fac'))
    params.perc_fac = args{2,strcmp(args(1,:),'perc_fac')};
else, params.perc_fac = 1.15;
end
% Percentile for excluding px due to low dF/F values
if ~isempty(args) && any(strcmp(args(1,:),'dFoF_prctile_1'))
    params.dFoF_prctile_1 = args{2,strcmp(args(1,:),'dFoF_prctile_1')};
else, params.dFoF_prctile_1 = 75;
end
if ~isempty(args) && any(strcmp(args(1,:),'dFoF_prctile_2'))
    params.dFoF_prctile_2 = args{2,strcmp(args(1,:),'dFoF_prctile_2')};
else, params.dFoF_prctile_2 = 50;
end
% % Range tested for percentile threshold value (defined by cut-off point in
% % cummulative distribution ab connected px)
% if ~isempty(args) && any(strcmp(args(1,:),'perc_range'))
%     params.perc_range = args{2,strcmp(args(1,:),'perc_range')};
% else, params.perc_range = [80 95];
% end
% % Cutoff to find appropriate percentile value
% if ~isempty(args) && any(strcmp(args(1,:),'perc_cutoff'))
%     params.perc_cutoff = args{2,strcmp(args(1,:),'perc_cutoff')};
% else, params.perc_cutoff = 0.99;
% end
% Minimum number of px that need to be connected to accept a connected component as ROI
if ~isempty(args) && any(strncmp(args(1,:),'co',2))
    params.connc_px_thresh = args{2,strncmp(args(1,:),'co', 2)};
else, params.connc_px_thresh = 10;
end
% Threshold above which px are included in the ROI during the gap-filling step
if ~isempty(args) && any(strncmp(args(1,:),'fi',2))
    params.fill_thresh = args{2,strncmp(args(1,:),'fi', 2)};
else, params.fill_thresh = 0.2;
end
% Threshold for combining overlapping ROIs 
if ~isempty(args) && any(strncmp(args(1,:),'ov',2))
    params.overlap_thresh = args{2,strncmp(args(1,:),'ov', 2)};
else, params.overlap_thresh = 0.5;
end

%% Event detection parameters
% Window for detecting peak after event onset 
if ~isempty(args) && any(strncmp(args(1,:),'peak',4))
    params.peak_win = args{2,strncmp(args(1,:),'peak', 4)};
else, params.peak_win_ms = 150;
end
% Percentile and factor for detecting events in ROI average trace 
if ~isempty(args) && any(strcmp(args(1,:),'ev_perc_thresh'))
    params.ev_perc_thresh = args{2,strncmp(args(1,:),'ev', 2)};
else, params.ev_perc_thresh = 97;
end
% if ~isempty(args) && any(strcmp(args(1,:),'ev_fac'))
%     params.ev_fac = args{2,strcmp(args(1,:),'ev_fac')};
% else, params.ev_fac = 3;
% end
% Factor defining multiples of max. neg amplitude above which event is classified as safe
if ~isempty(args) && any(strncmp(args(1,:),'sa',2))
    params.safe_fac = args{2,strncmp(args(1,:),'sa', 2)};
else, params.safe_fac = 1.5;
end
% Factor defining multiples of prctile above which event is detected
if ~isempty(args) && any(strncmp(args(1,:),'cr',2))
    params.critical_fac = args{2,strncmp(args(1,:),'cr', 2)};
else, params.critical_fac = 1.3;
end

end