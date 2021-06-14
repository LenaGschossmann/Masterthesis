function [params] = get_predefined_params(varargin)
%% Set predefined parameters for the automated analysis of spontaneous events in imaging data
% Possible input as name-value pair:
% ['dFoF_baseline_s', 'dFoF_baseline_frames', 'dFoF_baseline_shift',
% 'bg_perc', 'px_exclusion', 'responding_px_factor','connc_px_thresh', 'fill_thresh',
% 'ev_perc_thresh', 'safe_fac', 'critical_fac',
% 'overlap_thresh', 'subtract_perc', 'peak_win']

params = struct();

if numel(varargin) > 0
    args = reshape(varargin, [2 numel(varargin)/2]);
else
    args = [];
end

%% ROI detection parameters
% Rolling window for the calculation of baseline fluorescence for dF calculation
if ~isempty(args) && any(strcmp(args(1,:),'dFoF_baseline_s'))
    params.dFoF_baseline_s = args{2,strcmp(args(1,:),'dFoF_baseline_s')};
else, params.dFoF_baseline_s = 1;
end
if ~isempty(args) && any(strcmp(args(1,:),'dFoF_baseline_frames'))
    params.dFoF_baseline_frames = args{2,strcmp(args(1,:),'dFoF_baseline_frames')};
else, params.dFoF_baseline_frames = 10;
end
if ~isempty(args) && any(strcmp(args(1,:),'dFoF_baseline_shift'))
    params.dFoF_baseline_shift = args{2,strcmp(args(1,:),'dFoF_baseline_shift')};
else, params.dFoF_baseline_shift = 5;
end
% % Initial value for setting the background threshold
% if ~isempty(args) && any(strncmp(args(1,:),'bg',2))
%     params.bg_perc = args{2,strncmp(args(1,:),'bg',2)};
% else, params.bg_perc = 0.25;
% end
% Percentile of background values used for background correction
% if ~isempty(args) && any(strncmp(args(1,:),'su',2))
%     params.subtract_perc = args{2,strncmp(args(1,:),'su',2)};
% else, params.subtract_perc = 0.25;
% end
% Px exclusion
if ~isempty(args) && any(strcmp(args(1,:),'px_exclusion'))
    params.px_exclusion = args{2,strcmp(args(1,:),'px_exclusion')};
else, params.px_exclusion = 0.25;
end
% Classification of a px as responding
if ~isempty(args) && any(strcmp(args(1,:),'responding_px_factor'))
    params.responding_px_factor = args{2,strcmp(args(1,:),'responding_px_factor')};
else, params.responding_px_factor = 1.15;
end
if ~isempty(args) && any(strcmp(args(1,:),'responding_cc_px_factor'))
    params.responding_cc_px_factor = args{2,strcmp(args(1,:),'responding_cc_px_factor')};
else, params.responding_cc_px_factor = 3;
end
% Minimum number of px that need to be connected to accept a connected component as ROI
if ~isempty(args) && any(strcmp(args(1,:),'connc_px_thresh'))
    params.connc_px_thresh = args{2,strcmp(args(1,:),'connc_px_thresh')};
else, params.connc_px_thresh = 16;
end
% if ~isempty(args) && any(strcmp(args(1,:),'connc_px_thresh_high'))
%     params.connc_px_thresh_high = args{2,strcmp(args(1,:),'connc_px_thresh_high')};
% else, params.connc_px_thresh_high = 6;
% end
% Threshold above which px are included in the ROI during the gap-filling step
if ~isempty(args) && any(strncmp(args(1,:),'fi',2))
    params.fill_thresh = args{2,strncmp(args(1,:),'fi', 2)};
else, params.fill_thresh = 0.25;
end
% Threshold for combining overlapping ROIs 
% if ~isempty(args) && any(strcmp(args(1,:),'overlap_thresh'))
%     params.overlap_thresh = args{2,strcmp(args(1,:),'overlap_thresh')};
% else, params.overlap_thresh = 0.5;
% end
if ~isempty(args) && any(strcmp(args(1,:),'corr_threshold'))
    params.corr_threshold = args{2,strcmp(args(1,:),'corr_threshold')};
else, params.corr_threshold = 0.25;
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
else, params.critical_fac = 1.5;
end

end