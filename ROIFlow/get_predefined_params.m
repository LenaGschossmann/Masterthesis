function [params] = get_predefined_params(varargin)
%% Set predefined parameters for the automated analysis of spontaneous events in imaging data
% Possible input as name-value pair:
% ['dFoF_baseline_s', 'dFoF_baseline_frames', 'dFoF_baseline_shift',
% 'px_exclusion', 'responding_px_factor', 'connc_px_thresh', 'responding_cc_px_factor',
% 'lut_log_maxval', 'lut_log_growthrate', 'lut_log_midpoint',
% 'overlap_thresh', 'corr_threshold', ...
% 'ev_SD_factor', 'ev_SD_factor_safe']

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
% Minimum number of px that need to be connected to accept a connected component as ROI
if ~isempty(args) && any(strcmp(args(1,:),'connc_px_thresh'))
    params.connc_px_thresh = args{2,strcmp(args(1,:),'connc_px_thresh')};
else, params.connc_px_thresh = 16;
end
% Classification of a component as responding
if ~isempty(args) && any(strcmp(args(1,:),'responding_cc_px_factor'))
    params.responding_cc_px_factor = args{2,strcmp(args(1,:),'responding_cc_px_factor')};
else, params.responding_cc_px_factor = 3;
end
% LUT - Parameters for logistic function
if ~isempty(args) && any(strcmp(args(1,:),'lut_xlow'))
    params.lut_xlow = args{2,strcmp(args(1,:),'lut_xlow')};
else, params.lut_xlow = 1;
end
if ~isempty(args) && any(strcmp(args(1,:),'lut_log_minval'))
    params.lut_log_minval = args{2,strcmp(args(1,:),'lut_log_minval')};
else, params.lut_log_minval = 5;
end
if ~isempty(args) && any(strcmp(args(1,:),'lut_log_maxval'))
    params.lut_log_maxval = args{2,strcmp(args(1,:),'lut_log_maxval')};
else, params.lut_log_maxval = 6.5;
end
if ~isempty(args) && any(strcmp(args(1,:),'lut_log_growthrate'))
    params.lut_log_growthrate = args{2,strcmp(args(1,:),'lut_log_growthrate')};
else, params.lut_log_growthrate = 0.2;
end
if ~isempty(args) && any(strcmp(args(1,:),'lut_log_midpoint'))
    params.lut_log_midpoint = args{2,strcmp(args(1,:),'lut_log_midpoint')};
else, params.lut_log_midpoint = 25;
end
% Threshold above which px are included in the ROI during the gap-filling step
if ~isempty(args) && any(strcmp(args(1,:),'fill_thresh'))
    params.fill_thresh = args{2,strcmp(args(1,:),'fill_thresh')};
else, params.fill_thresh = 0.25;
end
% Threshold for combining overlapping ROIs 
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
% Factor defining multiples of SDnoise above which event is detected
if ~isempty(args) && any(strcmp(args(1,:),'ev_SD_factor'))
    params.ev_SD_factor = args{2,strcmp(args(1,:),'ev_SD_factor')};
else, params.ev_SD_factor = 2;
end
% Factor defining multiples of SD above which event is classified as safe
if ~isempty(args) && any(strcmp(args(1,:),'ev_SD_factor_safe'))
    params.ev_SD_factor_safe = args{2,strcmp(args(1,:),'ev_SD_factor_safe')};
else, params.ev_SD_factor_safe = 4;
end
% Window shown before and after cropped trace
if ~isempty(args) && any(strcmp(args(1,:),'crp_trc_pre_s'))
    params.crp_trc_pre_s = args{2,strcmp(args(1,:),'crp_trc_pre_s')};
else, params.crp_trc_pre_s = 0.4;
end
if ~isempty(args) && any(strcmp(args(1,:),'crp_trc_post_s'))
    params.crp_trc_post_s = args{2,strcmp(args(1,:),'crp_trc_post_s')};
else, params.crp_trc_post_s = 0.8;
end
end