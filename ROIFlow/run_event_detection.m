function [eventdata, roi_list, params] = run_event_detection(eventdata,roilist, params, ft, auto, getui)

%% Function to detect events in optical imaging data using iGluSnFR as sensor
n_frames = size(eventdata(1).filtrd_trace,2);
rand_pts = round(n_frames*0.15);

peakwin = ceil(params.peak_win_ms/(1000*ft));
inizone = params.dFoF_baseline_s;

SD_factor_FoF_1d = params.ev_SD_factor;
SD_factor_FoF = params.ev_SD_factor_safe;

ui = [];
if getui
    ui = inputdlg({'dF/F SD factor (slope):', 'F/F - SD factor (amp):', 'Specific ROIs (sep: ,)'} ,...
        'Set Values',[1 30], {num2str(SD_factor_FoF_1d), num2str(SD_factor_FoF), ''});
    if ~isempty(ui)
        SD_factor_FoF_1d = str2num(ui{1}); params.ev_SD_factor = SD_factor_FoF_1d;
        SD_factor_FoF = str2num(ui{2}); params.ev_SD_factor_safe = SD_factor_FoF;
    end
end

if isempty(ui), roi_list = 1:size(eventdata,2);
elseif isempty(ui{3}), roi_list = 1:size(eventdata,2);
else
    tmprois = split(ui{3},',');
    roi_list = NaN(numel(tmprois),1);
    for iRoi = 1:numel(tmprois),roi_list(iRoi)=roilist{str2num(tmprois{iRoi}),2}; end
end

nrois = numel(roi_list);

%% Define SD(noise)
FoF1d_std_noise = NaN(nrois,1);
FoF_std_noise = NaN(nrois,1);
for iRoi = 1:nrois
    tmpFoF = eventdata(iRoi).filtrd_FoFtrace;
    % SD of 1st derivative of FoF < 0
    FoF_filtrd_1d = diff(tmpFoF,  1, 2); % 1st derivative
    trc_neg = FoF_filtrd_1d; trc_neg = trc_neg(trc_neg < 0);
    trc_synt_symm = [trc_neg -(trc_neg)];
    FoF1d_std_noise(iRoi,1) = std(trc_synt_symm);
    % SD of FoF < 1
    trc_neg = tmpFoF; trc_neg = trc_neg(trc_neg < 1);
    trc_synt_symm = [trc_neg 1+(1-trc_neg)];
    FoF_std_noise(iRoi,1) = std(trc_synt_symm);
end

%% ********** Detect events **********
for iRoi = 1:nrois
    tmproi = roi_list(iRoi);
    tmpFoF_1d = diff(eventdata(iRoi).filtrd_FoFtrace, 1, 2);
    tmpFoF = eventdata(iRoi).filtrd_FoFtrace;
    tmptrc = eventdata(iRoi).filtrd_trace;
    stdFoF1d = FoF1d_std_noise(iRoi);
    stdFoF = FoF_std_noise(iRoi);
    ctrl_av_peak_FoF = [];
    
    % Criterion 1 - Slope: Exceeding F/F 1st derivative SDnoise threshold
    onset_idx = find(tmpFoF_1d > (stdFoF1d * SD_factor_FoF_1d) == 1);
    if ~isempty(onset_idx)
        % Delete events in very beginning of trace
        onset_idx(onset_idx.*ft <= inizone) = [];
        % Delete onset if directly preceded by another
        del = find(diff(onset_idx) == 1);
        onset_idx(del+1) = [];
    end
    num_ev = numel(onset_idx);
    
    % Peak Index
    peak_idx = zeros(1,size(onset_idx,2)); % Vector containing the indices of peak_idx points
    del = false(num_ev,1);
    for iEv = 1:num_ev
        if onset_idx(iEv)+peakwin > n_frames, tmpwin= n_frames-onset_idx(iEv);
        else, tmpwin = peakwin;
        end
        [~,idxP] = max(tmpFoF(onset_idx(iEv):onset_idx(iEv)+tmpwin));
        peak_idx(iEv) = onset_idx(iEv)+idxP-1;
        
        % Delete event if peak_idx too close to end of trace
        if (n_frames-peak_idx(iEv)) < peakwin, del(iEv) = true; end
        
        % Criterion 2 - F/F Amplitude: Delete if event too small
        if tmpFoF(peak_idx(iEv)) <= 1+(stdFoF * SD_factor_FoF), del(iEv) = true; end
    end
    % Update
    peak_idx(del) = [];
    onset_idx(del) = [];
    num_ev = numel(onset_idx);
    
    if num_ev > 0        
        %% Controls
        % Control 1: Detect control events in neg. range using 1st deriv. of F/F
        ctrl_FoF = -(tmpFoF - 1) + 1;
        ctrl_FoF_1d = diff(ctrl_FoF);
        ctrl_onset_idx = find(ctrl_FoF_1d > (stdFoF1d * SD_factor_FoF_1d) == 1);
        if ~isempty(ctrl_onset_idx)
            ctrl_onset_idx(ctrl_onset_idx.*ft <= inizone) = [];
            % Delete onset if directly preceded by another
            del = find(diff(onset_idx) == 1);
            onset_idx(del+1) = [];
        end
        ctrl_num_ev = numel(ctrl_onset_idx);
        
        ctrl_peak_idx = zeros(1,size(ctrl_onset_idx,2)); % Vector containing the indices of peak_idx points
        critical_ctrl_ev = false(num_ev,1);
        if ctrl_num_ev > 0
            for iEv = 1:ctrl_num_ev
                if ctrl_onset_idx(iEv)+peakwin > n_frames, tmpwin= n_frames-ctrl_onset_idx(iEv);
                else, tmpwin = peakwin;
                end
                [~,idxP] = max(ctrl_FoF(ctrl_onset_idx(iEv):ctrl_onset_idx(iEv)+tmpwin));
                ctrl_peak_idx(iEv) = ctrl_onset_idx(iEv)+idxP-1;
                % Delete event if peak_idx too close to end of trace
                if (n_frames-ctrl_peak_idx(iEv)) < peakwin, critical_ctrl_ev(iEv) = true; end
%                 % Criterion 2 - F/F Amplitude: Delete if event too small
%                 if ctrl_FoF(ctrl_peak_idx(iEv)) <= 1+(stdFoF * SD_factor_FoF), critical_ctrl_ev(iEv) = true; end
            end
            ctrl_peak_FoF = ctrl_FoF(ctrl_peak_idx);
            ctrl_peak_amp_FoF = ctrl_FoF(ctrl_peak_idx) - ctrl_FoF(ctrl_onset_idx);
            if isempty(ctrl_peak_FoF)
                ctrl_av_peak_FoF = [];
            else
                ctrl_av_peak_FoF = mean(ctrl_peak_FoF);
                ctrl_max_peak_FoF = max(ctrl_peak_FoF);
            end
        end
        % Discard all events that are within amplitude range of events detected in "neg." range
        del = tmpFoF(peak_idx) <= ctrl_max_peak_FoF;
        % Update
        peak_idx(del) = [];
        onset_idx(del) = [];
        num_ev = numel(onset_idx);
        
        % Control 2: Draw random points to assess by-chance detected peak
        % F/F values
        possiblepos = 1:n_frames-peakwin;
        randpos = randi(numel(possiblepos), rand_pts, 1);
        randpos = possiblepos(randpos);
        rnd_peak_FoF = NaN(1, rand_pts);
        for iP = 1:rand_pts % Calculate amplitude of maxima in the specified window after event onset_idx with randomly chosen points
            if randpos(iP)+peakwin > n_frames, tmpwin= n_frames-randpos(iP);
            else, tmpwin = peakwin;
            end
            rnd_peak_FoF(iP) = max(tmpFoF(randpos(iP):randpos(iP)+tmpwin));
        end
        
        %% Assess certainty
        %         ev_criterion1 = tmpFoF(onset_idx) > SD_factor_FoF * stdFoF;
        %         ev_criterion2 = tmpFoF(onset_idx) <= SD_factor_FoF * stdFoF;
        
        %         roi_criterion1 = sum(ev_criterion1) / sum(ev_criterion2);
        %         if all(~ev_criterion2) && ~all(~ev_criterion1), roi_criterion1 = 1; end
        %         roi_criterion2 = prctile(tmpFoF(tmpFoF>1),50)/prctile(tmpFoF(tmpFoF>1),99); %std(tmpdFoF(tmpdFoF>0))
        %         roi_criterion3 = sum(tmpdFoF < - SD_factor_FoF_1d * tmpSDnoise) / sum(ev_criterion2);
        
        
        % Check if difference between onset_idx & peak_idx is larger than between peak_idx and peak_idx+1
        rise = abs(diff([tmpFoF(onset_idx); tmpFoF(peak_idx)],1,1));
        decay = abs(diff([tmpFoF(peak_idx); tmpFoF(peak_idx+1)],1,1));
        auto_keep_ev = rise > 1.1 * decay;
        % Let amplitude override kinetics criterion
        amp_criterion = tmpFoF(peak_idx) >= 1+(2 * stdFoF * SD_factor_FoF);
        auto_keep_ev = auto_keep_ev | amp_criterion;
        if auto
            % Update
            onset_idx = onset_idx(auto_keep_ev);
            peak_idx = peak_idx(auto_keep_ev);
            num_ev = numel(onset_idx);
        end
        
        %         if roi_criterion1 >= 0.5 || (any(roi_criterion1) && roi_criterion2 <= 0.21 && roi_criterion3 <= 0.1)
%             % Accept all events above critical threshold
%             remark = 'perfect';
%             save_idx = onset_idx(ev_criterion1 | ev_criterion2);
%             revise_idx = [];
%         elseif any(ev_criterion1) || (roi_criterion2 <= 0.23 && roi_criterion3 <= 0.3)
%             remark = 'ok';
%             if auto
%                 save_idx = onset_idx(ev_criterion1 | ev_criterion3);
%                 revise_idx = [];
%             else
%                 save_idx = onset_idx(ev_criterion1);
%                 revise_idx = onset_idx(ev_criterion2);
%             end
%         else
%             if auto
%                 remark = 'excluded';
%                 save_idx = [];
%                 revise_idx = [];
%             else
%                 remark = 'revise';
%                 save_idx = onset_idx(ev_criterion1);
%                 revise_idx = onset_idx(ev_criterion2);
%             end
%         end
        
%         n_ev_class1 = sum(ev_criterion1);
%         n_ev_class2 = sum(ev_criterion2);
        
        %% Output
%         onset_binary = false(n_frames,1); onset_binary(onset_idx) = true; % Binarized vector indicating onsets
%         peak_binary = false(n_frames,1); peak_binary(peak_idx) = true; % Binarized vector indicating peaks
%         dFoF_events = tmpdFoF(onset_idx + 1);
        FoF_peak_amps = tmpFoF(peak_idx) - tmpFoF(onset_idx);
        peak_amps = tmptrc(peak_idx) - tmptrc(onset_idx);
        eventdata(tmproi).events_detected = true;
%         eventdata(tmproi).onset_binary = onset_binary';
%         eventdata(tmproi).peak_binary = peak_binary';
        eventdata(tmproi).onset_idx = onset_idx;
        eventdata(tmproi).peak_idx = peak_idx;
        eventdata(tmproi).auto_keep_ev = auto_keep_ev;
%         eventdata(tmproi).peak_amp = peak_amps';
%         eventdata(tmproi).FoF_peak_amp = FoF_peak_amps';
        eventdata(tmproi).FoF_SDnoise = stdFoF;
        eventdata(tmproi).FoF_1stDer_SDnoise = stdFoF1d;
        eventdata(tmproi).FoF_1stDer_threshold = stdFoF1d*SD_factor_FoF_1d;
        eventdata(tmproi).FoF_threshold = 1+stdFoF*SD_factor_FoF;
        eventdata(tmproi).ctrl_peak_FoF = ctrl_peak_FoF;
        eventdata(tmproi).ctrl_av_peak_FoF = ctrl_av_peak_FoF;
        eventdata(tmproi).rnd_peak_FoF = rnd_peak_FoF;
    else
        eventdata(tmproi).events_detected = false;
%         eventdata(tmproi).onset_binary = false(n_frames,1);
%         eventdata(tmproi).peak_binary = false(n_frames,1);
        eventdata(tmproi).onset_idx = [];
        eventdata(tmproi).peak_idx = [];
        eventdata(tmproi).auto_keep_ev = [];
%         eventdata(tmproi).peak_amp = [];
%         eventdata(tmproi).FoF_peak_amp = [];
        eventdata(tmproi).FoF_SDnoise = stdFoF;
        eventdata(tmproi).FoF_1stDer_SDnoise = stdFoF1d;
        eventdata(tmproi).FoF_1stDer_threshold = stdFoF1d*SD_factor_FoF_1d;
        eventdata(tmproi).FoF_threshold = 1+stdFoF*SD_factor_FoF;
        eventdata(tmproi).ctrl_peak_FoF = [];
        eventdata(tmproi).ctrl_av_peak_FoF = [];
        eventdata(tmproi).rnd_peak_FoF = [];
    end
end

end
