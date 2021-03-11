function [eventdata, selrois, params] = run_event_detection(eventdata,roilist, dFoFtraces, FoFtraces, traces,  params, ft, auto, getui)

n_frames = size(dFoFtraces,2);
peakwin = ceil(params.peak_win_ms/(1000*ft));
% threshold_fac = params.ev_fac;
perc_threshold = params.ev_perc_thresh;
inizone = params.fbase_winsize_s;
safe_fac = params.safe_fac;
critical_fac = params.critical_fac;
ui = [];
% ui = inputdlg({'Threshold Factor:', 'Save Factor', 'Specific ROIs (sep: ,)'} ,'Set Values',[1 30], {num2str(threshold_fac), num2str(save_fac), ''});
if getui
    ui = inputdlg({'Percentile Threshold:', 'Safe factor (multiple of max. neg. value):',...
        'Critical factor (multiple of percentile):', 'Specific ROIs (sep: ,)'} ,...
        'Set Values',[1 30], {num2str(perc_threshold),num2str(safe_fac), num2str(critical_fac),''});
    if ~isempty(ui)
        perc_threshold = str2num(ui{1}); params.ev_perc_thresh = perc_threshold;
        safe_fac = str2num(ui{2}); params.safe_fac = safe_fac;
        critical_fac = str2num(ui{3}); params.critical_fac = critical_fac;
    end
end
if isempty(ui), selrois = 1:size(dFoFtraces,1);
elseif isempty(ui{4}), selrois = 1:size(dFoFtraces,1);
else
    tmprois = split(ui{4},',');
    selrois = NaN(numel(tmprois),1);
    for iRoi = 1:numel(tmprois),selrois(iRoi)=roilist{str2num(tmprois{iRoi}),2}; end
end
nrois = numel(selrois);
rand_pts = round(n_frames*0.15);

for iRoi = 1:nrois
    tmproi = selrois(iRoi);
    tmpdFoF = dFoFtraces(tmproi,:);
    tmpFoF = FoFtraces(tmproi,:);
    tmptrc = traces(tmproi,:);
    
    %% All Negativ values
    neg_dFoF_vals = abs(tmpdFoF(tmpdFoF<0));
%     perc_neg_val = prctile(neg_dFoF_vals, perc_threshold);
%     max_neg_val = max(neg_dFoF_vals);
    %         mean_neg_amp = mean(neg_dFoF_vals);
    
    %% Detect event onset and peak based on sd of dFoF values
    perc_pos_val = prctile(tmpdFoF(tmpdFoF>0), perc_threshold);
    perc_val = perc_pos_val;
    max_val = perc_pos_val*critical_fac;
    [~, num_pos_ev, onset_idx, peak_idx] = get_ev_onset(tmpdFoF, n_frames, perc_val*critical_fac, peakwin, ft, inizone);
    
    if num_pos_ev > 0        
        %% Output
        onset_binary = false(n_frames,1); onset_binary(onset_idx) = true; % Binarized vector indicating onsets
        peak_binary = false(n_frames,1); peak_binary(peak_idx) = true; % Binarized vector indicating peaks
        dFoF_amps =  tmpdFoF(peak_idx);
        FoF_amps = tmpFoF(peak_idx);
        amps = tmptrc(peak_idx);
        num_pos_ev = numel(onset_idx);
        
        %% Classsify ROIs/ events for revision
%         % Detect events in negativ range as control
%         [~, num_neg_ev, ~, neg_peak_idx] = get_ev_onset(-tmpdFoF, n_frames, perc_val*critical_fac, peakwin, ft, inizone);
%         if num_neg_ev ~= 0, neg_dFoF_amps = abs(tmpdFoF(neg_peak_idx));
%         else, neg_dFoF_amps = 0;
%         end
        
        % Draw random points to quantify noise
        [rnd_amp] = get_rnd_ev(tmpdFoF, n_frames, rand_pts, peakwin);
        
        ev_criterion1 = dFoF_amps > max_val*safe_fac;
        ev_criterion2 = dFoF_amps <= max_val*safe_fac;
        
        roi_criterion1 = sum(ev_criterion1) / sum(ev_criterion2);
        if all(~ev_criterion2) && ~all(~ev_criterion1), roi_criterion1 = 1; end
%         roi_criterion2 = sum(neg_dFoF_amps> perc_neg_val*critical_fac) / sum(ev_criterion1 | ev_criterion2);
%         roi_criterion2 = perc_val / mean(dFoF_amps(ev_criterion1));
        roi_criterion2 = prctile(tmpdFoF(tmpdFoF>0),50)/prctile(tmpdFoF(tmpdFoF>0),99); %std(tmpdFoF(tmpdFoF>0))
        roi_criterion3 = sum(tmpdFoF<-perc_val*critical_fac) / sum(ev_criterion2);
        
        if auto
            test_idx = find(ev_criterion2 == 1);
            ev_criterion3 = false(1,numel(ev_criterion2));
            % Check if difference between onset & peak is larger than between peak and peak+1
            for iN = 1:numel(test_idx)
                if abs(diff([tmpdFoF(onset_idx(test_idx(iN))), dFoF_amps(test_idx(iN))])) >...
                        abs(diff([dFoF_amps(test_idx(iN)), tmpdFoF(peak_idx(test_idx(iN))+1)]))
                    test_mask(iN) = true;
                    ev_criterion3(test_idx(iN)) = true;
                end
            end
            %             rev_idx = test_idx(~criterion21);
        end
        
        if roi_criterion1 >= 0.5 || (any(roi_criterion1) && roi_criterion2 <= 0.21 && roi_criterion3 <= 0.1)
            % Accept all events above critical threshold
            remark = 'perfect';
            save_idx = onset_idx(ev_criterion1 | ev_criterion2);
            revise_idx = [];
        elseif any(ev_criterion1) || (roi_criterion2 <= 0.23 && roi_criterion3 <= 0.3)
            remark = 'ok';
            if auto
                save_idx = onset_idx(ev_criterion1 | ev_criterion3);
                revise_idx = [];
            else
                save_idx = onset_idx(ev_criterion1);
                revise_idx = onset_idx(ev_criterion2);
            end
        else
            if auto
                remark = 'excluded';
                save_idx = [];
                revise_idx = [];
            else
                remark = 'revise';
                save_idx = onset_idx(ev_criterion1);
                revise_idx = onset_idx(ev_criterion2);
            end
        end
        
        n_ev_class1 = sum(ev_criterion1);
        n_ev_class2 = sum(ev_criterion2);
        
        %% Output
        eventdata{tmproi,1} = onset_binary';
        eventdata{tmproi,2} = peak_binary';
        eventdata{tmproi,3} = onset_idx;
        eventdata{tmproi,4} = peak_idx;
        eventdata{tmproi,5} = [dFoF_amps' FoF_amps' amps'];
        eventdata{tmproi,6} = revise_idx;
        eventdata{tmproi,7} = save_idx;
        eventdata{tmproi,8} = rnd_amp;
        eventdata{tmproi,9} = [];
        eventdata{tmproi,10} = [perc_val max_val safe_fac critical_fac];
        eventdata{tmproi,11} = remark;
        eventdata{tmproi,12} = [n_ev_class1 n_ev_class2];
    else
        peak_binary = false(n_frames,1); % Binarized vector indicating peaks
        onset_binary = false(n_frames,1);
        eventdata{tmproi,1} = onset_binary;
        eventdata{tmproi,2} = peak_binary;
        eventdata{tmproi,3} = [];
        eventdata{tmproi,4} = [];
        eventdata{tmproi,5} = [];
        eventdata{tmproi,6} = [];
        eventdata{tmproi,7} = [];
        eventdata{tmproi,8} = [];
        eventdata{tmproi,9} = [];
        eventdata{tmproi,10} = [perc_val max_val safe_fac critical_fac];
        eventdata{tmproi,11} = 'No events';
        eventdata{tmproi,12} = [];
    end
end


%% Local functions
    function [perc_val, num_ev, onset, peak] = get_ev_onset(dFoFtrc, n_frames, perc_val, peakwin, ft, inizone)
%         if isempty(perc_neg_val), sd_val = std(dFoFtrc); end
        % Onset
        onset = find(dFoFtrc > (perc_val) == 1);
        onset = onset-1;
        if ~isempty(onset)
            onset(onset.*ft <= inizone) = [];
            test_dist = diff(onset); del = find(test_dist <= peakwin); onset(del+1) = [];
        end
        num_ev = numel(onset);
        
        % Peak
        peak = zeros(1,size(onset,2)); % Vector containing the indices of peak points
        del = false(num_ev,1);
        if num_ev > 0
            for iEv = 1:num_ev
                if onset(iEv)+peakwin > n_frames, tmpwin= n_frames-onset(iEv);
                else, tmpwin = peakwin;
                end
                [~,idxP] = max(dFoFtrc(onset(iEv):onset(iEv)+tmpwin));
                peak(iEv) = onset(iEv)+idxP-1;
                % Delete event if peak too close to end of trace
                if (n_frames-peak(iEv)) < peakwin, del(iEv) = true; end
            end
        end

        peak(del) = [];
        onset(del) = [];
        num_ev = numel(onset);
    end

    function [rnd_amp] = get_rnd_ev(dFoFtrc, n_frames, rand_pts, peakwin)
        possiblepos = 1:numel(dFoFtrc)-peakwin;
        randpos = randi(numel(possiblepos), rand_pts, 1);
        randpos = possiblepos(randpos);
        rnd_amp = NaN(rand_pts,1);
        for iP = 1:rand_pts % Calculate amplitude of maxima in the specified window after event onset with randomly chosen points
            if randpos(iP)+peakwin > n_frames, tmpwin= n_frames-randpos(iP);
            else, tmpwin = peakwin;
            end
            rnd_amp(iP) = max(dFoFtrc(randpos(iP):randpos(iP)+tmpwin));
        end
        rnd_amp = abs(rnd_amp);
    end

end
