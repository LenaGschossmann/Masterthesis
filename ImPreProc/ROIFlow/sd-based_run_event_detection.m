function [eventdata, selrois, params] = run_event_detection(eventdata,roilist, dFoFtraces, FoFtraces, traces,  params, ft, auto, getui)

n_frames = size(dFoFtraces,2);
peakwin = ceil(params.peak_win_ms/(1000*ft));
threshold_fac = params.ev_fac;
inizone = params.fbase_winsize_s;
safe_fac = params.safe_fac;
excl_prctile = params.excl_prctile;
critical_fac = params.critical_fac;
% save_fac = params.ev_save_fac;
ui = [];
% ui = inputdlg({'Threshold Factor:', 'Save Factor', 'Specific ROIs (sep: ,)'} ,'Set Values',[1 30], {num2str(threshold_fac), num2str(save_fac), ''});
if getui
    ui = inputdlg({'Threshold Factor:', 'Specific ROIs (sep: ,)'} ,'Set Values',[1 30], {num2str(threshold_fac), ''});
    if ~isempty(ui)
        threshold_fac = str2num(ui{1}); params.ev_fac = threshold_fac;
        %     save_fac = str2num(ui{2}); params.ev_save_fac = save_fac;
    end
end
if isempty(ui), selrois = 1:size(dFoFtraces,1);
elseif isempty(ui{2}), selrois = 1:size(dFoFtraces,1);
else
    tmprois = split(ui{2},',');
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
    
    %% Detect event onset and peak based on sd of dFoF values
    [sd_dFoF_val, num_pos_ev, onset_idx, peak_idx] = get_ev_onset(tmpdFoF, n_frames, [], threshold_fac, peakwin, ft, inizone);
%     [perc_neg_val, num_pos_ev, onset_idx, peak_idx] = get_ev_onset(tmpdFoF, n_frames, [], threshold_fac, peakwin, ft, inizone);
    
    if num_pos_ev > 0
        
        %% All Negativ values
        neg_dFoF_vals = abs(tmpdFoF(tmpdFoF<0));
        perc_neg_val = prctile(neg_dFoF_vals, excl_prctile);
%         mean_neg_amp = mean(neg_dFoF_vals);
        
        %% Detect events in negativ range as control
        [~, num_neg_ev, ~, neg_peak_idx] = get_ev_onset(-tmpdFoF, n_frames, sd_dFoF_val, threshold_fac, peakwin, ft, inizone);
        if num_neg_ev ~= 0
            neg_dFoF_amps = abs(tmpdFoF(neg_peak_idx));
            mean_neg_amp = mean(neg_dFoF_amps);
        else, neg_dFoF_amps = 0; mean_neg_amp = 0;
        end
        
        %% Delete any events with dFoF amplitude < perc_neg_val
%         tmp_amp = tmpdFoF(peak_idx);
%         del_ev = tmp_amp <= perc_neg_val;
%         onset_idx = onset_idx(~del_ev);
%         peak_idx = peak_idx(~del_ev);
        
        %% Draw random points to quantify noise
        [rnd_amp] = get_rnd_ev(tmpdFoF, n_frames, rand_pts, peakwin);
        
        %% Output
        onset_binary = false(n_frames,1); onset_binary(onset_idx) = true; % Binarized vector indicating onsets
        peak_binary = false(n_frames,1); peak_binary(peak_idx) = true; % Binarized vector indicating peaks
        dFoF_amps =  tmpdFoF(peak_idx);
        FoF_amps = tmpFoF(peak_idx);
        amps = tmptrc(peak_idx);
        min_dFoF_amp = min(dFoF_amps);
        num_pos_ev = numel(onset_idx);
        
        %% Classsify ROIs/ events for revision
        if auto
%             criterion1 = dFoF_amps > mean_neg_amp*safe_fac;
%             criterion2 = dFoF_amps <= mean_neg_amp*safe_fac & dFoF_amps > perc_neg_val*critical_fac;
            criterion1 = dFoF_amps > perc_neg_val*safe_fac;
            criterion2 = dFoF_amps <= perc_neg_val*safe_fac; % & dFoF_amps > perc_neg_val*critical_fac;
%             test_idx = find(criterion2 == 1);
%             test_n = numel(test_idx);
%             criterion21 = false(test_n,1);
%             % Check if difference between onset & peak is larger than
%             % between peak and peak+1
%             for iN = 1:test_n
%                 if abs(diff([tmpdFoF(onset_idx(test_idx(iN))), dFoF_amps(test_idx(iN))])) >...
%                         abs(diff([dFoF_amps(test_idx(iN)), tmpdFoF(peak_idx(test_idx(iN))+1)]))
%                     criterion21(iN) = true;
%                 end
%             end
% %             rev_idx = test_idx(~criterion21);
%             test_idx = test_idx(criterion21);
% %             test_amp = dFoF_amps(test_idx);
% %             test_amp = test_amp(criterion3);
%             %             [~,id] = sort(test_amp, 'descend');
%             %             nacc = round(numel(test_idx)/2);
%             %             keep_idx = test_idx(id(1:nacc));
%             criterion2 = false(1,numel(criterion1)); criterion2(test_idx) = true;
            save_idx = onset_idx(criterion1 | criterion2);
%             revise_idx = onset_idx(criterion2);
%             revise_idx = onset_idx(rev_idx);
            revise_idx = [];
            remark = 'auto';
            n_ev_class1 = sum(criterion1);
            n_ev_class2 = sum(criterion2);
        else
            perc_neg_val = NaN;
            factor1 = sum(neg_dFoF_amps>=min_dFoF_amp)/num_pos_ev;
            factor2 = sd_dFoF_val/mean(dFoF_amps);
            save_idx = onset_idx;
            if factor1 <= 0.15 && factor2 <= 0.15
                remark = 'perfect';
                revise_idx = [];
            elseif factor1 <= 0.25  && factor2 <= 0.25
                remark = 'good';
                % Mark single events as save
                revise_idx = onset_idx(dFoF_amps <= mean_neg_amp*safe_fac);
                %             revise_idx = onset_idx(dFoF_amps <= save_fac*sd_dFoF_val);
                for iRev = 1:numel(revise_idx), save_idx(revise_idx(iRev)==save_idx) = []; end
            else
                remark = 'revise';
                % Mark single events as save
                if num_neg_ev~=0, revise_idx = onset_idx(dFoF_amps <= mean_neg_amp*safe_fac);
                else,revise_idx = onset_idx;
                end
                %             revise_idx = onset_idx(dFoF_amps <= save_fac*sd_dFoF_val);
                for iRev = 1:numel(revise_idx), save_idx(revise_idx(iRev)==save_idx) = []; end
            end
            n_ev_class1 = numel(save_idx);
            n_ev_class2 = numel(revise_idx);
        end
        
        %% Output
        eventdata{tmproi,1} = onset_binary';
        eventdata{tmproi,2} = peak_binary';
        eventdata{tmproi,3} = onset_idx;
        eventdata{tmproi,4} = peak_idx;
        eventdata{tmproi,5} = [dFoF_amps' FoF_amps' amps'];
        eventdata{tmproi,6} = revise_idx;
        eventdata{tmproi,7} = save_idx;
        eventdata{tmproi,8} = rnd_amp;
        eventdata{tmproi,9} = neg_dFoF_vals';
%         eventdata{tmproi,10} = [threshold_fac mean_neg_amp*safe_fac perc_neg_val];
        eventdata{tmproi,10} = [sd_dFoF_val*threshold_fac perc_neg_val*safe_fac perc_neg_val*critical_fac];
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
        eventdata{tmproi,10} = [sd_dFoF_val*threshold_fac NaN NaN];
        eventdata{tmproi,11} = 'No events';
        eventdata{tmproi,12} = [];
    end
end

%% Select Events fully automatic
% if auto
%     template = get_crosscorr_template(eventdata,dFoFtraces, peakwin);
%     for iRoi = 1:nrois
%         tmproi = selrois(iRoi);
%         tmpdFoF = dFoFtraces(tmproi,:);
%         onset_idx = eventdata{tmproi,3};
%         peak_idx = eventdata{tmproi,4};
%         criterion1 = eventdata{tmproi,5}(:,1) > eventdata{tmproi,10}(2);
%         to_test = ~criterion1;
%         if any(to_test)
%             criterion3 = eval_pot_events(tmpdFoF, onset_idx, peak_idx, to_test, peakwin, template, 1, 0.9);
%         else, criterion3 = false(numel(criterion1),1);
%         end
%         eventdata{tmproi,7} = onset_idx(criterion1 | criterion3); % considers dynamics
%         eventdata{tmproi,6} = save_idx;
%     end
% end

%% Local functions
    function [sd_val, num_ev, onset, peak] = get_ev_onset(dFoFtrc, n_frames, sd_val, threshold_fac, peakwin, ft, inizone)
        if isempty(sd_val), sd_val = std(dFoFtrc); end
        % Onset
        onset = find(dFoFtrc > (sd_val*threshold_fac) == 1);
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
                if (n_frames-peak(iEv)) < peakwin, del(iEv) = true; end
            end
        end
        
        % Delete event if peak too close to end of trace
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



% %%
% neg_onset_idx = find((tmpdFoF < (sd_dFoF_val*(-threshold_fac))) == 1);
%     neg_onset_idx = neg_onset_idx-1;
%     if numel(neg_onset_idx) > 0
%         test_dist = diff(neg_onset_idx); del = find(test_dist <= peakwin); neg_onset_idx(del+1) = [];
%         num_neg_ev = numel(neg_onset_idx);
%         neg_dFoF_amps = NaN(num_neg_ev,1);
%         for iP = 1:num_neg_ev % Calculate amplitude of maxima in the specified window after event onset with randomly chosen points
%             if neg_onset_idx(iP)+peakwin > n_frames, tmpwin= n_frames-neg_onset_idx(iP);
%             else, tmpwin = peakwin;
%             end
%             neg_dFoF_amps(iP) = abs(min(tmpdFoF(neg_onset_idx(iP):neg_onset_idx(iP)+tmpwin)));
%         end
%         mean_neg_amp =  mean(neg_dFoF_amps);
%     else
%         mean_neg_amp = 0;
%     end