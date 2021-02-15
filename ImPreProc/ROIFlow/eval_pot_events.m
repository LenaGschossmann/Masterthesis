function evpter = eval_pot_events(tmpdFoF, onset_idx, peak_idx, to_test, peakwin, template, k, corr_thresh)

%% Set variables
ctrl_runs = 100;
template_def = 2;
n_frames = size(tmpdFoF,2);
% corr_thresh = 0.9;
test_idx = find(to_test == 1);
evpter = false(1, numel(onset_idx));

%% Control - scramble trace
ctrl_matched = zeros(ctrl_runs,1);
for iR = 1:ctrl_runs
    scrmb_idx = randperm(n_frames);
    scrmb_trc = tmpdFoF(scrmb_idx);
    ctrl_corr_trace = get_corr_trace(scrmb_trc, template);
    ctrl_matched(iR) = sum(ctrl_corr_trace.^k > corr_thresh);
end

%% Apply to trace
corr_trace = get_corr_trace(tmpdFoF, template);
% matched = corr_trace > corr_thresh;
% 
% % Check if number outside of random distrib
% if sum(matched) > prctile(ctrl_matched,95)
%     corr_matched_onset = find((corr_trace.^k > corr_thresh) == 1);
%     for iEv = 1:numel(test_idx)
%         if any(onset_idx(test_idx(iEv)) == corr_matched_onset), evpter(test_idx(iEv)) = true; end
%     end
% end
% 
figure, plot(tmpdFoF,'r');
for iN =1:numel(peak_idx), xline(peak_idx(idx(iN)), 'r'); end
figure, plot(corr_trace.^k);
for iN =1:numel(peak_idx), xline(peak_idx(idx(iN)), 'r'); end
% try, for iN =1:numel(corr_matched_onset), xline(corr_matched_onset(iN), 'k'); end, end

% evpter = evpter';



%% Local functions
    function corr_trace = get_corr_trace(trace, template)
        corr_trace = NaN(1,n_frames);
        for i = 1:n_frames-numel(template)+1
            sig_1 = tmpdFoF(i:i+numel(template)-1)-min(tmpdFoF(i:i+numel(template)-1));
%             sig_1 = trace(i:i+numel(template)-1);
            sig_2 = template;
            n = numel(template);
            corr_trace(i) = dot(sig_1,sig_2);
%             corr_trace(i) = (n * sum(sig_1.*sig_2) - sum(sig_1)*sum(sig_2)) /...
%                 sqrt((n*sum(sig_1.^2) - sum(sig_1)^2) * (n*sum(sig_2.^2) - sum(sig_2)^2));
        end
    end

end