function [shift_corr_traces] = bleach_correction(traces, n_frames, winsize)
% Applies simple exponential bleaching correction

% traces = double(traces);
% sm_traces = traces;
% % av_trace = mean(traces,1);
% % sm_traces = av_trace;
% 
% avsteps = NaN(n_frames,2);
% if mod(winsize,2) == 0, winsize = winsize+1;end
% halfwin = round(winsize/2);
% avsteps(1:halfwin,1) = 1; avsteps(1:halfwin,2) = halfwin:winsize;
% avsteps(halfwin+1:n_frames-halfwin,1) = 2:n_frames-winsize;
% avsteps(halfwin+1:n_frames-halfwin,2) = winsize+1:n_frames-1;
% avsteps(n_frames-halfwin+1:n_frames,1) = n_frames-winsize+1:n_frames-halfwin+1;
% avsteps(n_frames-halfwin+1:n_frames,2) = n_frames;
% for iAv = 1:n_frames
%     sm_traces(:,iAv) = mean(av_trace(:,avsteps(iAv,1):avsteps(iAv,2)),2);
% end
% 
% corr_traces = traces-sm_traces;
% shift_corr_traces = corr_traces + mean(traces(:,1:5),2);
% % shift_corr_traces = traces./(sm_traces./mean(sm_traces,2));
% % shift_corr_traces = corr_traces + mean(traces(:,1:25),2) - mean(sm_traces,2);
% 
Y = median(traces,1);
X = 1:numel(Y);
tbl = table(X', Y');
betas = [mean(Y) 0 -1];
mdl = @(b,x)b(1) + b(2)*x(:,1).^b(3);
% betas = [1 -1];
% mdl = @(b,x)b(1)*exp(b(2)*x(:,1));
% test=betas(1).*exp(betas(2).*X);
fitmdl = fitnlm(tbl, mdl, betas);
coeffs = table2array(fitmdl.Coefficients(:,1));
sm_trace = coeffs(1) + coeffs(2).* X .^coeffs(3);

shift_corr_traces = traces./(sm_trace./mean(sm_trace));

% shift_corr_traces = corr_traces;
% sm_traces = sm_trace;
% corr_traces = traces-sm_trace;
% shift_corr_traces = corr_traces + mean(traces(:,1:100),2);
end