function [shift_corr_traces,sm_traces] = bleach_correction(traces, n_frames, winsize)
% Applies simple exponential bleaching correction

traces = double(traces);
sm_traces = traces;
avsteps = NaN(n_frames,2);
if mod(winsize,2) == 0, winsize = winsize+1;end
halfwin = round(winsize/2);
avsteps(1:halfwin,1) = 1; avsteps(1:halfwin,2) = halfwin:winsize;
avsteps(halfwin+1:n_frames-halfwin,1) = 2:n_frames-winsize;
avsteps(halfwin+1:n_frames-halfwin,2) = winsize+1:n_frames-1;
avsteps(n_frames-halfwin+1:n_frames,1) = n_frames-winsize+1:n_frames-halfwin+1;
avsteps(n_frames-halfwin+1:n_frames,2) = n_frames;
for iAv = 1:n_frames
    sm_traces(:,iAv) = mean(traces(:,avsteps(iAv,1):avsteps(iAv,2)),2);
end
corr_traces = traces./(sm_traces./mean(sm_traces,2));
shift_corr_traces = corr_traces + (mean(traces(:,1:10),2)-mean(corr_traces,2));

corr_traces = traces-sm_traces;
shift_corr_traces = corr_traces + mean(traces(:,1:5),2);

end