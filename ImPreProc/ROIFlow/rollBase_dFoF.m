function [dFtraces, dFoFtraces] = rollBase_dFoF(lintraces,winsize,n_frames, mode)
%% Calculation of a rolling baseline dFoF stack
% mode 1: 'grand': normalize dF by average across all px and timepoints
% mode 2: 'px': normalize dF by pixelwise average across all timepoints
% mode 3: 'roll': normalize dF by rolling baseline average

lintraces = double(lintraces);
avtraces = lintraces;
avsteps = NaN(n_frames,2);
avsteps(1:winsize+1,1) = 1; avsteps(1:winsize+1,2) = winsize;
avsteps(winsize+2:n_frames,1) = 2:n_frames-winsize;
avsteps(winsize+2:n_frames,2) = winsize+1:n_frames-1;

for iAv = 1:n_frames
    avtraces(:,iAv) = mean(lintraces(:,avsteps(iAv,1):avsteps(iAv,2)),2);
end

dFtraces = lintraces-avtraces;

if strcmp(mode,'grand')
    px_av_val = mean(lintraces,'all');
    dFoFtraces = dFtraces./px_av_val;
elseif strcmp(mode,'px')
    px_av_val = mean(lintraces,2);
    dFoFtraces = dFtraces./px_av_val;
else
    dFoFtraces = dFtraces./avtraces;
end
end