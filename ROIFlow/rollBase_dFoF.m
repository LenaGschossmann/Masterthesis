function [dFtraces, FoFtraces, dFoFtraces] = rollBase_dFoF(lintraces,winsize,shift, mode)
%% Calculation of a rolling baseline dF/F stack
% Baseline is of length winsize and ends at t-shift if t is the position for
% which the dF/F value is calculated
% mode options:
% 'grand': normalize dF by average across all px and timepoints
% 'px': normalize dF by pixelwise average across all timepoints
% 'roll': normalize dF by rolling baseline average

lintraces = double(lintraces);
n_frames = size(lintraces,2);
avtraces = lintraces;
avsteps = NaN(n_frames,2);
avsteps(1:winsize+shift,1) = 1; avsteps(1:winsize+shift,2) = winsize;
avsteps(winsize+shift+1:n_frames,1) = 2:n_frames-winsize-shift+1;
avsteps(winsize+shift+1:n_frames,2) = winsize+1:n_frames-shift;

for iAv = 1:n_frames
    avtraces(:,iAv) = mean(lintraces(:,avsteps(iAv,1):avsteps(iAv,2)),2);
end

dFtraces = lintraces-avtraces;

if strcmp(mode,'grand')
    px_av_val = mean(lintraces,'all');
    dFoFtraces = dFtraces./px_av_val;
    FoFtraces = lintraces./px_av_val;
elseif strcmp(mode,'px')
    px_av_val = mean(lintraces,2);
    dFoFtraces = dFtraces./px_av_val;
    FoFtraces = lintraces./px_av_val;
else
    dFoFtraces = dFtraces./avtraces;
    FoFtraces = lintraces./avtraces;
end
end