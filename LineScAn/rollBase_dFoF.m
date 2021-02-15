function [dFtraces, dFoFtraces] = rollBase_dFoF(input)
%% Calculation of a rolling baseline dFoF stack
global dFWIN FTIME

winsize = round(dFWIN/FTIME);
lintraces = double(input)';
n_frames = size(lintraces,2);
avtraces = lintraces;
avsteps = NaN(n_frames,2);
avsteps(1:winsize+1,1) = 1; avsteps(1:winsize+1,2) = winsize;
avsteps(winsize+2:n_frames,1) = 2:n_frames-winsize;
avsteps(winsize+2:n_frames,2) = winsize+1:n_frames-1;

for iAv = 1:n_frames
    avtraces(:,iAv) = mean(lintraces(:,avsteps(iAv,1):avsteps(iAv,2)),2);
end

dFtraces = lintraces-avtraces;
px_av_val = mean(lintraces,2);
dFoFtraces = dFtraces./px_av_val;

dFtraces = dFtraces';
dFoFtraces = dFoFtraces';

end