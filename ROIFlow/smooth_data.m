function [smoothed] = smooth_data(data, kernel)
% Smooth input data in moving window style;
% mean of fields covered by window is filled into middle position of window
% Input:
% data: 1D vector of data to be smoothed
% win: Integer defining the size of window for smoothing
% Output:
% smoothed: 1D vector of smoothed data

% Declare globally shared variables
% 
% if win > 0
%     hw = round(0.5*win/FTIME);
%     iPos = hw+1;
%     smoothed = data;
%     
%     while iPos < numel(data)- hw+1
%         smoothed(iPos) = mean(data(iPos-hw:iPos+hw));
%         iPos = iPos+1;
%     end
%     smoothed(1:hw) = mean(data(1:hw));
%     smoothed(end-hw+1:end) = mean(data(end-hw+1:end));
% else
%     smoothed = data;
% end

if ~isempty(kernel)
    win = numel(kernel);
    if mod(win,2) == 0, win = win+1; end
    data = double(data);
    n_frames = size(data,2);
    startpos = round(win/2);
    stoppos = n_frames-round(win/2)+1;
    smoothed = data;
    smsteps = NaN(n_frames,2);
    smsteps(1:startpos,1) = 1; smsteps(1:startpos,2) = win;
    smsteps(startpos+1:stoppos,1) = 2:n_frames-win+1;
    smsteps(startpos+1:stoppos,2) = win+1:n_frames;
    smsteps(stoppos+1:n_frames,1) = n_frames-win+1; smsteps(stoppos+1:n_frames,2) = n_frames;
    for iSm = 1:n_frames
        smoothed(:,iSm) = sum(data(:,smsteps(iSm,1):smsteps(iSm,2)) .* kernel, 2) ./ sum(kernel);
    end
else
    smoothed = data;
end



end
