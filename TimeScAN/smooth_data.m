function [smoothed] = smooth_data(data, win)
% Smooth input data in moving window style;
% mean of fields covered by window is filled into middle position of window
% Input:
% data: 1D vector of data to be smoothed
% win: Integer defining the size of window for smoothing
% Output:
% smoothed: 1D vector of smoothed data

hw = round(0.5*win);
iPos = hw+1;
smoothed = data;

while iPos < numel(data)- hw+1
    smoothed(iPos) = mean(data(iPos-hw:iPos+hw));
    iPos = iPos+1;
end
smoothed(1:hw) = mean(data(1:hw));
smoothed(end-hw+1:end) = mean(data(end-hw+1:end));
end
