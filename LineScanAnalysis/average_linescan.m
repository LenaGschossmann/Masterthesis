function [binned] = average_linescan(input, wsz)
% Bin input matrix along column dimension (all fields in binning window
% will take on the windows mean value)
% Input:
% input: 2D matrix of which the columns should be binned
% wsz: binning window size
% Output:
% binned: 2D binned matrix (same size as input matrix)

tmphw = size(input);
hw = round(wsz/2);
binned = zeros(tmphw);
iAv = hw;                                                                   % Start position at half of binning window size

while iAv <= tmphw(2)-hw
    binned(:,iAv-hw+1:iAv+hw) = repmat(mean(input(:,iAv-hw+1:iAv+hw), 2), [1 wsz]);
    iAv = iAv+wsz;
end

% Fill in at start and end of each row
binned(:,1:hw) = repmat(mean(input(:,1:hw), 2), [1 hw]);
if mod(tmphw(2),wsz) ~= 0, binned(:,end-(iAv-wsz):end) = repmat(mean(input(:,end-(iAv-wsz):end), 2), [1 hw]);
else, binned(:,end-hw+1:end) = repmat(mean(input(:,end-hw+1:end), 2), [1 hw]);
end

end