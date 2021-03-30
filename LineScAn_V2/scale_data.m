function [scaled] = scale_data(input, newscale)
% Linearly transform scaling of input matrix values into a new scale
% Input:
% input: 2D matrix
% newscale: [min max] new scale
% Output:
% scaled: 2D matrix

newdelta = diff(newscale);
oldscale = [min(input,[],'all') max(input,[],'all')];
input = input-oldscale(1); % Shift and set minimum to zero
oldscale = [0 max(input,[],'all')];
scaled = input.*(100/oldscale(2)); %...express as % of max value
scaled = newscale(1) + scaled.*(newdelta/100);
end