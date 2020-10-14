function [deltaFovF] = deltaFovF_linescan(input)
% Calculate delta F over F value for a 2D matrix (the baseline F is taken
% as overall mean of input matrix)
% Input:
% input: 2D matrix
% Output:
% deltaFovF: 2D matrix (same size as input matrix) where each value is
% given as normalized deviation from overall mean

fix_avF = mean(input,1);
deltaFovF = input-fix_avF;
deltaFovF = deltaFovF./fix_avF;
end