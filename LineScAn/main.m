%% ** Script for line scan analysis **

% Add path to directory with scripts for reading bioformat files
% wd = pwd(); p = split(wd,'\'); p=p(1:end-1); pp = []; for ii = 1:numel(p), pp = fullfile(pp, p{ii}); end
% pp = fullfile(pp,'OpenBFiles');
% addpath(pp); clear('p', 'pp');

addpath(genpath('G:\Lena\Masterthesis\In vivo\LenaMatlabCode\OpenBFiles'));

close all;

% Initiate variables
set_global_vars();
ini_ctrl_box_linescan();
