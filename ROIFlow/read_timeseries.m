function [result] = read_timeseries(filepointer, mode)
%% Function to read timeseries file
% Depends on OpenBFiles folder
% Input:
% filepointer: (char) full filename
% mode: 'load' (to load complete stack), 'aip' (to load stack generate
% average intensity projection and clear stack from memory


%% Add paths
wd = pwd(); p = split(wd,'\'); p = p(1:end-2); pp = []; for ii=1:numel(p), pp = fullfile(pp,p{ii}); end
pp=fullfile(pp,'OpenBFiles'); addpath(genpath(pp)), clear('p','pp');

%% Read file
if strcmp(filepointer(end-3:end),'.lsm') || strcmp(filepointer(end-3:end),'.nd2')
    [inputim,~] = load_bf_file(filepointer,false);
else
    disp('Reading .tif file...');
    Info=imfinfo(filepointer);
    imsize=[Info(1).Width,Info(1).Height];
    n_frames = length(Info);
    
    % Read timeseries
    inputim = NaN(imsize(2), imsize(1),n_frames);
    parfor frame=1:n_frames
        inputim(:,:,frame) = double(imread(filepointer,frame));
    end
end

if strcmp(mode,'load')
    result = inputim;
elseif strcmp(mode,'aip')
    result =  mean(inputim,3);
    clear('inputim');
end

disp('File loaded.');

end