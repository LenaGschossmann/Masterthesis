%% **** Timeseries preprocessing 2 - Delta F series ****

%% Add paths
wd = pwd(); p = split(wd,'\'); p = p(1:end-1); pp = []; for ii=1:numel(p), pp = fullfile(pp,p{ii}); end
pp=fullfile(pp,'OpenBFiles'); addpath(genpath(pp)), clear('p','pp');

%% Parameters
% BASEPOINTS = 10; % Frames used for baseline
winsize_s = 0.5;
frametime = 0.06144;

%% Batch selection (.tiff multipage file)
[files, paths] = uigetfile({'*.tif*'; '*.lsm'}, [],'Multiselect', 'on');
if ~isa(files, 'cell'), files = {files}; end
if ~isa(paths, 'cell'), paths = {paths}; end

%% Calculate Delta F stack
for iF = 1:size(files,2)
    tic;
    fullfilename = strcat(paths{1},files{iF});
    
    if strcmp(fullfilename(end-3:end),'.lsm')
        [inputim,~] = load_bf_file(fullfilename,false);
        imsize = [size(inputim,2) size(inputim,1)];
        num_px = imsize(1)*imsize(2);
        nframes = size(inputim,3);
        destpath = strcat(paths{1}, 'dF_', files{iF}(1:end-3),"tif");
    else
        Info=imfinfo(fullfilename);
        imsize=[Info(1).Width,Info(1).Height];
        num_px = imsize(1)*imsize(2);
        nframes = length(Info);
        destpath = strcat(paths{1}, 'dF_', files{iF});
        
        % Read timeseries
        inputim = zeros(imsize(2), imsize(1),nframes);
        for frame=1:nframes
            inputim(:,:,frame) = single(imread(fullfilename,frame));
        end
    end
    
    %% Rolling baseline
    winsize = round(winsize_s/frametime);
    startpos = winsize+1;
    deltafim = zeros(size(inputim));
    for iC = 1:imsize(1)
        for iR = 1:imsize(2)
            for iSl = startpos:nframes
                basef = mean(inputim(iR,iC,iSl-winsize:iSl-1));
                deltafim(iR,iC,iSl) = inputim(iR,iC,iSl)-basef;
            end
            % Deal with start
            basef = mean(inputim(iR,iC,1:iSl-1));
            deltafim(iR,iC,1:iSl-1) = inputim(iR,iC,1:iSl-1)-basef;
        end
    end
    
    for idx = 1:size(deltafim,3)
        imwrite(uint16(deltafim(:,:,idx)), destpath, 'tiff', 'WriteMode', 'append');
    end
    toc
end
