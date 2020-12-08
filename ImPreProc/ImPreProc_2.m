%% **** Timeseries preprocessing 2 - Delta F series ****

%% Parameters
% BASEPOINTS = 10; % Frames used for baseline
winsize_s = 0.5;
frametime = 0.051;

%% Batch selection (.tiff multipage file)
% Select .nd2 file from NikonM
[files, paths] = uigetfile('*.tif*', [],'Multiselect', 'on');
if ~isa(files, 'cell'), files = {files}; end
if ~isa(paths, 'cell'), paths = {paths}; end

%% Calculate Delta F stack
for iF = 1:size(files,2)
    tic;
    fullfilename = strcat(paths{1},files{iF});
    Info=imfinfo(fullfilename);
    imsize=[Info(1).Width,Info(1).Height];
    num_px = imsize(1)*imsize(2);
    nframes = length(Info);
    
    % Read timeseries
    inputim = zeros(imsize(2), imsize(1),nframes);
    for frame=1:nframes
        inputim(:,:,frame) = single(imread(fullfilename,frame));
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
    
    destpath = strcat(paths{1}, '\','dF_', files{iF});
    for idx = 1:size(deltafim,3)
        imwrite(uint16(deltafim(:,:,idx)), destpath, 'tiff', 'WriteMode', 'append');
    end
    toc
end
