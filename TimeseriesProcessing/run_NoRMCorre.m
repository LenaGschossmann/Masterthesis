%% **** Timeseries preprocessing 1 - Motion Correction ****
% Pnevmatikakis, E. A., & Giovannucci, A. (2017). NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data. Journal of neuroscience methods, 291, 83-94.

% Add paths
wd = pwd(); p = split(wd,'\'); p = p(1:end-1); pp = []; for ii=1:numel(p), pp = fullfile(pp,p{ii}); end
pp1=fullfile(pp,'OpenBFiles'); addpath(genpath(pp1)), clear('p','pp1');
pp2=fullfile(pp,'NoRMCCorre'); addpath(genpath(pp2)), clear('p','pp2');
pp3=fullfile(pp,'CaImAn_scripts'); addpath(genpath(pp3)), clear('p','pp3');

%% Batch selection
% Select .nd2 file from NikonM
[FileName, PathName] = uigetfile('*.lsm*', [],'G:\Lena\Masterthesis\In vivo\tmp', 'Multiselect', 'on');
% [FileName, PathName] = uigetfile('*.tif*', [],'G:\Katja\iglusnfrKatja', 'Multiselect', 'on');
if ~isa(FileName, 'cell'), FileName = {FileName}; end
files = FileName;

%% Motion correction
% Readout and cropping
num2read = [] ;                                  %% num of frames to be read
sframe = 1;                                      %% starting frame
numchan = 1;                                     %% num of channels in TIFF file
crop = [10 10 10 10];

% Set some parameters
K = 10;                                         % number of components to be found
tau = 1;                                        % std of gaussian kernel (size of neuron in pixels) [was []]
p = 2;                                          % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = .8;                                 % merging threshold of overlapping components
refine = false;                                 % Manually refine components ???
options = CNMFSetParms(...                          'd1',d1,'d2',d2,...                           % dimensions of datasets
    'init_method','HALS',...                      % initialize algorithm with plain NMF
    'max_iter_hals_in',50,...                     % maximum number of iterations
    'search_method','dilate',...                  % search locations when updating spatial components
    'temporal_iter',2,...                         % number of block-coordinate descent steps
    'merge_thr',merge_thr,...                           % merging threshold
    'conn_comp',false,...                         % do not limit to largest connected component for each found component
    'maxthr',0.05,...                              % for every component set pixels that are below max_thr*max_value to 0
    'plot_df',0 ... 5
    );
options.numchan  = numchan;
options.sframe   = sframe;
options.num2read = num2read;
options.crop = crop;
options.tau = tau;
options.fs = 10;                                 % Set frequency of timeseries

% NoRMCorre
tic;
for i = 1:size(files, 2)                         % not included
    close all                                   % reading files
    disp(['Working on ' files{i}])
    nam = [PathName files{i}];
    [Data,~] = readdata(nam,options);
    
    options_nonrigid = NoRMCorreSetParms('d1',size(Data,1),'d2',size(Data,2),...
        'grid_size',[32,32],'mot_uf',4,'bin_width',200,...
        'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',1000); % init_batch: 200
    
    [Data,shifts,template,options_nonrigid] = normcorre_batch(Data,options_nonrigid); %channel used to perform the movement correction
    
    % Apply shifts to the red channel (channel to which movement correction withh be applied)
    Data = apply_shifts(Data,shifts,options_nonrigid);
    
    % Save
    tmpName = regexp(files{i}, '\', 'split');
    tmpName = tmpName{end};
    tmpName = tmpName(1:end-4);
    name = strcat(tmpName, 'MoCor','.tif');
    name = fullfile(PathName, name);
    imwrite(uint16(Data(:, :, 1)), name, 'tif', 'WriteMode', 'overwrite');
    for j = 2: size(Data,3)
        imwrite(uint16(Data(:, :, j)), name, 'tif', 'WriteMode', 'append');
    end
end
toc