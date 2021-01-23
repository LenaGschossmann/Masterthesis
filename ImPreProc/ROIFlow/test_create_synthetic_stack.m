%% ***** Create synthetic tetst timeseries *****

%% Load stack
[f,p] = uigtfile({'*.tif*'; '*.lsm'; '*.nd2'}, 'Select file for test stack','Multiselect', 'off');
filepointer = fullfile(p,f);
data = read_timeseries(filepointer, 'load');
n_px = size(data,1)*size(data,2);
n_frames = size(data,3);
n_elem = n_frames*n_px;

%% Test stack 1: AIP + multiplicative gauss noise
aip = mean(data,3);
mu =
sigma = 
syndata2 = randn(size(1),size(2),size(3));

%% Test stack 2: Scrambled stack
rnd_data = reshape(data,[n_elem 1]);
rnd_data = rnd_data(randperm(n_elem));
syndata2 = reshape(rnd_data,[size(data,1), size(data,2), n_frames]);

% for idx = 1:size(dFoFim,3)
%     imwrite(uint16(dFoFim(:,:,idx)), destpath, 'tiff', 'WriteMode', 'append');
% end