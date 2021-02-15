%% ***** Create synthetic tetst timeseries *****

%% Load stack
[f,p] = uigetfile({'*.tif*'; '*.lsm'; '*.nd2'}, 'Select file for test stack','Multiselect', 'off');
filepointer = fullfile(p,f);
data = read_timeseries(filepointer, 'load');
n_px = size(data,1)*size(data,2);
n_frames = size(data,3);
n_elem = n_frames*n_px;
datatraces = reshape(data,[n_px n_frames]);
ft = 0.06144;

%% Output path
% savepath1 = 'C:\Users\lena_\Projects\iGluSnFR\Analysis\Tests\LenaMatlabCode\DD_ROI_detection\NegCtrl\noised_data_conf.tif';
savepath1 = strcat(filepointer(1:end-4),'_NOISED.tif');
savepath2 = 'C:\Users\lena_\Projects\iGluSnFR\Analysis\Tests\LenaMatlabCode\DD_ROI_detection\NegCtrl\scrambled_stack_conf.tif';

%% Test stack 1: AIP + multiplicative gauss noise
aip = mean(data,3);
aipstack = repmat(aip,[1 1 n_frames]);
% [~, dFoFtraces] = rollBase_dFoF(datatraces,round(0.5/ft),n_frames);
mu = 1;
sigma = 0.5;
noise = mu+randn(size(data,1),size(data,2),n_frames).*sigma;
syndata1 = aipstack.*noise;

n=noise-min(noise,[],'all'); n=n./max(n,[],'all'); n=n.*255;
for idx = 1:n_frames
    imwrite(uint16(syndata1(:,:,idx)), savepath1, 'tiff', 'WriteMode', 'append');
end

%% Test stack 2: Scrambled stack
rnd_data = reshape(data,[n_elem 1]);
rnd_data = rnd_data(randperm(n_elem));
syndata2 = reshape(rnd_data,[size(data,1), size(data,2), n_frames]);

for idx = 1:n_frames
    imwrite(uint16(syndata2(:,:,idx)), savepath2, 'tiff', 'WriteMode', 'append');
end