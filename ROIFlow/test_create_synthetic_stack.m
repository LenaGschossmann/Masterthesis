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

%% Create synthetic control trace
% % 2PM
% % A896_ts_2: roi6, roi10, roi11, roi12, roi17
% mu = [2.868 2.441 3.312 3.407 4.139];
% sigma = [0.322 0.318 0.311 0.797 2.323];
% nf = 3000;
% p='C:\Users\lena_\Projects\iGluSnFR\Analysis\Tests\LenaMatlabCode\DD_ROI_detection\ev_detection\test_thresholds\Ctrl_traces\2pm_A896.xlsx';

% Conf
% A906_sl_4_ts_1: roi3, roi4, roi5, roi6, roi1_8
mu = [0.424 0.425 0.431 0.443 0.546];
sigma = [0.048 0.019 0.07 0.475 0.222];
nf = 1200;
p='C:\Users\lena_\Projects\iGluSnFR\Analysis\Tests\LenaMatlabCode\DD_ROI_detection\ev_detection\test_thresholds\Ctrl_traces\conf_A906.xlsx';

rndtrc = NaN(numel(mu), nf);
for ii =1:numel(mu), rndtrc(ii,:) = mu(ii)+randn(1,nf).*sigma(ii); end

writematrix(rndtrc, p);
