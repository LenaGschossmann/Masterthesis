%% ***** Create synthetic tetst timeseries *****

clear all;

% Add paths
wd = pwd(); p = split(wd,'\'); p = p(1:end-1); pp = []; for ii=1:numel(p), pp = fullfile(pp,p{ii}); end
addpath(genpath(pp)), clear('p','pp');

%% Load stack
[f_ori,p_ori] = uigetfile({'*.tif*'; '*.lsm'; '*.nd2'}, 'Select time series','Multiselect', 'off');
[f_FoF,p_FoF] = uigetfile({'*.xlsx*'}, 'Select FoF mean and sd files','Multiselect', 'on', p_ori);
if size(f_FoF,2) ~= 2
    error('Please select exactly 2 files!');
else
    if ~isempty(regexp(f_FoF{1},'mean')), fp_mean_FoF = fullfile(p_FoF,f_FoF{1}); fp_sd_FoF = fullfile(p_FoF,f_FoF{2});
    else, fp_mean_FoF = fullfile(p_FoF,f_FoF{2}); fp_sd_FoF = fullfile(p_FoF,f_FoF{1});
    end
end

fp_ori = fullfile(p_ori,f_ori);

data_ori = read_timeseries(fp_ori, 'load');
data_aip = mean(data_ori,3);
mean_FoF = readmatrix(fp_mean_FoF);
sd_FoF = readmatrix(fp_sd_FoF);

n_px = size(data_ori,1)*size(data_ori,2);
n_frames = size(data_ori,3);
n_elem = n_frames*n_px;

%% Output path
outname_ctrl1 = strcat(p_ori, 'CTRL1_', f_ori(1:end-4),'.tif'); % AIP with multiplicative gauss noise
outname_ctrl1_noise = strcat(p_ori, 'CTRL1_noise_', f_ori(1:end-4),'.xlsx'); % AIP with multiplicative gauss noise
outname_ctrl2 = strcat(p_ori, 'CTRL2_', f_ori(1:end-4),'.tif'); % Scrambled stack

%% Ctrl 1: AIP + multiplicative gauss noise
aipseries = repmat(data_aip,[1 1 n_frames]);
% mu = data_FoF.mean;
% % mu = 1;
% sigma = data_FoF.sd;

noise = mean_FoF + randn(size(data_aip,1),size(data_aip,2),n_frames) .* sd_FoF;
outdata_ctrl1 = aipseries.*noise;

% % Write ctrl series and noise
% for idx = 1:n_frames
%     imwrite(uint16(outdata_ctrl1(:,:,idx)), outname_ctrl1, 'tiff', 'WriteMode', 'append');
% end

% Scale to uint8
% scaled_noise = noise - min(noise,[],'all');
% scaled_noise = scaled_noise./max(scaled_noise,[],'all');
% scaled_noise = scaled_noise.*255; 
% for idx = 1:n_frames
%     imwrite(uint8(scaled_noise(:,:,idx)), outname_ctrl1_noise, 'tiff', 'WriteMode', 'append');
% end

% Histogram of noise
% edges = -50:0.1:50;
% [cnts, ~] = histcounts(noise, edges);
% writematrix([edges(2:end)' cnts'], outname_ctrl1_noise);

% Save noise histogram of individual px
px_b_y = [65 82 370 117]; px_b_x = [261 477 443 18]; 
px_s_y = [217 169 29 421]; px_s_x = [380 446 270 279];

% Histogram of px-wise noise
edges = -50:0.1:50;

for ipx = 1:numel(px_b_y)
    tmpoutname_b = strcat(p_ori, 'CTRL1_noise_', f_ori(1:end-4),'_B_',num2str(px_b_y(ipx)), '_',num2str(px_b_x(ipx)), '.xlsx');
    tmpoutname_s = strcat(p_ori, 'CTRL1_noise_', f_ori(1:end-4),'_S_',num2str(px_s_y(ipx)), '_',num2str(px_s_x(ipx)), '.xlsx');
    
    [cnts, ~] = histcounts(noise(px_b_y(ipx), px_b_x(ipx),:), edges);
    writematrix([edges(2:end)' cnts'], tmpoutname_b);
    
    [cnts, ~] = histcounts(noise(px_s_y(ipx), px_s_x(ipx),:), edges);
    writematrix([edges(2:end)' cnts'], tmpoutname_s);    
end

% % Ctrl 2: Scrambled stack
% rnd_data = reshape(data_ori,[n_elem 1]);
% rnd_data = rnd_data(randperm(n_elem));
% outdata_ctrl2 = reshape(rnd_data,[size(data_ori,1), size(data_ori,2), n_frames]);
% 
% for idx = 1:n_frames
%     imwrite(uint16(outdata_ctrl2(:,:,idx)), outname_ctrl2, 'tiff', 'WriteMode', 'append');
% end

% %% **** Create synthetic control trace ****
% % 2PM
% % A896_ts_2: roi6, roi10, roi11, roi12, roi17
% % mu = [2.868 2.441 3.312 3.407 4.139];
% % sigma = [0.322 0.318 0.311 0.797 2.323];
% mu = [0.0210159]; % dF/F: ts_2 roi7
% sigma = [0.684875]; %dF/F: ts2 roi7
% nf = 3000;
% p='C:\Users\lena_\Projects\iGluSnFR\Analysis\Tests\LenaMatlabCode\DD_ROI_detection\ev_detection\test_thresholds\Ctrl_traces\2pm_A896_ts7.xlsx';
% 
% % %Conf
% % A906_sl_4_ts_1: roi3, roi4, roi5, roi6, roi1_8
% % mu = [0.424 0.425 0.431 0.443 0.546];
% % sigma = [0.048 0.019 0.07 0.475 0.222];
% % nf = 1200;
% % p='C:\Users\lena_\Projects\iGluSnFR\Analysis\Tests\LenaMatlabCode\DD_ROI_detection\ev_detection\test_thresholds\Ctrl_traces\conf_A906.xlsx';
% 
% rndtrc = NaN(numel(mu), nf);
% for ii =1:numel(mu), rndtrc(ii,:) = mu(ii)+randn(1,nf).*sigma(ii); end
% 
% writematrix(rndtrc, p);
