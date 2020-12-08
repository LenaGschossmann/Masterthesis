addpath(genpath('/media/2Photon/Matlab code/Martin Source Code/ca_source_extraction'));
name = '/media/2Photon/imaging Negar & Martin/M103/07.12.16/M103.071216.1325.tif';
options.crop = [20 10 10 10];
[Data,Datared] = readdata(name,options);
%%
amplitudeThreshold = [11.5 12.3];
win = [32 15];
% RippleNoiseRemoval(Data,amplitudeThreshold(1),win(1),1,[]);
% RippleNoiseRemoval(Datared,amplitudeThreshold(2),win(2),1,[]);
%%
[Data,~] = RippleNoiseRemoval(Data,amplitudeThreshold(1),win(1),0,[]);
% [Datared,~] = RippleNoiseRemoval(Datared,amplitudeThreshold(2),win(2),0,[]);

%%
Y = Data;
clear Data Datared

%%

TiffData(Y,[], '/media/2Photon/imaging Negar & Martin/M103/07.12.16','M103.071216.1325.tif','M103.071216.1325_Y.tif');