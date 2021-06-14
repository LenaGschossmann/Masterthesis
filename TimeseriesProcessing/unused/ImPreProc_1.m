%% Run Motion correction - NoRMCorre
% Pnevmatikakis, E. A., & Giovannucci, A. (2017). NoRMCorre: An online algorithm for piecewise rigid motion correction of calcium imaging data. Journal of neuroscience methods, 291, 83-94.

%% Batch selection
% Select .nd2 file from NikonM
[FileName, PathName] = uigetfile('*.nd2*', [],'G:\Lena\Masterthesis\In vivo\tmp', 'Multiselect', 'on');
% [FileName, PathName] = uigetfile('*.tif*', [],'G:\Katja\iglusnfrKatja', 'Multiselect', 'on');
if ~isa(FileName, 'cell'), FileName = {FileName}; end
files = FileName;

run_mocorr(PathName, files,false);