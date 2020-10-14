function [] = open_file(filepointer)

global IMHW COMPOSITE2D FTIME SPATH FNAME IMMETA PLOTRANGE

[SPATH, FNAME, ~] = fileparts(filepointer);
[tmpstack, imInfo] = load_bf_file(filepointer, true);

% Frametime from metadata
IMMETA =imInfo.metadata;
timestrings = flipud(IMMETA(strncmp(IMMETA, 'TimeStamp',9)));
timesplit = cellfun(@(x) regexp(x,'= ','split'), timestrings,'UniformOutput', false);
sortsplit = cellfun(@(x) x{1}, timesplit, 'UniformOutput', false);
timesplit = cellfun(@(x) str2double(x{2}), timesplit, 'UniformOutput', false);
sortsplit = cellfun(@(x) regexp(x,'#','split'), sortsplit,'UniformOutput', false);
sortsplit = cellfun(@(x) str2double(x{2}), sortsplit, 'UniformOutput', false);
sortsplit = cell2mat(sortsplit);
[~,sortidx] = sort(sortsplit, 'ascend');
timeseries = cell2mat(timesplit); timeseries = timeseries(sortidx);

end