function [] = open_file(filepointer, rangein1, rangein2)

global IMHW COMPOSITE2D FTIME SPATH FNAME IMMETA PLOTRANGE

[SPATH, FNAME, ~] = fileparts(filepointer);
[tmpstack, imInfo] = load_bf_file(filepointer, true);

% Collapse 3rd dimension
dims = size(tmpstack);
COMPOSITE2D = zeros(dims(1)*dims(3), dims(2));
catrange = 1:dims(1);
for frame = 1:dims(3)
    COMPOSITE2D(catrange,:) = tmpstack(:,:,frame);
    catrange = catrange+dims(1);
end
IMHW = size(COMPOSITE2D);

% Frametime from metadata
IMMETA =imInfo.metadata;
timestrings = flipud(IMMETA(strncmp(IMMETA, 'timestamp',9)));
ftimeseries = cellfun(@(x) regexp(x,'= ','split'), timestrings,'UniformOutput', false);
ftimeseries = cellfun(@(x) str2double(x{2}), ftimeseries, 'UniformOutput', false);
FTIME = mean(diff(cell2mat(ftimeseries)))/dims(1); %[s]

% Update name in control box figure
figures = get(groot,'Children');
tmpfig = figures(strncmp({figures(:).Name}, 'Line',4));
namebox = findobj(tmpfig, 'type', 'uicontrol', 'style', 'text');
set(namebox(1), 'string', FNAME);

PLOTRANGE = [1 size(COMPOSITE2D,1)];
set(rangein1,'String', num2str(PLOTRANGE(1)));
set(rangein2,'String', num2str(PLOTRANGE(2)));

end