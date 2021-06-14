function run_rollBase_dFoF(varargin)
%% Function to generate a (uint8) dF/F timeseries with a rolling baseline window
% Input arguments: [winsize, mode, output]
% winsize: length of rolling baseline window in [s]
% mode: mode of normalization ('roll', 'grand', 'px')
% 'grand': normalize dF by average across all px and timepoints
% 'px': normalize dF by pixelwise average across all timepoints
% 'roll': normalize dF by rolling baseline average
% output: type of output ('dFoF', 'dF', 'FoF')

%% Add paths
wd = pwd(); p = split(wd,'\'); p = p(1:end-1); pp = []; for ii=1:numel(p), pp = fullfile(pp,p{ii}); end
addpath(genpath(pp)), clear('p','pp');

if isempty(varargin)
%     winsize = 0.5; % [s]
    winsize = 10; % [frames]
    mode = 'roll';
    output = 'dFoF';
else
    if isempty(varargin{1}), winsize = 10; %winsize = 0.5;
    else, winsize = varargin{1};  end
    if isempty(varargin{2}), mode = 'roll';
    else, mode = varargin{2};  end
    if isempty(varargin{3}), output = 'dFoF';
    else, output = varargin{3};  end
end

[filenames, pathname] = uigetfile({'*.tif*'; '*.lsm'; '*.nd2'}, 'Select files for dF/F calculation','Multiselect', 'on');
if isa(filenames,'char'), filenames = {filenames}; end
nfiles = size(filenames,2);

frametimes = NaN(nfiles,1);
for iF = 1:nfiles
    answer = input(sprintf('Frametime (ms) for %s: ',filenames{iF}));
    frametimes(iF) = answer/1000; % [s]
end

for iF = 1:nfiles
    fprintf('** Start with %s **\n', filenames{iF});
%     tmpwin = round(winsize/frametimes(iF)); % [s]
    tmpwin = winsize;
    tmpim = read_timeseries(fullfile(pathname,filenames{iF}), 'load');
    tmpim = double(tmpim);
    
    imsize = size(tmpim);
    n_px = imsize(1)*imsize(2);
    n_frames = imsize(3);
    lintraces = reshape(tmpim,[n_px n_frames]);
    
    disp(sprintf('Calculating %s values...', output));
    [dFtraces, FoFtraces, dFoFtraces] = rollBase_dFoF(lintraces, tmpwin, 5, mode);
    
    if strcmp(output,'dF')
        finalim = reshape(dFtraces,[imsize(1) imsize(2) n_frames]);
        tmpname = 'dF_';
    elseif strcmp(output,'dFoF')
        finalim = reshape(dFoFtraces,[imsize(1) imsize(2) n_frames]);
        tmpname = 'dFoF_';
    else
        finalim = reshape(FoFtraces,[imsize(1) imsize(2) n_frames]);
        tmpname = 'FoF_';
    end
    
    %% Output
    fprintf('Writing dF/F timeseries for %s ...\n', filenames{iF});
    cutfilename = split(filenames{iF},'.');
    cutfilename = cutfilename(1:end-1);
    if size(cutfilename,1) > 1
        cutfilename_new = [];
        for ii = 1:size(cutfilename,1), cutfilename_new = [cutfilename_new '.' cutfilename{ii}]; end
        cutfilename = cutfilename_new(2:end);
    end
    tmpp = strcat(pathname, tmpname, cutfilename, '.tif');
    if isa(tmpp,'cell'), tmpp = tmpp{1}; end
    
    % Scale
    scaled_finalim = finalim - min(finalim,[],'all');
    scaled_finalim = scaled_finalim./max(scaled_finalim,[],'all');
    scaled_finalim = scaled_finalim.*255;
    finalim = scaled_finalim;

    % Write output stack
    for idx = 1:n_frames
        imwrite(uint8(finalim(:,:,idx)), tmpp, 'tiff', 'WriteMode', 'append');
    end
    
    % Print unscaled value parameters and histogram
    fprintf('Writing unscaled dF/F parameters for %s ...\n', filenames{iF});
    tmppmean = strcat(tmpp(1:end-4), '_mean.xlsx');
    tmppsd = strcat(tmpp(1:end-4), '_sd.xlsx');
    if strcmp(output,'dF')
%         tbl.mean = mean(dFtraces, 'all');
%         tbl.sd = std(dFtraces,[], 'all');
    elseif strcmp(output,'dFoF')
%         tbl.mean = mean(dFoFtraces, 'all');
%         tbl.sd = std(dFoFtraces,[], 'all');
    else
        mean_mtrx = reshape(mean(FoFtraces, 2), [imsize(1) imsize(2)]);
        sd_mtrx = reshape(std(FoFtraces,[], 2), [imsize(1) imsize(2)]);
%         writematrix(mean_mtrx, tmppmean);
%         writematrix(sd_mtrx, tmppsd);
        
        % Save noise histogram of individual px
        px_b_y = [65 82 370 117]; px_b_x = [261 477 443 18]; 
        px_s_y = [217 169 29 421]; px_s_x = [380 446 270 279];
        lin_b_idx = NaN(numel(px_b_y), 1); lin_s_idx = NaN(numel(px_b_y), 1);
        
        for ipx = 1:numel(px_b_y)
            lin_b_idx(ipx) = sub2ind([imsize(1) imsize(2)], px_b_y(ipx), px_b_x(ipx));
            lin_s_idx(ipx) = sub2ind([imsize(1) imsize(2)], px_s_y(ipx), px_s_x(ipx));
        end
        edges = -50:0.1:50;
        for ipx = 1:numel(px_b_y)
            tmpoutname_b = strcat(tmpp(1:end-4),'_B_',num2str(px_b_y(ipx)), '_',num2str(px_b_x(ipx)), '.xlsx');
            tmpoutname_s = strcat(tmpp(1:end-4),'_S_',num2str(px_s_y(ipx)), '_',num2str(px_s_x(ipx)), '.xlsx');
            
            [cnts, ~] = histcounts(FoFtraces(lin_b_idx(ipx),:), edges);
            writematrix([edges(2:end)' cnts'], tmpoutname_b);
            
            [cnts, ~] = histcounts(FoFtraces(lin_s_idx(ipx),:), edges);
            writematrix([edges(2:end)' cnts'], tmpoutname_s);
        end
%         tmpphist = strcat(tmpp(1:end-4),'_hist.xlsx');
%         [cnts, ~] = histcounts(FoFtraces, edges);
%         writematrix([edges(2:end)' cnts'], tmpphist);
    end
    disp('Finished.');
end

end