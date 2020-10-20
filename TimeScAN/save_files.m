function save_files(traceindices)

global FULLFILENAMES ROILIST evINFO FFORMAT SAVEPATH FTIMEVEC PREPOINTS POSTPOINTS TOTP TRACEDATA SMTHWIN

numtr = numel(traceindices);

[~,n,~] = fileparts(FULLFILENAMES);
pn = fullfile(SAVEPATH,n);

croppedtraces = [];

if numtr == 1
    names = {'_binary_traces', '_event_info', '_cropped_traces'};
    trname = strcat('_',ROILIST{traceindices});
    names = strcat(repelem({pn},numel(names)), repelem({trname}, numel(names)), names, repelem({FFORMAT},numel(names)));
    
    % Binary traces
    tbl = table();
    tbl.timepoints = FTIMEVEC';
    tbl.crossing = uint8(evINFO(traceindices).binarycross);
    tbl.peaks = uint8(evINFO(traceindices).binarypeaks);
    writetable(tbl, names{1});
    
    % Event info
    tbl = table();
    tbl.crossing_idx = evINFO(traceindices).crossidx;
    tbl.peaks_idx = evINFO(traceindices).peakidx;
    tbl.amplitude = evINFO(traceindices).amps;
    tbl.amplitude_dFoF = evINFO(traceindices).amps_dFoF;
    tbl.risetime = evINFO(traceindices).risetimes;
    tbl.decaytime = evINFO(traceindices).decaytimes;
    tbl.iei = [evINFO(traceindices).ieis; NaN];
    tbl.baseline_values = evINFO(traceindices).baselinevalues;
    writetable(tbl, names{2});
    
    % Cropped traces
    croppedtraces = crop_trace(traceindices);
    writematrix(croppedtraces, names{3});
    
    
else
    names = {'_crossings_binary',...
        '_peaks_binary',...
        '_crossing_idx',...
        '_peaks_idx',...
        '_amplitude',...
        '_amplitude_dFoF',...
        '_risetimes',...
        '_decaytimes',...
        '_iei',...
        '_baseline',...
        '_cropped_traces'};
    names = strcat(repelem({pn},numel(names)), names, repelem({FFORMAT},numel(names)));
    
    tbl_cross_binary = struct();
    tbl_peaks_binary = struct();
    tbl_crossidx = struct();
    tbl_peakidx = struct();
    tbl_amps = struct();
    tbl_amps_dFoF = struct();
    tbl_rise = struct();
    tbl_decay = struct();
    tbl_ieis = struct();
    tbl_baselines = struct();
    maxelem = 0;
    
    for iTr = 1:numtr, maxelem = max([maxelem numel(evINFO(traceindices(iTr)).crossidx)]); end
    
    for iTr = 1:numtr
        tmpidx = traceindices(iTr);
        tmpelem = numel(evINFO(tmpidx).crossidx);
        
        % Binary crossings
        tbl_cross_binary.timepoints = FTIMEVEC';
        tbl_peaks_binary.timepoints = FTIMEVEC';
        if isnan(evINFO(tmpidx).binarycross)
            tbl_cross_binary.(ROILIST{tmpidx}) = zeros(size(tbl_cross_binary.timepoints));
            tbl_peaks_binary.(ROILIST{tmpidx}) = zeros(size(tbl_cross_binary.timepoints));
        else
            tbl_cross_binary.(ROILIST{tmpidx}) = uint8(evINFO(tmpidx).binarycross);
            tbl_peaks_binary.(ROILIST{tmpidx}) = uint8(evINFO(tmpidx).binarypeaks);
            % Crop traces
            if isempty(croppedtraces), croppedtraces = crop_trace(tmpidx);
            else, croppedtraces = [croppedtraces crop_trace(tmpidx)];
            end
        end
        
        % Event info
        if tmpelem < maxelem
            fill = maxelem-tmpelem;
            tbl_crossidx.(ROILIST{tmpidx}) = [evINFO(tmpidx).crossidx; (repelem(0,fill))'];
            tbl_peakidx.(ROILIST{tmpidx}) = [evINFO(tmpidx).peakidx; (repelem(0,fill))'];
            tbl_amps.(ROILIST{tmpidx}) = [evINFO(tmpidx).amps; (repelem(0,fill))'];
            tbl_amps_dFoF.(ROILIST{tmpidx}) = [evINFO(tmpidx).amps_dFoF; (repelem(0,fill))'];
            tbl_rise.(ROILIST{tmpidx}) = [evINFO(tmpidx).risetimes; (repelem(0,fill))'];
            tbl_decay.(ROILIST{tmpidx}) = [evINFO(tmpidx).decaytimes; (repelem(0,fill))'];
            tbl_ieis.(ROILIST{tmpidx}) = [evINFO(tmpidx).ieis; (repelem(0,maxelem-1-numel(evINFO(tmpidx).ieis)))'];
            tbl_baselines.(ROILIST{tmpidx}) = [evINFO(tmpidx).baselinevalues; (repelem(0,fill))'];
       else       
            tbl_crossidx.(ROILIST{tmpidx}) = evINFO(tmpidx).crossidx;
            tbl_peakidx.(ROILIST{tmpidx}) = evINFO(tmpidx).peakidx;
            tbl_amps.(ROILIST{tmpidx}) = evINFO(tmpidx).amps;
            tbl_amps_dFoF.(ROILIST{tmpidx}) = evINFO(tmpidx).amps_dFoF;
            tbl_rise.(ROILIST{tmpidx}) = evINFO(tmpidx).risetimes;
            tbl_decay.(ROILIST{tmpidx}) = evINFO(tmpidx).decaytimes;
            tbl_ieis.(ROILIST{tmpidx}) = [evINFO(tmpidx).ieis];
            tbl_baselines.(ROILIST{tmpidx}) = evINFO(tmpidx).baselinevalues;
        end
      
    end

    disp('Writing files...');
    writetable(struct2table(tbl_cross_binary), names{1});
    writetable(struct2table(tbl_peaks_binary), names{2});
    writetable(struct2table(tbl_crossidx), names{3});
    writetable(struct2table(tbl_peakidx), names{4});
    writetable(struct2table(tbl_amps), names{5});
    writetable(struct2table(tbl_amps_dFoF), names{6});
    writetable(struct2table(tbl_rise), names{7});
    writetable(struct2table(tbl_decay), names{8});
    writetable(struct2table(tbl_ieis), names{9});
    writetable(struct2table(tbl_baselines), names{10});
    writematrix(croppedtraces, names{11});
    
    % Save struct()
    strname = strcat(pn,'_event_struct.m');
    save(strname, 'evINFO', '-v7.3');
    disp('Writing finished!');
    
    if numtr == numel(ROILIST)
        close all;
        clear all;
    end
end

    function [tmpcroptrace] = crop_trace(tmpidx)
        %% Cropped traces
        tmpelem = numel(evINFO(tmpidx).crossidx);
        tmpcroptrace = NaN(tmpelem,PREPOINTS+POSTPOINTS+1);
        worktrace = TRACEDATA(:,tmpidx);
        if SMTHWIN ~= 0, worktrace = smooth_data(worktrace,SMTHWIN); end
        
        for iE = 1:tmpelem
            crossidx = evINFO(tmpidx).crossidx(iE);
            baseline = evINFO(tmpidx).baselinevalues(iE);
            % Event preceding part
            if crossidx-PREPOINTS < 1
                tmpcroptrace(iE,PREPOINTS-crossidx+2:PREPOINTS+1) = worktrace(1:crossidx)';
            else
                tmpcroptrace(iE,1:PREPOINTS+1) = worktrace(crossidx-PREPOINTS:crossidx)';
            end
            % Event following part
            if crossidx+POSTPOINTS > TOTP
                tmpcroptrace(iE,PREPOINTS+2:PREPOINTS+1+numel(worktrace)-crossidx) = worktrace(crossidx+1:end);
            else
                tmpcroptrace(iE,PREPOINTS+2:end) = worktrace(crossidx+1:crossidx+POSTPOINTS);
            end
            % Delta F over F
            tmpcroptrace(iE,:) = (tmpcroptrace(iE,:)-baseline)./baseline;
        end
        tmpcroptrace = tmpcroptrace'; % columns: Events, rows: Timepoints
    end

end