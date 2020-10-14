function save_files(traceindices)

global FULLFILENAMES ROILIST evINFO FFORMAT SAVEPATH

numtr = numel(traceindices);

[~,n,~] = fileparts(FULLFILENAMES);
pn = fullfile(SAVEPATH,n);

if numtr == 1
    names = {'_binary_traces', '_event_info'};
    trname = strcat('_',ROILIST{traceindices});
    names = strcat(repelem({pn},numel(names)), repelem({trname}, numel(names)), names, repelem({FFORMAT},numel(names)));
    
    % Binary traces
    tbl = table();
    tbl.timepoints = FTIMEVEC;
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
    
else
    names = {'_crossings_binary',...
        '_peaks_binary',...
        'crossing_idx',...
        'peaks_idx',...
        '_amplitude',...
        '_amplitude_dFoF',...
        '_risetimes',...
        '_decaytimes',...
        '_iei',...
        '_baseline'};
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
        tmpelem = numel(evINFO(traceindices(iTr)).crossidx);
        
        % Binary crossings
        tbl_cross_binary.timepoints = FTIMEVEC;
        tbl_peaks_binary.timepoints = FTIMEVEC;
        tbl_cross_binary.(ROILIST{tmpidx}) = uint8(evINFO(tmpidx).binarycross);
        tbl_peaks_binary.(ROILIST{tmpidx}) = uint8(evINFO(tmpidx).binarypeaks);
        
        % Event info
        if tmpelem < maxelem
            fill = maxelem-tmpelem;
            tbl_crossidx.(ROILIST{tmpidx}) = [evINFO(tmpidx).crossidx; (repelem(0,fill))'];
            tbl_peakidx.(ROILIST{tmpidx}) = [evINFO(tmpidx).peakidx; (repelem(0,fill))'];
            tbl_amps.(ROILIST{tmpidx}) = [evINFO(tmpidx).amps; (repelem(0,fill))'];
            tbl_amps_dFoF.(ROILIST{tmpidx}) = [evINFO(tmpidx).amps_dFoF; (repelem(0,fill))'];
            tbl_rise.(ROILIST{tmpidx}) = [evINFO(tmpidx).risetimes; (repelem(0,fill))'];
            tbl_decay.(ROILIST{tmpidx}) = [evINFO(tmpidx).decaytimes; (repelem(0,fill))'];
            tbl_ieis.(ROILIST{tmpidx}) = [evINFO(tmpidx).ieis; (repelem(0,fill+1))'];
            tbl_baselines.(ROILIST{tmpidx}) = [evINFO(tmpidx).baselinevalues; (repelem(0,fill))'];
       else       
            tbl_crossidx.(ROILIST{tmpidx}) = evINFO(tmpidx).crossidx;
            tbl_peakidx.(ROILIST{tmpidx}) = evINFO(tmpidx).peakidx;
            tbl_amps.(ROILIST{tmpidx}) = evINFO(tmpidx).amps;
            tbl_amps_dFoF.(ROILIST{tmpidx}) = evINFO(tmpidx).amps_dFoF;
            tbl_rise.(ROILIST{tmpidx}) = evINFO(tmpidx).risetimes;
            tbl_decay.(ROILIST{tmpidx}) = evINFO(tmpidx).decaytimes;
            tbl_ieis.(ROILIST{tmpidx}) = [evINFO(tmpidx).ieis;NaN];
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
    
    % Save struct()
    strname = strcat(pn,'_event_struct.m');
    save(strname, 'evINFO', '-v7.3');
    disp('Writing finished!');
end

end