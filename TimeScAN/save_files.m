function save_files(traceindices)

global FULLFILENAMES ROILIST evINFO FFORMAT SAVEPATH FTIMEVEC PREPOINTS...
    POSTPOINTS TOTP TRACEDATA SMTHWIN

warning('off','MATLAB:xlswrite:AddSheet');

numtr = numel(traceindices);

[~,n,~] = fileparts(FULLFILENAMES);
pn = fullfile(SAVEPATH,n);

croppedtraces = [];

if numtr == 1
    names = {'_binary_traces', '_event_info', '_event_summary', '_cropped_traces'};
    trname = strcat('_',ROILIST{traceindices});
    names = strcat(repelem({pn},numel(names)), repelem({trname}, numel(names)), names, repelem({FFORMAT},numel(names)));
    
    % Binary traces
    tbl = table();
    tbl.timepoints = FTIMEVEC';
    tbl.onset = uint8(evINFO(traceindices).binaryonset);
    tbl.peaks = uint8(evINFO(traceindices).binarypeaks);
    writetable(tbl, names{1});
    
    % Individual Event info
    tbl = table();
    tbl.onset_idx = evINFO(traceindices).onsetidx;
    tbl.peaks_idx = evINFO(traceindices).peakidx;
    tbl.amplitude = evINFO(traceindices).amps;
    tbl.amplitude_dFoF = evINFO(traceindices).amps_dFoF;
    tbl.risetime = evINFO(traceindices).risetimes;
    tbl.decaytime = evINFO(traceindices).decaytimes;
    tbl.iei = [evINFO(traceindices).ieis; NaN];
    tbl.baseline_values = evINFO(traceindices).baselinevalues;
    %!!tbl.roi_area_um = repmat(AREADATA(traceindices),[numel(tbl.onset_idx) 1]);
    writetable(tbl, names{2});
    
    % Summarized Event info
    tbl = table();
    tbl.av_amplitude = evINFO(traceindices).avamp;
    tbl.av_amplitude_dFoF = evINFO(traceindices).avamp_dFoF;
    tbl.trace_range = strcat(num2str(evINFO(traceindices).rangein),' - ', num2str(evINFO(traceindices).rangeout));
    tbl.trace_av = evINFO(traceindices).average; % Within given range
    tbl.trace_sd = evINFO(traceindices).sd; % Within given range
    tbl.av_iei = evINFO(traceindices).aviei;
    tbl.cv_iei = evINFO(traceindices).cviei;
    tbl.eventrate = evINFO(traceindices).eventrate;
    %!!tbl.roi_area_um = AREADATA(traceindices);
    writetable(tbl,names{3});
    
    % Cropped traces
    [croppedtraces, traceids] = crop_trace(traceindices);
    writematrix(croppedtraces, names{4}, 'Sheet', 'Traces');
    writecell(traceids, names{4}, 'Sheet', 'Trace IDs');
    
    % Save struct()
    strname = strcat(pn,'_event_struct.m');
    save(strname, 'evINFO', '-v7.3');
    disp('Writing finished!');
    
else
    names = {'onsets',...
        'peaks',...
        'onset_idx',...
        'peaks_idx',...
        'amp',...
        'amp_dFoF',...
        'risetimes',...
        'decaytimes',...
        'iei',...
        'baseline',...
        'roi_area_um'};

    binaryname = strcat(pn, '_binary_events', '.xlsx');
    eventinfoname = strcat(pn, '_event_info', '.xlsx');
    eventsummaryname = strcat(pn, '_event_summary', '.xlsx');
    croptrname = strcat(pn, '_cropped_traces', '.xlsx');
    
    tbl_onset_binary = struct();
    tbl_peaks_binary = struct();
    tbl_onsetidx = struct();
    tbl_peakidx = struct();
    tbl_amps = struct();
    tbl_amps_dFoF = struct();
    tbl_rise = struct();
    tbl_decay = struct();
    tbl_ieis = struct();
    tbl_baselines = struct();
    tbl_areas = struct();
    maxelem = 0;
    
    for iTr = 1:numtr, maxelem = max([maxelem numel(evINFO(traceindices(iTr)).onsetidx)]); end
    
    for iTr = 1:numtr
        tmpidx = traceindices(iTr);
        tmpelem = numel(evINFO(tmpidx).onsetidx);
        
        % Binary onsets
        tbl_onset_binary.timepoints = FTIMEVEC';
        tbl_peaks_binary.timepoints = FTIMEVEC';
        if all(isnan(evINFO(tmpidx).binaryonset)) || isempty(evINFO(tmpidx).binaryonset)
            tbl_onset_binary.(ROILIST{tmpidx}) = zeros(size(tbl_onset_binary.timepoints));
            tbl_peaks_binary.(ROILIST{tmpidx}) = zeros(size(tbl_onset_binary.timepoints));
        else
            tbl_onset_binary.(ROILIST{tmpidx}) = uint8(evINFO(tmpidx).binaryonset);
            tbl_peaks_binary.(ROILIST{tmpidx}) = uint8(evINFO(tmpidx).binarypeaks);
            % Crop traces
            if isempty(croppedtraces)
                [croppedtraces, traceids] = crop_trace(tmpidx);
            else
                [tmptrace, tmpids] = crop_trace(tmpidx);
                croppedtraces = [croppedtraces tmptrace]; traceids = [traceids, tmpids];
            end
        end
        
        % Event info
        if tmpelem < maxelem
            fill = maxelem-tmpelem;
            tbl_onsetidx.(ROILIST{tmpidx}) = [evINFO(tmpidx).onsetidx; (repelem(NaN,fill))'];
            tbl_peakidx.(ROILIST{tmpidx}) = [evINFO(tmpidx).peakidx; (repelem(NaN,fill))'];
            tbl_amps.(ROILIST{tmpidx}) = [evINFO(tmpidx).amps; (repelem(NaN,fill))'];
            tbl_amps_dFoF.(ROILIST{tmpidx}) = [evINFO(tmpidx).amps_dFoF; (repelem(NaN,fill))'];
            tbl_rise.(ROILIST{tmpidx}) = [evINFO(tmpidx).risetimes; (repelem(NaN,fill))'];
            tbl_decay.(ROILIST{tmpidx}) = [evINFO(tmpidx).decaytimes; (repelem(NaN,fill))'];
            tbl_ieis.(ROILIST{tmpidx}) = [evINFO(tmpidx).ieis; (repelem(NaN,maxelem-1-numel(evINFO(tmpidx).ieis)))'];
            tbl_baselines.(ROILIST{tmpidx}) = [evINFO(tmpidx).baselinevalues; (repelem(NaN,fill))'];
        else       
            tbl_onsetidx.(ROILIST{tmpidx}) = evINFO(tmpidx).onsetidx;
            tbl_peakidx.(ROILIST{tmpidx}) = evINFO(tmpidx).peakidx;
            tbl_amps.(ROILIST{tmpidx}) = evINFO(tmpidx).amps;
            tbl_amps_dFoF.(ROILIST{tmpidx}) = evINFO(tmpidx).amps_dFoF;
            tbl_rise.(ROILIST{tmpidx}) = evINFO(tmpidx).risetimes;
            tbl_decay.(ROILIST{tmpidx}) = evINFO(tmpidx).decaytimes;
            tbl_ieis.(ROILIST{tmpidx}) = [evINFO(tmpidx).ieis];
            tbl_baselines.(ROILIST{tmpidx}) = evINFO(tmpidx).baselinevalues;
        end
        %!!tbl_areas.(ROILIST{tmpidx}) = AREADATA(tmpidx);
        
        % Summarized Event info
        tbl = struct();
        tbl.av_amplitude = evINFO(tmpidx).avamp;
        tbl.av_amplitude_dFoF = evINFO(tmpidx).avamp_dFoF;
        tbl.trace_range = strcat(num2str(evINFO(tmpidx).rangein),' - ', num2str(evINFO(tmpidx).rangeout));
        tbl.trace_av = evINFO(tmpidx).average; % Within given range
        tbl.trace_sd = evINFO(tmpidx).sd; % Within given range
        tbl.av_iei = evINFO(tmpidx).aviei;
        tbl.cv_iei = evINFO(tmpidx).cviei;
        tbl.eventrate = evINFO(tmpidx).eventrate;
        %!!tbl.roi_area_um = AREADATA(tmpidx);
        writetable(struct2table(tbl), eventsummaryname, 'Sheet', ROILIST{tmpidx});
    end

    disp('Writing files...');    
    writetable(struct2table(tbl_onset_binary), binaryname, 'Sheet', names{1});
    writetable(struct2table(tbl_peaks_binary), binaryname, 'Sheet', names{2});
    writetable(struct2table(tbl_onsetidx), eventinfoname, 'Sheet', names{3});
    writetable(struct2table(tbl_peakidx), eventinfoname, 'Sheet', names{4});
    writetable(struct2table(tbl_amps), eventinfoname, 'Sheet', names{5});
    writetable(struct2table(tbl_amps_dFoF), eventinfoname, 'Sheet', names{6});
    writetable(struct2table(tbl_rise), eventinfoname, 'Sheet', names{7});
    writetable(struct2table(tbl_decay), eventinfoname, 'Sheet', names{8});
    writetable(struct2table(tbl_ieis), eventinfoname, 'Sheet', names{9});
    writetable(struct2table(tbl_baselines), eventinfoname, 'Sheet', names{10});
    writetable(struct2table(tbl_areas), eventinfoname, 'Sheet', names{11});
    writematrix(croppedtraces, croptrname, 'Sheet', 'Traces');
    writecell(traceids, croptrname, 'Sheet', 'Trace IDs'); % Info: Roi ID (row 1), Event ID (row 2), Onset index (row 3), Frametime (row 4)
    
    % Save struct()
    strname = strcat(pn,'_event_struct.m');
    save(strname, 'evINFO', '-v7.3');
    disp('Writing finished!');
    
    if numtr == numel(ROILIST)
        close all;
        clear all;
    end
end

roidistmtrx = ML_position(traceindices);
roinames = ROILIST(traceindices);
tbl = table();
tbl.Roi = roinames;
tbl.Distance_to_somata = roidistmtrx;
writetable(tbl, strcat(pn,'Roi_to_soma_dist.xlsx'));

    function [tmpcroptrace, tmptraceids] = crop_trace(tmpidx)
        %% Cropped traces
        tmpelem = numel(evINFO(tmpidx).onsetidx);
        tmpcroptrace = zeros(tmpelem,PREPOINTS+POSTPOINTS+1);
        tmptraceids = cell(4,tmpelem);
        worktrace = TRACEDATA(:,tmpidx);
        if SMTHWIN ~= 0, worktrace = smooth_data(worktrace,SMTHWIN); end
        
        for iE = 1:tmpelem
            onsetidx = evINFO(tmpidx).onsetidx(iE);
            baseline = evINFO(tmpidx).baselinevalues(iE);
            % Event preceding part
            if onsetidx-PREPOINTS < 1
                tmpcroptrace(iE,PREPOINTS-onsetidx+2:PREPOINTS+1) = worktrace(1:onsetidx)';
            else
                tmpcroptrace(iE,1:PREPOINTS+1) = worktrace(onsetidx-PREPOINTS:onsetidx)';
            end
            % Event following part
            if onsetidx+POSTPOINTS > TOTP
                tmpcroptrace(iE,PREPOINTS+2:PREPOINTS+1+numel(worktrace)-onsetidx) = worktrace(onsetidx+1:end);
            else
                tmpcroptrace(iE,PREPOINTS+2:end) = worktrace(onsetidx+1:onsetidx+POSTPOINTS);
            end
            % Delta F over F
            tmpcroptrace(iE,:) = (tmpcroptrace(iE,:)-baseline)./baseline;
            
            tmptraceids{1,iE} = ROILIST{tmpidx};
            tmptraceids{2,iE} = iE;
            tmptraceids{3,iE} = evINFO(tmpidx).onsetidx(iE);
            tmptraceids{4,iE} = evINFO(tmpidx).ftime;
        end
        tmpcroptrace = tmpcroptrace'; % columns: Events, rows: Timepoints
    end

end