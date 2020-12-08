function save_files_2(traceindices)

global FULLFILENAMES ROILIST evINFO SAVEPATH FTIMEVEC PREPOINTS...
    POSTPOINTS TOTP TRACEDATA SMTHWIN FTIME LOCAREAINFO SOMAINFO PXCLASSCNT

warning('off','MATLAB:xlswrite:AddSheet');

numtr = numel(traceindices);

[~,n,~] = fileparts(FULLFILENAMES);
pn = fullfile(SAVEPATH,n);

% Check if files with results already exist
tmpfiles = dir(SAVEPATH);
if size(tmpfiles,1) > 2
    prevfiles = [];
    for iF = 3:size(tmpfiles), prevfiles = [prevfiles; {tmpfiles(iF).name}]; end
    
    if any(~cellfun(@isempty ,regexp(prevfiles, '_binary_events', 'match')))
        tmpid = find(~cellfun(@isempty ,regexp(prevfiles, '_binary_events', 'match')) == 1);
        v = []; for iN = 1:numel(tmpid), p = split(prevfiles{tmpid(iN)}(1:end-5),'_'); v = [v str2num(p{end})]; end
        if isempty(v), v=0; end
        binaryname = strcat(pn, '_binary_events_', num2str(max(v)+1),'.xlsx'); % Table with binary traces
    else, binaryname = strcat(pn, '_binary_events', '.xlsx'); % Table with binary traces
    end
    if any(~cellfun(@isempty ,regexp(prevfiles, '_cropped_traces', 'match')))
        tmpid = find(~cellfun(@isempty ,regexp(prevfiles, '_cropped_traces', 'match')) == 1);
        v = []; for iN = 1:numel(tmpid), p = split(prevfiles{tmpid(iN)}(1:end-5),'_'); v = [v str2num(p{end})]; end
        if isempty(v), v=0; end
        croppedtrname = strcat(pn, '_cropped_traces_', num2str(max(v)+1),'.xlsx'); % Table with binary traces
    else, croppedtrname = strcat(pn, '_cropped_traces', '.xlsx'); % Table with cropped events
    end
    if any(~cellfun(@isempty ,regexp(prevfiles, '_event_info', 'match')))
        tmpid = find(~cellfun(@isempty ,regexp(prevfiles, '_event_info', 'match')) == 1);
        v = []; for iN = 1:numel(tmpid), p = split(prevfiles{tmpid(iN)}(1:end-5),'_'); v = [v str2num(p{end})]; end
        if isempty(v), v=0; end
        eventinfoname1 = strcat(pn, '_event_info_', num2str(max(v)+1),'.xlsx'); % Table with binary traces
    else, eventinfoname1 = strcat(pn, '_event_info', '.xlsx'); % Table with all events
    end
    if any(~cellfun(@isempty ,regexp(prevfiles, '_iei_info', 'match')))
        tmpid = find(~cellfun(@isempty ,regexp(prevfiles, '_iei_info', 'match')) == 1);
        v = []; for iN = 1:numel(tmpid), p = split(prevfiles{tmpid(iN)}(1:end-5),'_'); v = [v str2num(p{end})]; end
        if isempty(v), v=0; end
        eventinfoname2 = strcat(pn, '_iei_info_', num2str(max(v)+1),'.xlsx'); % Table with binary traces
    else, eventinfoname2 = strcat(pn, '_iei_info', '.xlsx'); % Table with all events' ieis
    end
    if any(~cellfun(@isempty ,regexp(prevfiles, '_event_summary', 'match')))
        tmpid = find(~cellfun(@isempty ,regexp(prevfiles, '_event_summary', 'match')) == 1);
        v = []; for iN = 1:numel(tmpid), p = split(prevfiles{tmpid(iN)}(1:end-5),'_'); v = [v str2num(p{end})]; end
        if isempty(v), v=0; end
        eventsummaryname = strcat(pn, '_event_summary_', num2str(max(v)+1),'.xlsx'); % Table with binary traces
    else, eventsummaryname = strcat(pn, '_event_summary', '.xlsx'); % Table with event summary statistics
    end
    if any(~cellfun(@isempty ,regexp(prevfiles, '_event_struct', 'match')))
        tmpid = find(~cellfun(@isempty ,regexp(prevfiles, '_event_struct', 'match')) == 1);
        v = []; for iN = 1:numel(tmpid), p = split(prevfiles{tmpid(iN)}(1:end-3),'_'); v = [v str2num(p{end})]; end
        if isempty(v), v=0; end
        strname = strcat(pn,'_event_struct_', num2str(max(v)+1),'.xlsx'); % Table with binary traces
    else, strname = strcat(pn,'_event_struct.m'); % Save struct as .mat
    end
else
    binaryname = strcat(pn, '_binary_events', '.xlsx'); % Table with binary traces
    croppedtrname = strcat(pn, '_cropped_traces', '.xlsx'); % Table with cropped events
    eventinfoname1 = strcat(pn, '_event_info', '.xlsx'); % Table with all events
    eventinfoname2 = strcat(pn, '_iei_info', '.xlsx'); % Table with all events' ieis
    eventsummaryname = strcat(pn, '_event_summary', '.xlsx'); % Table with event summary statistics
    strname = strcat(pn,'_event_struct.m'); % Save struct as .mat
end

% Initialize tables
strct_onset_binary = struct();
strct_peaks_binary = struct();
tbl_events1 = table();
tbl_events2 = table();
tbl_summary = table();
cell_cropped_trc = cell(0);
cell_cropped_trc{1,1} = 'ROI#'; cell_cropped_trc{2,1} = 'Event#';
cell_cropped_trc{3,1} = 'Onset Index'; cell_cropped_trc{4,1} = 'Frametime';
icroptr = 1;
evcnter = 1;
ieicnter = 1;
sumcnter = 1;

if ~isempty(SOMAINFO)
    disp('Calculating distance information...');
    roidistmtrx = ML_position(traceindices);
end

disp('Preparing files for saving...');
for iTr = 1:numtr
    tmpidx = traceindices(iTr);
    tmpelem = numel(evINFO(tmpidx).onsetidx);
    
    if ~isempty(LOCAREAINFO), tmparea = LOCAREAINFO(iTr,1); else, tmparea = NaN; end
    if~isempty(SOMAINFO), tmpdist = roidistmtrx(iTr); else, tmpdist = NaN; end
    
    % Binary information
    strct_onset_binary.timepoints = FTIMEVEC';
    strct_peaks_binary.timepoints = FTIMEVEC';
    if all(isnan(evINFO(tmpidx).binaryonset)) || isempty(evINFO(tmpidx).binaryonset)
        strct_onset_binary.(ROILIST{tmpidx}) = zeros(size(strct_onset_binary.timepoints));
        strct_peaks_binary.(ROILIST{tmpidx}) = zeros(size(strct_onset_binary.timepoints));
        tbl_summary.Roi(sumcnter) = ROILIST(tmpidx);
        tbl_summary.Tot_time_s(sumcnter) = TOTP*FTIME;
        tbl_summary.Tot_dend_area_sqr_um(sumcnter) = PXCLASSCNT(2,2);
        tbl_summary.Tot_nonDend_area_sqr_um(sumcnter) = PXCLASSCNT(2,1);
        tbl_summary.Num_events(sumcnter) = 0;
        tbl_summary.Roi_area_um(sumcnter) = tmparea;
        tbl_summary.Roi_dist_soma_um(sumcnter) = tmpdist;
    else
        strct_onset_binary.(ROILIST{tmpidx}) = uint8(evINFO(tmpidx).binaryonset);
        strct_peaks_binary.(ROILIST{tmpidx}) = uint8(evINFO(tmpidx).binarypeaks);
        
        % Crop traces
        for iEv = 1:tmpelem
            [tmptrace] = crop_trace(tmpidx,iEv);
            icroptr = icroptr + 1;
            cell_cropped_trc{1, icroptr} = ROILIST{tmpidx};
            cell_cropped_trc{2, icroptr} = iEv;
            cell_cropped_trc{3, icroptr} = evINFO(tmpidx).onsetidx(iEv);
            cell_cropped_trc{4, icroptr} = FTIME;
            for iC = 1:size(tmptrace,1), cell_cropped_trc{iC+4, icroptr} = tmptrace(iC); end
        end
        
        for iEv = 1:tmpelem
            tbl_events1.Tot_eventidx(evcnter) = evcnter;
            tbl_events1.Roi(evcnter) = ROILIST(tmpidx);
            tbl_events1.Eventidx(evcnter) = iEv;
            tbl_events1.Onsetidx(evcnter) = evINFO(tmpidx).onsetidx(iEv);
            tbl_events1.Peakidx(evcnter) =  evINFO(tmpidx).peakidx(iEv);
            tbl_events1.Amplitude_raw(evcnter) = evINFO(tmpidx).amps(iEv);
            tbl_events1.Amplitude_dFoF(evcnter) = evINFO(tmpidx).amps_dFoF(iEv);
            tbl_events1.Baseline(evcnter) = evINFO(tmpidx).baselinevalues(iEv);
            tbl_events1.Baseline_s(evcnter) = evINFO(tmpidx).baseframes*FTIME;
            %         tbl_events1.Risetime = ;
            %         tbl_events1.Decaytime = ;
            evcnter = evcnter +1;
        end
        
        % Separate table for IEIs
        tmpieis = [evINFO(tmpidx).ieis];
        tbl_events2.Roi(ieicnter:ieicnter+numel(tmpieis)-1) = repelem(ROILIST(tmpidx), numel(tmpieis));
        tbl_events2.IEIs(ieicnter:ieicnter+numel(tmpieis)-1) = tmpieis;
        ieicnter = ieicnter+numel(tmpieis);
        
        % Summarized Event info
        tbl_summary.Roi(sumcnter) = ROILIST(tmpidx);
        tbl_summary.Tot_time_s(sumcnter) = TOTP*FTIME;
        tbl_summary.Tot_dend_area_sqr_um(sumcnter) = PXCLASSCNT(2,2);
        tbl_summary.Tot_nonDend_area_sqr_um(sumcnter) = PXCLASSCNT(2,1);
        tbl_summary.Num_events(sumcnter) = tmpelem;
        tbl_summary.Event_av_amplitude(sumcnter) = evINFO(tmpidx).avamp;
        tbl_summary.Event_av_amplitude_dFoF(sumcnter) = evINFO(tmpidx).avamp_dFoF;
        tbl_summary.Trace_range(sumcnter) = {strcat(num2str(evINFO(tmpidx).rangein),' - ', num2str(evINFO(tmpidx).rangeout))};
        tbl_summary.Trace_av_raw_intensity(sumcnter) = evINFO(tmpidx).average; % Within given range
        tbl_summary.Trace_sd_raw_intensity(sumcnter) = evINFO(tmpidx).sd; % Within given range
        tbl_summary.Event_av_iei(sumcnter) = evINFO(tmpidx).aviei;
        tbl_summary.Event_CV_iei(sumcnter) = evINFO(tmpidx).cviei;
        tbl_summary.Eventrate(sumcnter) = evINFO(tmpidx).eventrate;
        tbl_summary.Roi_area_um(sumcnter) = tmparea;
        tbl_summary.Roi_dist_soma_um(sumcnter) = tmpdist;
    end
    sumcnter = sumcnter+1;
end

% Add framecount
for iC = 1:size(cell_cropped_trc,1)-4, cell_cropped_trc{iC+4,1} = FTIME*iC; end

%% Write to files
disp('Writing files...');
writetable(struct2table(strct_onset_binary), binaryname, 'Sheet', 'Onset');
writetable(struct2table(strct_peaks_binary), binaryname, 'Sheet', 'Peaks');
writetable(tbl_events1, eventinfoname1);
writetable(tbl_events2, eventinfoname2);
writetable(tbl_summary, eventsummaryname);
xlswrite(croppedtrname, cell_cropped_trc);
save(strname, 'evINFO', '-v7.3');
disp('Writing finished!');

%% Cropped traces
    function [tmpcroptrace] = crop_trace(tmpidx, ev)
        tmpcroptrace = zeros(1,PREPOINTS+POSTPOINTS+1);
        worktrace = TRACEDATA(:,tmpidx);
        if SMTHWIN ~= 0, worktrace = smooth_data(worktrace,SMTHWIN); end
        onsetidx = evINFO(tmpidx).onsetidx(ev);
        baseline = evINFO(tmpidx).baselinevalues(ev);
        % Event preceding part
        if onsetidx-PREPOINTS < 1
            tmpcroptrace(1,PREPOINTS-onsetidx+2:PREPOINTS+1) = worktrace(1:onsetidx)';
        else
            tmpcroptrace(1,1:PREPOINTS+1) = worktrace(onsetidx-PREPOINTS:onsetidx)';
        end
        % Event following part
        if onsetidx+POSTPOINTS > TOTP
            tmpcroptrace(1,PREPOINTS+2:PREPOINTS+1+numel(worktrace)-onsetidx) = worktrace(onsetidx+1:end);
        else
            tmpcroptrace(1,PREPOINTS+2:end) = worktrace(onsetidx+1:onsetidx+POSTPOINTS);
        end
        % Delta F over F
        tmpcroptrace = (tmpcroptrace-baseline)./baseline;
        tmpcroptrace = tmpcroptrace'; % columns: Events, rows: Timepoints
    end
end