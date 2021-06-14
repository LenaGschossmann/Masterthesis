function save_files(sp, mode)

% Declare globally shared variables
global POSITIONFIG FONTSIZE PLOTRANGE WINSZ SMTHWIN...
    roiINFO traceINFO IMHW FTIME COMPOSITE2D EVDATA dFSHIFT dFWIN CRPPRESEC CRPPOSTSEC

diapos = [POSITIONFIG(1)+POSITIONFIG(3)/2 POSITIONFIG(2)+POSITIONFIG(4)/2 200 100];
diafig = figure('Position',diapos,'Name', 'Saving','toolbar', 'none', 'menu', 'none');
set(gca, 'Color', 'none'); axis off;
diatxtinfo = {'Files are saved', 'this may take a while...'};
diatxt = uicontrol('parent', diafig, 'style', 'text', 'position', [10 5 180 80], 'string', diatxtinfo, 'FONTSIZE', FONTSIZE+1);
pause(0.5);

% Prepare variables
prepts = round(CRPPRESEC/FTIME);
postpts = round(CRPPOSTSEC/FTIME);
totpts = prepts+postpts+1;
totevents = 0;
tblcnt = 1;
cntRoi = 1;
all_crp_dFoF = [];
all_crp_info = [];
deltawinsz  = round(dFWIN/FTIME);
timepoints = (1:size(COMPOSITE2D,1))'.*FTIME;
tbl_onset_binary.Time_s = timepoints;
tbl_peaks_binary.Time_s = timepoints;
dFoF_traces.Time_s =  timepoints;
event_info = struct('ROI',[],'Event',[],'Onset_frame',[],'Onset_time', [], 'Peak_frame',[],...
    'Peak_time',[],'Amp',[],'Amp_FoF',[],'Amp_dFoF',[],'IEI',[]);
roi_summary = struct('ROI',[],'Num_Events',[],'Amp_Mean',[],'Amp_FoF_Mean',[],...
    'Amp_dFoF_Mean',[], 'Amp_SD',[],'Amp_FoF_SD',[],'Amp_dFoF_SD',[],'Trc_Mean',[],...
    'Trc_SD',[],'Threshold',[],'IEI_Mean',[], 'CV_IEI',[], 'EvRate',[]);
rec_summary = struct('Num_ROIs',[],'ToT_Num_Events',[],'Amp_Mean',[],'Amp_FoF_Mean',[],'Amp_dFoF_Mean',[],...
    'Amp_SD',[],'Amp_FoF_SD',[],'Amp_dFoF_SD',[], 'IEI_Mean',[],'EvRate_Mean',[]);
trace_xls = strcat(sp, '\TRACES.xlsx'); if isa(trace_xls,'cell'), trace_xls=trace_xls{1};end
crptrace_xls = strcat(sp, '\CROPTRC.xlsx'); if isa(crptrace_xls,'cell'), crptrace_xls=crptrace_xls{1};end
binary_xls = strcat(sp, '\BINARY.xlsx'); if isa(binary_xls,'cell'), binary_xls=binary_xls{1};end
eventinfo_xls = strcat(sp, '\EVENTS.xlsx'); if isa(eventinfo_xls,'cell'), eventinfo_xls=eventinfo_xls{1};end
summary_xls = strcat(sp, '\SUMMARY.xlsx'); if isa(summary_xls,'cell'), summary_xls=summary_xls{1};end
ev_struct = strcat(sp, '\EVSTRCT.mat'); if isa(ev_struct,'cell'), ev_struct=ev_struct{1};end

answer = questdlg('Do you want to save the 2D matrix as .tif ?', 'Saving', 'Yes','No','Yes');
if strcmp(answer, 'Yes')
    % Save ROI-unrelated stuff
    set(diatxt,'string',[diatxtinfo,{'','...create .tif files.'}]);
    pause(0.2);
    vals = COMPOSITE2D;
    imwrite(uint16(vals), strcat(sp,'raw.tif'), 'tif');
    
    averaged = average_linescan(COMPOSITE2D, WINSZ);
    vals = averaged;
    imwrite(uint16(vals), strcat(sp,'AV_','bin_',num2str(WINSZ),'.tif'), 'tif');
    
    lin_averaged = averaged';
    [~, ~,alldFoF] = rollBase_dFoF(lin_averaged,deltawinsz,dFSHIFT, 'roll');
    writematrix(alldFoF', strcat(sp,'dFoF_','bin_',num2str(WINSZ),'_range_', num2str(PLOTRANGE(1)),'-',num2str(PLOTRANGE(2)), '.csv'));
end

% Determine ROIs to save
if strcmp(mode, 'all')
    saverois = 1:size(roiINFO,2);
    if isempty(roiINFO(end).position), saverois = saverois(1:end-1); end
else
    saveroisidx = [traceINFO([traceINFO(:).save] == 1).roiID];
    saverois = [];
    for iE = 1:numel(saveroisidx), saverois = [saverois find([roiINFO(:).ID] == saveroisidx(iE))]; end
end

% Check if data already exist
traceidx = zeros(size(saverois));

for iRoi = 1:numel(saverois)
    existbin = [traceINFO(:).roiID] == roiINFO(saverois(iRoi)).ID;
    if any(existbin)
        existidx = find(existbin);
        iEx = 1;
        while iEx <= numel(existidx)
            if traceINFO(existidx(iEx)).params{2,1} == WINSZ &&...
                  traceINFO(existidx(iEx)).params{3,1} == SMTHWIN
                traceidx(iRoi) = existidx(iEx);
                break;
            end
            iEx = iEx+1;
        end
    end
    if traceidx(iRoi) == 0
        [traceidx] = create_trc(saverois, iRoi, traceidx);
    end
end


for iRoi = 1:numel(saverois)
    tmproi = saverois(iRoi);
    save_idx = EVDATA{tmproi,7};
    roiid = roiINFO(tmproi).ID;
    
    set(diatxt,'string',[diatxtinfo,{'','...save ROI # ', num2str(roiid)}]);
    pause(0.2);
    
    avvals = traceINFO(traceidx(iRoi)).roi_av{1};
    dFoFvals = traceINFO(traceidx(iRoi)).dFoF_roi_av{1};
    FoFvals = traceINFO(traceidx(iRoi)).FoF_roi_av{1};
    dFvals = traceINFO(traceidx(iRoi)).dF_roi_av{1};
    timestamp = (1*FTIME:FTIME:IMHW(1)*FTIME)';
    
    scmarkedroi = traceINFO(traceidx(iRoi)).plotmarked{1};
    nframes = size(avvals,1);
    
    writeMtrx = [timestamp avvals FoFvals dFvals dFoFvals];
    writeMtrx = array2table(writeMtrx, 'VariableNames', {'frametime_s', 'ROI_average', 'ROI_FoF', 'ROI_dF','ROI_dFoF'});
    imwrite(uint16(scmarkedroi), strcat(sp,'ROI_',num2str(roiid),'.tif'), 'tif');
    writetable(writeMtrx, strcat(sp,'ROI_',num2str(roiid),'.csv'));
    
    roiname = sprintf('ROI_%i',saverois(iRoi));
    roi_summary(cntRoi).ROI = roiname;
    roi_summary(cntRoi).Trc_Mean = mean(avvals);
    roi_summary(cntRoi).Trc_SD = std(avvals);
    roi_summary(cntRoi).Trc_dFoF_SD = std(dFoFvals);
    roi_summary(cntRoi).Threshold = EVDATA{tmproi,10}(1);
    roi_summary(cntRoi).Trc_SNR = roi_summary(cntRoi).Amp_dFoF_Mean/roi_summary(cntRoi).Trc_dFoF_SD;
    
    %% Save events
    if ~isempty(save_idx)
        nevents = numel(save_idx);
        totevents = totevents+nevents;
        save_ptr = NaN(nevents,1);
        for iE = 1:nevents, save_ptr(iE) = find(save_idx(iE)==EVDATA{tmproi,3}); end
        
        dFoF_croptrc = NaN(nevents,totpts);
        for iE = 1:numel(save_idx)
            % Event preceding part
            if save_idx(iE)-prepts < 1
                dFoF_croptrc(iE,prepts-save_idx(iE)+2:prepts+1) = dFoFvals(1:save_idx(iE))';
            else
                dFoF_croptrc(iE,1:prepts+1) = dFoFvals(save_idx(iE)-prepts:save_idx(iE))';
            end
            % Event following part
            if save_idx(iE)+postpts > nframes
                dFoF_croptrc(iE,prepts+2:prepts+1+size(dFoFvals,1)-save_idx(iE)) = dFoFvals(save_idx(iE)+1:end);
            else
                dFoF_croptrc(iE,prepts+2:end) = dFoFvals(save_idx(iE)+1:save_idx(iE)+postpts);
            end
            crop_info{1,iE} = roiname;
            crop_info{2,iE} = iE;
            crop_info{3,iE} = save_idx(iE);
            crop_info{4,iE} = save_idx(iE)*FTIME;
        end
        all_crp_dFoF = [all_crp_dFoF dFoF_croptrc'];
        all_crp_info = [all_crp_info crop_info];
        
        % Binary traces
        tmponset = zeros(nframes,1); tmponset(save_idx,1) = 1;
        tmppeak = zeros(nframes,1); tmppeak(EVDATA{tmproi,4}(save_ptr),1) = 1;
        tbl_onset_binary.(roiname) = uint8(tmponset);
        tbl_peaks_binary.(roiname) = uint8(tmppeak);
        
        % Event info
        % Calculate IEIs
        tmpieis = diff(reshape(EVDATA{tmproi,7},[nevents 1]).*FTIME); tmpieis = [0; tmpieis];
        for iE =1:nevents
            event_info(tblcnt).ROI = tmproi;
            event_info(tblcnt).Event = iE;
            event_info(tblcnt).Onset_frame = EVDATA{tmproi,7}(iE);
            event_info(tblcnt).Onset_time = EVDATA{tmproi,7}(iE)*FTIME;
            event_info(tblcnt).Peak_frame = EVDATA{tmproi,4}(save_ptr(iE));
            event_info(tblcnt).Peak_time = EVDATA{tmproi,4}(save_ptr(iE))*FTIME;
            event_info(tblcnt).Amp = EVDATA{tmproi,5}(save_ptr(iE),3);
            event_info(tblcnt).Amp_FoF = EVDATA{tmproi,5}(save_ptr(iE),2);
            event_info(tblcnt).Amp_dFoF = EVDATA{tmproi,5}(save_ptr(iE),1);
            event_info(tblcnt).IEI = tmpieis(iE);
            tblcnt=tblcnt+1;
        end
        
        % ROI Event Summary
        roi_summary(cntRoi).Num_Events = nevents;
        roi_summary(cntRoi).Amp_Mean = mean(EVDATA{tmproi,5}(save_ptr,3));
        roi_summary(cntRoi).Amp_FoF_Mean = mean(EVDATA{tmproi,5}(save_ptr,2));
        roi_summary(cntRoi).Amp_dFoF_Mean = mean(EVDATA{tmproi,5}(save_ptr,1));
        roi_summary(cntRoi).Amp_SD = std(EVDATA{tmproi,5}(save_ptr,3));
        roi_summary(cntRoi).Amp_FoF_SD = std(EVDATA{tmproi,5}(save_ptr,2));
        roi_summary(cntRoi).Amp_dFoF_SD = std(EVDATA{tmproi,5}(save_ptr,1));
        if isempty(roi_summary(cntRoi).Amp_dFoF_Mean), roi_summary(cntRoi).Trc_SNR = NaN; end
        roi_summary(cntRoi).Event_confidence = EVDATA{tmproi,12}(1)/sum(EVDATA{tmproi,12});
        roi_summary(cntRoi).IEI_Mean = mean(tmpieis);
        roi_summary(cntRoi).CV_IEI = std(tmpieis)/mean(tmpieis);
        roi_summary(cntRoi).EvRate = nevents/(nframes*FTIME);
        cntRoi = cntRoi +1;
    else
        % ROI Event Summary
        roi_summary(cntRoi).Num_Events = 0;
        roi_summary(cntRoi).Amp_Mean = NaN;
        roi_summary(cntRoi).Amp_FoF_Mean = NaN;
        roi_summary(cntRoi).Amp_dFoF_Mean = NaN;
        roi_summary(cntRoi).Amp_SD = NaN;
        roi_summary(cntRoi).Amp_FoF_SD = NaN;
        roi_summary(cntRoi).Amp_dFoF_SD = NaN;
        if isempty(roi_summary(cntRoi).Amp_dFoF_Mean), roi_summary(cntRoi).Trc_SNR = NaN; end
        roi_summary(cntRoi).Event_confidence = NaN;
        roi_summary(cntRoi).IEI_Mean =NaN;
        roi_summary(cntRoi).CV_IEI = NaN;
        roi_summary(cntRoi).EvRate = NaN;
    end
    
end

%% Recoring Summary
rec_summary.Num_ROIs = numel(saverois);
rec_summary.ToT_Num_Events = totevents;
rec_summary.Num_Events_Mean = totevents/numel(saverois);
rec_summary.Amp_Mean = mean([roi_summary(:).Amp_Mean]);
rec_summary.Amp_FoF_Mean = mean([roi_summary(:).Amp_FoF_Mean]);
rec_summary.Amp_dFoF_Mean = mean([roi_summary(:).Amp_dFoF_Mean]);
rec_summary.Amp_SD = std([roi_summary(:).Amp_Mean]);
rec_summary.Amp_FoF_SD = std([roi_summary(:).Amp_FoF_Mean]);
rec_summary.Amp_dFoF_SD = std([roi_summary(:).Amp_dFoF_Mean]);
rec_summary.Mean_SNR = mean([roi_summary(:).Trc_SNR]);
rec_summary.Event_confidence = mean([roi_summary(:).Event_confidence]);
rec_summary.IEI_Mean = mean([roi_summary(:).IEI_Mean]);
rec_summary.EvRate_Mean = mean([roi_summary(:).EvRate]);
rec_summary.REC_time_s = nframes.*FTIME;

%% Save
% Cropped traces
if ~isempty(save_idx)
    writecell(all_crp_info, crptrace_xls, 'Sheet', 'TraceInfo'); % Info: Roi ID (row 1), Event ID (row 2), Onset index (row 3), Frametime (row 4)
    writematrix(all_crp_dFoF, crptrace_xls,'Sheet', 'dFoF_traces');
    
    % Save binary traces
    writetable(struct2table(tbl_onset_binary), binary_xls, 'Sheet', 'Onset');
    writetable(struct2table(tbl_peaks_binary), binary_xls, 'Sheet', 'Peak');
    
    % Event info
    writetable(struct2table(event_info), eventinfo_xls, 'Sheet', 'Event');
end

% Summary
writetable(struct2table(roi_summary), summary_xls, 'Sheet', 'ROIs');
writetable(struct2table(rec_summary), summary_xls, 'Sheet', 'Recording');

close(diafig);
okbox = msgbox('Files saved', '', 'modal');
end