function save_files(sp, fi, mode)

% Declare globally shared variables
global POSITIONFIG FONTSIZE CLIMRAW THRESHOLD PEAKTHRESHOLD SMTHWIN...
    figINFO roiINFO traceINFO IMHW FTIME COMPOSITE2D EVDATA dFWIN FNAME

diapos = [POSITIONFIG(1)+POSITIONFIG(3)/2 POSITIONFIG(2)+POSITIONFIG(4)/2 200 100];
diafig = figure('Position',diapos,'Name', 'Saving','toolbar', 'none', 'menu', 'none');
set(gca, 'Color', 'none'); axis off;
diatxtinfo = {'Files are saved', 'this may take a while...'};
diatxt = uicontrol('parent', diafig, 'style', 'text', 'position', [10 5 180 80], 'string', diatxtinfo, 'FONTSIZE', FONTSIZE+1);
pause(0.5);

% Prepare variables
pr = figINFO(fi).plotrange;
wsz = figINFO(fi).avwinsize;
prepts = round(0.5/FTIME);
postpts = round(1.5/FTIME);
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

if ~figINFO(fi).saved
    answer = questdlg('Do you want to save the 2D matrix as .tif ?', 'Saving', 'Yes','No','Yes');
    if strcmp(answer, 'Yes')
        % Save ROI-unrelated stuff
        averaged = average_linescan(COMPOSITE2D, wsz);
        lin_averaged = averaged';
        [~, alldFoF] = rollBase_dFoF(lin_averaged,deltawinsz,size(lin_averaged,2), 'roll');
        alldFoF = alldFoF';
        %         % .csv
        %         set(diatxt,'string',[diatxtinfo,{'','...create .csv files.'}]);
        %         pause(0.2);
        %         writematrix(COMPOSITE2D, strcat(sp,'raw.csv'));
        %         writematrix(averaged, strcat(sp,'AV_bin_',num2str(wsz), '.csv'));
        %         writematrix(dFoF, strcat(sp,'dFoF_bin_',num2str(wsz), '.csv'));
        % .tifs
        set(diatxt,'string',[diatxtinfo,{'','...create .tif files.'}]);
        pause(0.2);
        vals = COMPOSITE2D;
        imwrite(uint16(vals), strcat(sp,'raw_range.tif'), 'tif');
        vals = averaged;
        imwrite(uint16(vals), strcat(sp,'AV_','bin_',num2str(wsz),'.tif'), 'tif');
        vals = alldFoF;
        tic
        writematrix(vals, strcat(sp,'dFovF_','bin_',num2str(wsz),'_range_', num2str(pr(1)),'-',num2str(pr(2)), '.csv'));
        toc
        %         fid = fopen(strcat(sp,'dFovF_','bin_',num2str(wsz),'_range_', num2str(pr(1)),'-',num2str(pr(2)), '.txt'), 'w');
        %         fwrite(fid, vals, 'double');
        %         fclose(fid);
        %         imwrite(vals, strcat(sp,'dFovF_','bin_',num2str(wsz),'_range_', num2str(pr(1)),'-',num2str(pr(2)),'.tif'), 'tif');
        figINFO(fi).saved = true;
    end
else
    answer = 'No';
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
isline = false(size(saverois));

for iRoi = 1:numel(saverois)
    if roiINFO(iRoi).mode == 1, isline(iRoi) = true; end
    figid = figINFO(fi).IDs;
    existidx1 = [traceINFO(:).figID] == figid;
    existidx2 = [traceINFO(:).roiID] == roiINFO(saverois(iRoi)).ID;
    if any(existidx1 & existidx2)
        existidx = find(existidx1 & existidx2);
        iEx = 1;
        while iEx <= numel(existidx)
            if all(traceINFO(existidx(iEx)).fig_params{1,1} == pr) &&...
                    traceINFO(existidx(iEx)).fig_params{2,1} == wsz &&...
                    traceINFO(existidx(iEx)).fig_params{3,1} == SMTHWIN
                traceidx(iRoi) = existidx(iEx);
                break;
            end
            iEx = iEx+1;
        end
    end
    if traceidx(iRoi) == 0
        if isline(iRoi), tmppr = [1 IMHW(1)]; else, tmppr = pr; end
        [traceidx] = create_trc(figid, saverois, iRoi, traceidx, tmppr, wsz);
    end
end


for iRoi = 1:numel(saverois)
    tmproi = saverois(iRoi);
    save_idx = EVDATA{tmproi,7};
    if ~isempty(save_idx)
        roiid = roiINFO(tmproi).ID;
        roiname = sprintf('ROI_%i',saverois(iRoi));
        
        set(diatxt,'string',[diatxtinfo,{'','...save ROI # ', num2str(roiid)}]);
        pause(0.2);
        
        if isline(iRoi)
            avvals = traceINFO(traceidx(iRoi)).tot_binned_roi_av{1};
            dFoFvals = traceINFO(traceidx(iRoi)).tot_dFoF_roi_av{1};
            timestamp = traceINFO(traceidx(iRoi)).tot_timestamp{1};
        else
            avvals =  traceINFO(traceidx(iRoi)).binned_roi_av{1};
            dFoFvals = traceINFO(traceidx(iRoi)).dFoF_roi_av{1};
            timestamp = traceINFO(traceidx(iRoi)).timestamp{1};
        end
        scmarkedroi = traceINFO(traceidx(iRoi)).plotmarked{1};
        nframes = size(avvals,1);
        
        writeMtrx = [timestamp avvals dFoFvals];
        writeMtrx = array2table(writeMtrx, 'VariableNames', {'frametime_s', 'ROI_average', 'ROI_dFoF'});
        imwrite(uint16(scmarkedroi), strcat(sp,'AV_ROI_',num2str(roiid),'_binning_', num2str(wsz), '.tif'), 'tif');
        writetable(writeMtrx, strcat(sp,'ROI_',num2str(roiid),'.csv'));
        
        %% Save events
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
        roi_summary(cntRoi).ROI = roiname;
        roi_summary(cntRoi).Num_Events = nevents;
        roi_summary(cntRoi).Amp_Mean = mean(EVDATA{tmproi,5}(save_ptr,3));
        roi_summary(cntRoi).Amp_FoF_Mean = mean(EVDATA{tmproi,5}(save_ptr,2));
        roi_summary(cntRoi).Amp_dFoF_Mean = mean(EVDATA{tmproi,5}(save_ptr,1));
        roi_summary(cntRoi).Amp_SD = std(EVDATA{tmproi,5}(save_ptr,3));
        roi_summary(cntRoi).Amp_FoF_SD = std(EVDATA{tmproi,5}(save_ptr,2));
        roi_summary(cntRoi).Amp_dFoF_SD = std(EVDATA{tmproi,5}(save_ptr,1));
        roi_summary(cntRoi).Trc_Mean = mean(avvals);
        roi_summary(cntRoi).Trc_SD = std(avvals);
        roi_summary(cntRoi).Trc_dFoF_SD = std(dFoFvals);
        roi_summary(cntRoi).Threshold = EVDATA{tmproi,10}(1);
        roi_summary(cntRoi).Trc_SNR = roi_summary(cntRoi).Amp_dFoF_Mean/roi_summary(cntRoi).Trc_dFoF_SD;
        if isempty(roi_summary(cntRoi).Amp_dFoF_Mean), roi_summary(cntRoi).Trc_SNR = NaN; end
        roi_summary(cntRoi).Event_confidence = EVDATA{tmproi,12}(1)/sum(EVDATA{tmproi,12});
        roi_summary(cntRoi).IEI_Mean = mean(tmpieis);
        roi_summary(cntRoi).CV_IEI = std(tmpieis)/mean(tmpieis);
        roi_summary(cntRoi).EvRate = nevents/(nframes*FTIME);
        cntRoi = cntRoi +1;
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
writecell(all_crp_info, crptrace_xls, 'Sheet', 'TraceInfo'); % Info: Roi ID (row 1), Event ID (row 2), Onset index (row 3), Frametime (row 4)
writematrix(all_crp_dFoF, crptrace_xls,'Sheet', 'dFoF_traces');

% Save binary traces
writetable(struct2table(tbl_onset_binary), binary_xls, 'Sheet', 'Onset');
writetable(struct2table(tbl_peaks_binary), binary_xls, 'Sheet', 'Peak');

% Event info
writetable(struct2table(event_info), eventinfo_xls, 'Sheet', 'Event');

% Summary
writetable(struct2table(roi_summary), summary_xls, 'Sheet', 'ROIs');
writetable(struct2table(rec_summary), summary_xls, 'Sheet', 'Recording');

close(diafig);
okbox = msgbox('Files saved', '', 'modal');
end