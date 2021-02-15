function save_files(sp, fi, mode)

% Declare globally shared variables
global POSITIONFIG FONTSIZE CLIMRAW THRESHOLD PEAKTHRESHOLD SMTHWIN...
    figINFO roiINFO traceINFO IMHW FTIME COMPOSITE2D

diapos = [POSITIONFIG(1)+POSITIONFIG(3)/2 POSITIONFIG(2)+POSITIONFIG(4)/2 200 100];
diafig = figure('Position',diapos,'Name', 'Saving','toolbar', 'none', 'menu', 'none');
set(gca, 'Color', 'none'); axis off;
diatxtinfo = {'Files are saved', 'this may take a while...'};
diatxt = uicontrol('parent', diafig, 'style', 'text', 'position', [10 5 180 80], 'string', diatxtinfo, 'FONTSIZE', FONTSIZE+1);
pause(0.5);

% Get parameters
pr = figINFO(fi).plotrange;
wsz = figINFO(fi).avwinsize;

if ~figINFO(fi).saved
    answer = questdlg('Do you want to save the 2D matrix as .tif ?', 'Saving', 'Yes','No','Yes');
    if strcmp(answer, 'Yes')
        % Save ROI-unrelated stuff
        averaged = average_linescan(COMPOSITE2D, wsz);
        [~, dFoF] = rollBase_dFoF(averaged);
%         % .csv
%         set(diatxt,'string',[diatxtinfo,{'','...create .csv files.'}]);
%         pause(0.2);
%         writematrix(COMPOSITE2D, strcat(sp,'raw.csv'));
%         writematrix(averaged, strcat(sp,'AV_bin_',num2str(wsz), '.csv'));
%         writematrix(dFoF, strcat(sp,'dFoF_bin_',num2str(wsz), '.csv'));
        % .tifs
        set(diatxt,'string',[diatxtinfo,{'','...create .tif files.'}]);
        pause(0.2);
        vals = COMPOSITE2D(pr(1):pr(2),:);
        imwrite(uint16(vals), strcat(sp,'raw_range_', num2str(pr(1)),'-',num2str(pr(2)),'.tif'), 'tif');
        vals = averaged(pr(1):pr(2),:);
        imwrite(uint16(vals), strcat(sp,'AV_','bin_',num2str(wsz),'_range_', num2str(pr(1)),'-',num2str(pr(2)),'.tif'), 'tif');
        vals = dFoF(pr(1):pr(2),:);
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
trcexists = false(size(saverois));
isline = false(size(saverois));

for iRoi = 1:numel(saverois)
    figid = figINFO(fi).IDs;
    existidx1 = [traceINFO(:).figID] == figid;
    existidx2 = [traceINFO(:).roiID] == roiINFO(saverois(iRoi)).ID;
    if any(existidx1 & existidx2)
        existidx = find(existidx1 & existidx2);
        iEx = 1;
        while iEx <= numel(existidx)
            if all(traceINFO(existidx(iEx)).fig_params{1,1} == pr) &&...
                    traceINFO(existidx(iEx)).fig_params{2,1} == wsz &&...
                    traceINFO(existidx(iEx)).fig_params{3,1} == SMTHWIN &&...
                    traceINFO(existidx(iEx)).fig_params{4,1} == THRESHOLD &&...
                    traceINFO(existidx(iEx)).fig_params{5,1} == PEAKTHRESHOLD % check plotrange & binning
                trcexists(iRoi) = true;
                traceidx(iRoi) = existidx(iEx);
                break;
            end
            iEx = iEx+1;
        end
    end
    if roiINFO(iRoi).mode == 1, isline(iRoi) = true; end
end

if any(~trcexists)
    averaged = average_linescan(COMPOSITE2D, wsz);
end

for iRoi = 1:numel(saverois)
    roiid = roiINFO(saverois(iRoi)).ID;
    set(diatxt,'string',[diatxtinfo,{'','...save ROI # ', num2str(roiid)}]);
    pause(0.2);
    if ~trcexists(iRoi)
        if isline(iRoi), tmppr = [1 IMHW(1)]; else, tmppr = pr; end
        val = averaged(tmppr(1):tmppr(2),:);
        scmarkedroi = mark_rois(val,tmppr, saverois(iRoi));
        [avvals,dFoFvals,yrange,allvals,~] = calc_roi_av_trace(saverois(iRoi), averaged, tmppr);
        [evinfo, smoothed] = get_trc_params(allvals,[tmppr(1)+yrange(1)-1 tmppr(1)+yrange(2)-1], [],[]);
        crossings = evinfo.crossings;
        supraT = evinfo.suprathreshold;
        peaks = evinfo.peaks;
        timestamp = ((tmppr(1)+yrange(1)-1)*FTIME:FTIME:(tmppr(1)+yrange(2)-1)*FTIME)';
    else
        if isline(iRoi)
            avvals = traceINFO(traceidx(iRoi)).tot_binned_roi_av{1};
            dFoFvals = traceINFO(traceidx(iRoi)).tot_dFoF_roi_av{1};
            smoothed = traceINFO(traceidx(iRoi)).tot_smoothed{1};
            crossings = traceINFO(traceidx(iRoi)).tot_events.crossings;
            supraT = traceINFO(traceidx(iRoi)).tot_events.suprathreshold;
            peaks = traceINFO(traceidx(iRoi)).tot_events.peaks;
            timestamp = traceINFO(traceidx(iRoi)).tot_timestamp{1};
            evinfo =  traceINFO(traceidx(iRoi)).tot_events;
        else
            avvals =  traceINFO(traceidx(iRoi)).binned_roi_av{1};
            dFoFvals = traceINFO(traceidx(iRoi)).dFoF_roi_av{1};
            smoothed = traceINFO(traceidx(iRoi)).smoothed{1};
            crossings = traceINFO(traceidx(iRoi)).events.crossings;
            supraT = traceINFO(traceidx(iRoi)).events.suprathreshold;
            peaks = traceINFO(traceidx(iRoi)).events.peaks;
            timestamp = traceINFO(traceidx(iRoi)).timestamp{1};
            evinfo =  traceINFO(traceidx(iRoi)).events;
        end
        scmarkedroi = traceINFO(traceidx(iRoi)).plotmarked{1};
    end
    writeMtrx = [timestamp avvals smoothed dFoFvals supraT crossings peaks];
    writeMtrx = array2table(writeMtrx, 'VariableNames', {'frametime_s', 'ROI_average', 'ROI_smth_average','ROI_dFoF', 'threshold_binary','crossing_binary', 'peak_binary'});
    imwrite(uint16(scmarkedroi), strcat(sp,'AV_ROI_',num2str(roiid),'_binning_', num2str(wsz), '.tif'), 'tif');
    writetable(writeMtrx, strcat(sp,'ROI_',num2str(roiid),'_values_binning_', num2str(wsz), '_smth_', num2str(SMTHWIN), '_thresh_', num2str(THRESHOLD),'_threshpk_', num2str(PEAKTHRESHOLD),'.csv'));
    
    % Save event info as table
    evTable = table();
    if isempty(evinfo.crossidx)
        evTable.crossing_idx = NaN;
        evTable.peaks_idx = NaN;
        evTable.peak_amps = NaN;
        evTable.ieis = NaN;
    else
        evTable.crossing_idx = evinfo.crossidx;
        evTable.peaks_idx = evinfo.peakidx;
        evTable.peak_amps = evinfo.amps;
        evTable.ieis = [evinfo.ieis; NaN];
    end
    lt = size(evinfo.crossidx,1);
    if lt > 1, filler = repelem({''},lt-1)'; else, filler = []; end
    evTable.roi_start_x = [roiINFO(saverois(iRoi)).position(1); filler];
    evTable.roi_end_x = [roiINFO(saverois(iRoi)).position(1)+roiINFO(saverois(iRoi)).position(3); filler];
    evTable.roi_start_y = [roiINFO(saverois(iRoi)).position(2); filler];
    evTable.roi_end_y = [roiINFO(saverois(iRoi)).position(2)+roiINFO(saverois(iRoi)).position(4); filler];
    evTable.threshold = [evinfo.threshold; filler];
    evTable.event_type = [evinfo.eventtype; filler];
    evTable.av_eventrate_Hz = [evinfo.eventrate; filler];
    evTable.av_intereventinterval_s = [evinfo.aviei; filler];
    evTable.av_peak_amp = [evinfo.avamp;filler];
    evTable.cv_iei = [evinfo.cviei; filler];
    writetable(evTable, strcat(sp,'ROI_',num2str(roiid),'_eventinfo_binning_', num2str(wsz), '_smth_', num2str(SMTHWIN), '_thresh_', num2str(THRESHOLD),'_threshpk_', num2str(PEAKTHRESHOLD),'.csv'));
end
close(diafig);
okbox = msgbox('Files saved', '', 'modal');
end