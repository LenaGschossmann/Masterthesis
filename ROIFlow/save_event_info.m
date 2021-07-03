function [summary, roi_summary, synchro_dist_tbl] = save_event_info(eventdata, roiselection, roidata, params, path, filename, batchmode)

%% Unpack & Initialize
% roiselection = reshape(roiselection,[numel(roiselection) 1]);
save_rois = find(roiselection == 1);

nrois = numel(save_rois);
roinames = cell(nrois,1);
for i=1:nrois, roinames{i} = sprintf('ROI_%i',save_rois(i)); end
ft = roidata.frametime_s;
prepts = round(params.crp_trc_pre_s/ft);
postpts = round(params.crp_trc_post_s/ft);
totpts = prepts+postpts+1;
tot_events = 0;
tblcnt = 1;
all_crp_dFoF = [];
all_crp_FoF = [];
all_crp_info = [];
event_info = struct('ROI',[],'Event',[],'Onset_frame',[],'Peak_frame',[],'Peak_time',[],...
    'Amp',[],'Amp_FoF',[],'Peak_dFoF',[],'Baseline',[],'IEI',[],...
    'Subtract_value',[],'Offset',[],'Pmt_gain',[]);
roi_summary = struct('ROI',[], 'Num_Events',[],'Size_FiltKernel',[],...
    'Amp_Mean',[],'Amp_FoF_Mean',[],'Peak_dFoF_Mean', [],'Amp_SD',[],'Amp_FoF_SD',[],'Peak_dFoF_SD',[],...
    'Trc_Mean',[],'Trc_FoF_Mean',[],'Trc_dFoF_Mean',[],'Trc_SD',[],'Trc_FoF_SD',[],'Trc_dFoF_SD',[],...
    'Threshold_Slope',[],'Threshold_Amp',[],'Trc_SNR',[],'IEI_Mean',[],'CV_IEI',[],'EvRate_Hz',[],...
    'Area_um',[],'Area_ROI_PercOf_Struct',[],'Ctr_X_um',[],'Ctr_Y_um',[],...
    'Subtract_value',[],'Offset',[],'Pmt_gain',[]);
rec_summary = struct('Num_ROIs',[], 'ToT_Num_Events',[], 'Num_Events_Mean',[],...
    'Amp_Mean',[],'Amp_FoF_Mean',[],'Peak_dFoF_Mean', [],'Amp_SD',[],'Amp_FoF_SD',[],'Peak_dFoF_SD',[],...
    'Mean_SNR',[],'IEI_Mean',[],'EvRate_Hz_Mean',[],'REC_time_s',[],...
    'Area_um_Mean',[],'Area_ROI_PercOf_Struct',[],'Dend_Area_um',[],'Bg_Area_um',[],'FoV_Area_um',[],...
    'Events_per_Dend_Area',[],'Subtract_value',[],'Offset',[],'Pmt_gain',[]);
synchro_tbl = struct();
roi_names_cell = cell(nrois,1);
figsize = [30 30 1200 700];
cmap_n = max(cellfun(@max, roidata.roi_bounds(:,2)));
cmap = get_colormap([1 0 0],[1 1 0],[0 1 1],cmap_n);
uint8_im = (roidata.aip./max(roidata.aip,[],'all'))*255;

%% Outputnames
warning('off','MATLAB:xlswrite:AddSheet');
filename=filename(9:end-4);
savepath = fullfile(path,filename);
roipath = strcat(savepath,'\ROIs\'); if isa(roipath,'cell'), roipath=roipath{1};end
if ~isfolder(savepath), mkdir(savepath); end
if ~isfolder(roipath), mkdir(roipath); end

trace_xls = strcat(savepath, '\TRACES_', filename,'.xlsx'); if isa(trace_xls,'cell'), trace_xls=trace_xls{1};end
crptrace_xls = strcat(savepath, '\CROPTRC_', filename,'.xlsx'); if isa(crptrace_xls,'cell'), crptrace_xls=crptrace_xls{1};end
binary_xls = strcat(savepath, '\BINARY_', filename,'.xlsx'); if isa(binary_xls,'cell'), binary_xls=binary_xls{1};end
eventinfo_xls = strcat(savepath, '\EVENTS_', filename,'.xlsx'); if isa(eventinfo_xls,'cell'), eventinfo_xls=eventinfo_xls{1};end
summary_xls = strcat(savepath, '\SUMMARY_', filename,'.xlsx'); if isa(summary_xls,'cell'), summary_xls=summary_xls{1};end
dia_mtrx_xls = strcat(savepath, '\DIAMTRX_', filename,'.xlsx'); if isa(dia_mtrx_xls,'cell'), dia_mtrx_xls=dia_mtrx_xls{1};end
synchro_dist_xls = strcat(savepath, '\SYNCHRO_DIST_', filename,'.xlsx'); if isa(synchro_dist_xls,'cell'), synchro_dist_xls=synchro_dist_xls{1};end
synchro_tbl_xls = strcat(savepath, '\SYNCHRO_', filename,'.xlsx'); if isa(synchro_tbl_xls,'cell'), synchro_tbl_xls=synchro_tbl_xls{1};end
ev_struct = strcat(savepath, '\EVSTRCT_', filename,'.mat'); if isa(ev_struct,'cell'), ev_struct=ev_struct{1};end
aip_tiff =  strcat(savepath, '\AIP_', filename,'.tiff'); if isa(aip_tiff,'cell'), aip_tiff=aip_tiff{1};end
bg_tiff = strcat(savepath, '\BGMASK_', filename,'.tiff'); if isa(bg_tiff,'cell'), bg_tiff=bg_tiff{1};end
roi_ovw_png = strcat(savepath, '\ROIOVW_', filename,'.png'); if isa(roi_ovw_png,'cell'), roi_ovw_png=roi_ovw_png{1};end

%% Event-related
n_frames = size(roidata.traces,2);
timepoints = (1:n_frames)'.*ft;
tbl_onset_binary.Time_s = timepoints;
tbl_peaks_binary.Time_s = timepoints;
traces.Time_s = timepoints;
dFoF_traces.Time_s =  timepoints;
FoF_traces.Time_s =  timepoints;
syn_mtrx_Amp = NaN(nrois, n_frames);
syn_mtrx_Amp_FoF = NaN(nrois, n_frames);
syn_mtrx_Peak_dFoF = NaN(nrois, n_frames);
synchro_dist_tbl = table();
    
if ~isempty(save_rois)
    for iRoi = 1:nrois
        tmproi = save_rois(iRoi);
        curr_dFoFtrace =  eventdata(tmproi).filtrd_dFoFtrace';
        curr_FoFtrace = eventdata(tmproi).filtrd_FoFtrace';
        curr_trace = eventdata(tmproi).filtrd_trace';
        save_on_idx = eventdata(tmproi).onset_idx;
        save_peak_idx = eventdata(tmproi).peak_idx;
        n_events = numel(save_on_idx);    
        tot_events = tot_events+n_events;
        save_ptr = 1:numel(save_on_idx);
%         save_ptr = NaN(n_events,1);
%         for iE = 1:n_events, save_ptr(iE) = find(save_on_idx(iE)==eventdata{tmproi,3}); end
        n_baseline = roidata.baseline_frames;
        
        %% Traces
        traces.(roinames{iRoi}) = curr_trace;
        dFoF_traces.(roinames{iRoi}) = curr_dFoFtrace;
        FoF_traces.(roinames{iRoi}) =  curr_FoFtrace;
        
        %% Cropped traces
        dFoF_croptrc = NaN(totpts,n_events);
        FoF_croptrc = NaN(totpts,n_events);
        crop_info = cell(4,n_events);
        for iE = 1:numel(save_on_idx)
            % Event preceding part
            if save_on_idx(iE)-prepts < 1
                dFoF_croptrc(prepts-save_on_idx(iE)+2:prepts+1,iE) = curr_dFoFtrace(1:save_on_idx(iE));
                FoF_croptrc(prepts-save_on_idx(iE)+2:prepts+1,iE) = curr_FoFtrace(1:save_on_idx(iE));
            else
                dFoF_croptrc(1:prepts+1,iE) = curr_dFoFtrace(save_on_idx(iE)-prepts:save_on_idx(iE));
                FoF_croptrc(1:prepts+1,iE) = curr_FoFtrace(save_on_idx(iE)-prepts:save_on_idx(iE));
            end
            % Event following part
            if save_on_idx(iE)+postpts > n_frames
                dFoF_croptrc(prepts+2:prepts+1+size(curr_dFoFtrace,1)-save_on_idx(iE),iE) = curr_dFoFtrace(save_on_idx(iE)+1:end);
                FoF_croptrc(prepts+2:prepts+1+size(curr_dFoFtrace,1)-save_on_idx(iE),iE) = curr_FoFtrace(save_on_idx(iE)+1:end);
            else
                dFoF_croptrc(prepts+2:end,iE) = curr_dFoFtrace(save_on_idx(iE)+1:save_on_idx(iE)+postpts);
                FoF_croptrc(prepts+2:end,iE) = curr_FoFtrace(save_on_idx(iE)+1:save_on_idx(iE)+postpts);
            end
            crop_info{1,iE} = roinames{iRoi};
            crop_info{2,iE} = iE;
            crop_info{3,iE} = save_on_idx(iE);
            crop_info{4,iE} = save_on_idx(iE)*ft;
        end
        all_crp_dFoF = [all_crp_dFoF dFoF_croptrc];
        all_crp_FoF = [all_crp_FoF FoF_croptrc];
        all_crp_info = [all_crp_info crop_info];
        
        % Binary traces
        curr_on_binary = zeros(n_frames,1); curr_on_binary(save_on_idx(save_ptr),1) = 1;
        curr_peak_binary = zeros(n_frames,1); curr_peak_binary(save_peak_idx(save_ptr),1) = 1;
        tbl_onset_binary.(roinames{iRoi}) = uint8(curr_on_binary);
        tbl_peaks_binary.(roinames{iRoi}) = uint8(curr_peak_binary);
        
        % Event info
        % Calculate IEIs
        tmpieis = diff(reshape(eventdata(tmproi).onset_idx,[n_events 1]).*ft); tmpieis = [0; tmpieis];
        for iE =1:n_events
            tmpon = save_on_idx(save_ptr(iE));
            tmppk = save_peak_idx(save_ptr(iE));
            event_info(tblcnt).ROI = tmproi;
            event_info(tblcnt).Event = iE;
            event_info(tblcnt).Onset_frame = tmpon;
            event_info(tblcnt).Onset_time = tmpon*ft;
            event_info(tblcnt).Peak_frame = tmppk;
            event_info(tblcnt).Peak_time = tmppk*ft;
            event_info(tblcnt).Baseline = mean(curr_trace(tmpon - n_baseline:tmpon - 1));
            event_info(tblcnt).Amp = curr_trace(tmppk) - event_info(tblcnt).Baseline;
            event_info(tblcnt).Amp_FoF = curr_FoFtrace(tmppk) - 1;
%             event_info(tblcnt).Peak_dFoF = curr_dFoFtrace(tmppk);
            event_info(tblcnt).IEI = tmpieis(iE);
            event_info(tblcnt).Subtract_value = roidata.subtract_value;
            event_info(tblcnt).Offset = roidata.offset;
            event_info(tblcnt).Pmt_gain = roidata.pmt_gain;
            tblcnt=tblcnt+1;
            
%             syn_mtrx_Peak_dFoF(iRoi,tmpon) = curr_dFoFtrace(tmppk);
            syn_mtrx_Amp_FoF(iRoi,tmpon) = curr_FoFtrace(tmppk) - 1; 
            syn_mtrx_Amp(iRoi,tmpon) = curr_trace(tmppk) - mean(curr_trace(tmpon - n_baseline:tmpon - 1));
        end
        
        % ROI Event Summary
        roi_summary(iRoi).ROI = roinames{iRoi};
        roi_summary(iRoi).Num_Events = n_events;
        roi_summary(iRoi).Size_FiltKernel = eventdata(tmproi).size_FiltKernel;
        roi_summary(iRoi).Amp_Mean = mean(syn_mtrx_Amp(iRoi,:), 'omitnan');
        roi_summary(iRoi).Amp_FoF_Mean = mean(syn_mtrx_Amp_FoF(iRoi,:), 'omitnan');
%         roi_summary(iRoi).Peak_dFoF_Mean = mean(syn_mtrx_Peak_dFoF(iRoi,:),'omitnan');
        roi_summary(iRoi).Amp_SD = std(syn_mtrx_Amp(iRoi,:),'omitnan');
        roi_summary(iRoi).Amp_FoF_SD = std(syn_mtrx_Amp_FoF(iRoi,:), 'omitnan');
        roi_summary(iRoi).Peak_dFoF_SD = std(syn_mtrx_Peak_dFoF(iRoi,:), 'omitnan');
        roi_summary(iRoi).Trc_Mean = mean(curr_trace);
        roi_summary(iRoi).Trc_FoF_Mean = mean(curr_FoFtrace);
        roi_summary(iRoi).Trc_dFoF_Mean = mean(curr_dFoFtrace);
        roi_summary(iRoi).Trc_SD = std(curr_trace);
        roi_summary(iRoi).Trc_FoF_SD = std(curr_FoFtrace);
        roi_summary(iRoi).Trc_dFoF_SD = std(curr_dFoFtrace);
        roi_summary(iRoi).Threshold_Slope = eventdata(tmproi).FoF_1stDer_threshold;
        roi_summary(iRoi).Threshold_Amp = eventdata(tmproi).FoF_threshold;
        roi_summary(iRoi).Trc_SNR = mean(curr_FoFtrace(save_peak_idx(save_ptr)))/std(curr_FoFtrace);
        if isempty(roi_summary(iRoi).Amp_FoF_Mean), roi_summary(iRoi).Trc_SNR = NaN; end
        roi_summary(iRoi).IEI_Mean = mean(tmpieis);
        roi_summary(iRoi).CV_IEI = std(tmpieis)/mean(tmpieis);
        roi_summary(iRoi).EvRate_Hz = n_events*60/(n_frames*ft);
        roi_summary(iRoi).Area_um = roidata.roi_area_um(tmproi);
        roi_summary(iRoi).Area_ROI_PercOf_Struct = roidata.roi_area_um(tmproi)*100/roidata.dendritic_area_um;
        roi_summary(iRoi).Ctr_X_um = roidata.roi_centroids_um(tmproi,1);
        roi_summary(iRoi).Ctr_Y_um = roidata.roi_centroids_um(tmproi,2);
        roi_summary(iRoi).Subtract_value = roidata.subtract_value;
        roi_summary(iRoi).Offset = roidata.offset;
        roi_summary(iRoi).Pmt_gain = roidata.pmt_gain;
        
        % ROI snap
        f1=figure('Position',figsize,'Visible','off');
        imagesc(uint8(uint8_im));colormap('gray');axis('off');daspect([1 1 1]);
        boundary = roidata.roi_bounds{save_rois(iRoi),1}; boundary = boundary{1};
        hold on, plot(boundary(:,2), boundary(:,1), 'Color',cmap(:,roidata.roi_bounds{save_rois(iRoi),2})', 'LineWidth', 1)
        roi_png = strcat(savepath, '\ROI_',num2str(tmproi),'_', filename,'.png'); if isa(roi_png,'cell'), roi_png=roi_png{1};end
        saveas(f1,roi_png); close(f1);
        
        % Save ROI coordinates as .roi file
        roi_name = strcat('roi_',num2str(iRoi));
        writeImageJROI(boundary, 3, 0,0,0, roipath, roi_name);
        roi_names_cell{iRoi} = strcat(roipath, roi_name,'.roi');
    end
    
    %% Calculate ROI Distance and Synchronicity Matrix
    disp('Calculate distance matrix...');
    synchro_mtrx = ones(nrois,nrois);
    distance_mtrx = ones(nrois,nrois);
    fill_mtrx = true(nrois,nrois);
    roi_mtrx = cell(nrois,nrois);
    synchro_dist_mtrx = NaN(sum(1:nrois-1),6);
    cnt = 1;
    for iRoi1 = 1:nrois-1
        tmproi1 = save_rois(iRoi1);
        tmpbin1 = zeros(n_frames,1);
        save_ptr1 = 1:numel(eventdata(tmproi1).onset_idx);
        tmpbin1(eventdata(tmproi1).onset_idx(save_ptr1)) = 1;
        tmpx1 = roidata.roi_centroids_um(tmproi1,1);
        tmpy1 = roidata.roi_centroids_um(tmproi1,2);
        evr1 = roi_summary(iRoi1).EvRate_Hz;
        for iRoi2 = iRoi1+1:nrois
            % Distance
            tmproi2 = save_rois(iRoi2);
            tmpx2 = roidata.roi_centroids_um(tmproi2,1);
            tmpy2 = roidata.roi_centroids_um(tmproi2,2);
            dist = norm([tmpx1 tmpy1] -[tmpx2 tmpy2]);
            evr2 = roi_summary(iRoi2).EvRate_Hz;
            distance_mtrx(iRoi1,iRoi2) = dist;
            % Synchronicity expressed as Pearsons correlation coefficient
            save_ptr2 = 1:numel(eventdata(tmproi2).onset_idx);
            tmpbin2 = zeros(n_frames,1); tmpbin2(eventdata(tmproi2).onset_idx(save_ptr2)) = 1;
            av1 = mean(tmpbin1); av2 = mean(tmpbin2);
            r = sum((av1-tmpbin1).*(av2-tmpbin2)) / sqrt(sum((av1-tmpbin1).^2)*sum((av2-tmpbin2).^2));
            synchro_mtrx(iRoi1,iRoi2) =  r;
            fill_mtrx(iRoi1,iRoi2) = false;
            roi_mtrx{iRoi1,iRoi2} = strcat(roinames{iRoi1},'-',roinames{iRoi2});
            synchro_dist_mtrx(cnt,1) = tmproi1; synchro_dist_mtrx(cnt,2) = tmproi2;
            synchro_dist_mtrx(cnt,3) = r; synchro_dist_mtrx(cnt,4) = dist;
            synchro_dist_mtrx(cnt,5) = evr1; synchro_dist_mtrx(cnt,6) = evr2;
            cnt = cnt+1;
        end
    end
    trsp_mtrx = distance_mtrx';
    distance_mtrx(fill_mtrx) = trsp_mtrx(fill_mtrx);
    trsp_mtrx = synchro_mtrx';
    synchro_mtrx(fill_mtrx) = trsp_mtrx(fill_mtrx);
    trsp_mtrx = roi_mtrx';
    roi_mtrx(fill_mtrx) = trsp_mtrx(fill_mtrx);
    synchro_dist_tbl.ROI_1 = synchro_dist_mtrx(:,1);
    synchro_dist_tbl.ROI_2 = synchro_dist_mtrx(:,2);
    synchro_dist_tbl.Dist_um = synchro_dist_mtrx(:,4);
    synchro_dist_tbl.PearsonR = synchro_dist_mtrx(:,3);
    synchro_dist_tbl.EvRate_Hz_ROI_1 = synchro_dist_mtrx(:,5);
    synchro_dist_tbl.EvRate_Hz_ROI_2 = synchro_dist_mtrx(:,6);
    
    %% Calculate Synchronicity Histogram
    for iT = 1:n_frames
        evid = find(isnan(syn_mtrx_Amp(:,iT))==0);
        synchro_tbl(iT).Time_s = timepoints(iT);
        synchro_tbl(iT).Eventcount = numel(evid);
        synchro_tbl(iT).Amp_Mean = mean(syn_mtrx_Amp(evid,iT));
        synchro_tbl(iT).Amp_FoF_Mean = mean(syn_mtrx_Amp_FoF(evid,iT));
%         synchro_tbl(iT).Peak_dFoF_Mean  = mean(syn_mtrx_Peak_dFoF(evid,iT));
        synchro_tbl(iT).Amp_SD = std(syn_mtrx_Amp(evid,iT));
        synchro_tbl(iT).Amp_FoF_SD = std(syn_mtrx_Amp_FoF(evid,iT));
        synchro_tbl(iT).Peak_dFoF_SD = std(syn_mtrx_Peak_dFoF(evid,iT));
    end
    
    %% Recoring Summary
    rec_summary.Num_ROIs = nrois;
    rec_summary.ToT_Num_Events = tot_events;
    rec_summary.Num_Events_Mean = tot_events/nrois;
    rec_summary.Amp_Mean = mean([roi_summary(:).Amp_Mean], 'omitnan');
    rec_summary.Amp_FoF_Mean = mean([roi_summary(:).Amp_FoF_Mean], 'omitnan');
%     rec_summary.Peak_dFoF_Mean = mean([roi_summary(:).Peak_dFoF_Mean], 'omitnan');
    rec_summary.Amp_SD = std([roi_summary(:).Amp_Mean], 'omitnan');
    rec_summary.Amp_FoF_SD = std([roi_summary(:).Amp_FoF_Mean], 'omitnan');
%     rec_summary.Peak_dFoF_SD = std([roi_summary(:).Peak_dFoF_Mean], 'omitnan');
    rec_summary.Mean_SNR = mean([roi_summary(:).Trc_SNR], 'omitnan');
    rec_summary.IEI_Mean = mean([roi_summary(:).IEI_Mean], 'omitnan');
    rec_summary.EvRate_Hz_Mean = mean([roi_summary(:).EvRate_Hz], 'omitnan');
    rec_summary.REC_time_s = n_frames.*ft;
    rec_summary.Area_um_Mean = mean([roi_summary(:).Area_um], 'omitnan');
    rec_summary.Area_ROI_PercOf_Struct = mean([roi_summary(:).Area_ROI_PercOf_Struct], 'omitnan');
    rec_summary.Dend_Area_um = roidata.dendritic_area_um;
    rec_summary.Bg_Area_um = roidata.background_area_um;
    rec_summary.FoV_Area_um = roidata.fov_area_um;
    rec_summary.Events_per_Dend_Area = tot_events/roidata.dendritic_area_um;
    rec_summary.Subtract_value = roidata.subtract_value;
    rec_summary.Offset = roidata.offset;
    rec_summary.Pmt_gain = roidata.pmt_gain;
    
    %% Save
    % Traces
    writetable(struct2table(dFoF_traces), trace_xls, 'Sheet', 'dFoF_traces');
    writetable(struct2table(FoF_traces), trace_xls, 'Sheet', 'FoF_traces');
    writetable(struct2table(traces), trace_xls, 'Sheet', 'Traces');
    
    % Cropped traces
    writecell(all_crp_info, crptrace_xls, 'Sheet', 'TraceInfo'); % Info: Roi ID (row 1), Event ID (row 2), Onset index (row 3), Frametime (row 4)
    writematrix(all_crp_dFoF, crptrace_xls,'Sheet', 'dFoF_traces');
    writematrix(all_crp_FoF, crptrace_xls,'Sheet', 'FoF_traces');
    
    % Save binary traces
    writetable(struct2table(tbl_onset_binary), binary_xls, 'Sheet', 'Onset');
    writetable(struct2table(tbl_peaks_binary), binary_xls, 'Sheet', 'Peak');
    
    % Event info
    writetable(struct2table(event_info), eventinfo_xls, 'Sheet', 'Event');
    
    % Summary
    writetable(struct2table(roi_summary), summary_xls, 'Sheet', 'ROIs');
    writetable(struct2table(rec_summary), summary_xls, 'Sheet', 'Recording');
    
    % Diagonal Matrices (Distance, Synchronicity, ROI Names)
    writecell(roi_mtrx, dia_mtrx_xls,'Sheet', 'ROIs');
    writematrix(distance_mtrx, dia_mtrx_xls, 'Sheet', 'Dist_um');
    writematrix(synchro_mtrx, dia_mtrx_xls, 'Sheet', 'PearsonR');
    
    % Synchronicity Table
    writetable(struct2table(synchro_tbl), synchro_tbl_xls, 'Sheet', 'Synchronicity');
    writetable(synchro_dist_tbl, synchro_dist_xls);
    
    % Zip .roi files
    roi_zip = strcat(savepath,'\ROIs.zip'); if isa(roi_zip,'cell'), roi_zip=roi_zip{1};end
    zip(roi_zip,roi_names_cell);
    roipath = roipath(1:end-1);
    rmdir(roipath, 's');
    
end

% Write AIP w/o ROIs as .tiff
imwrite(uint8(uint8_im), aip_tiff);

% Write Background mask as .tiff (background = black/0)
imsize = size(roidata.aip);
bg_im = zeros(imsize(1), imsize(2),1);
bg_im(roidata.background_mask==0) = 255;
imwrite(uint8(bg_im), bg_tiff);

% Write ROI overview as .png and event cell as struct
if ~isempty(save_rois)
    % ROI overview
    f1=figure('Position',figsize,'Visible','off');
    imagesc(uint8(uint8_im));colormap('gray');axis('off');daspect([1 1 1]);
    for iRoi = 1:nrois
        boundary = roidata.roi_bounds{save_rois(iRoi),1}; boundary = boundary{1};
        hold on, plot(boundary(:,2), boundary(:,1), 'Color',cmap(:,roidata.roi_bounds{save_rois(iRoi),2})', 'LineWidth', 1)
    end
    saveas(f1,roi_ovw_png); close(f1);
    
    % Cell with events
    eventinfo = struct();
    eventinfo.eventdata = eventdata;
    save(ev_struct, 'eventinfo', '-v7.3');
end

disp('Files saved');

%% Output
if batchmode
    summary = rec_summary;
else
    summary = [];
end
end