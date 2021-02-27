function save_event_info(eventdata, roidata, path, filename)

%% Unpack & Initialize
save_rois = find(cellfun(@isempty, eventdata(:,7))==0);

nrois = numel(save_rois);
roinames = cell(nrois,1);
n_frames = size(roidata.traces,2);
for i=1:nrois, roinames{i} = sprintf('ROI_%i',save_rois(i)); end
ft = roidata.frametime_s;
prepts = round(0.5/ft);
postpts = round(1.5/ft);
totpts = prepts+postpts+1;
totevents = 0;
tblcnt = 1;
all_crp_dFoF = [];
all_crp_FoF = [];
all_crp_info = [];
event_info = struct('ROI',[],'Event',[],'Onset_frame',[],'Onset_time', [], 'Peak_frame',[],...
    'Peak_time',[],'Amp',[],'Amp_FoF',[],'Amp_dFoF',[],'IEI',[]);
roi_summary = struct('ROI',[],'Num_Events',[],'Amp_Mean',[],'Amp_FoF_Mean',[],...
    'Amp_dFoF_Mean',[], 'Amp_SD',[],'Amp_FoF_SD',[],'Amp_dFoF_SD',[],'Trc_Mean',[],...
    'Trc_SD',[],'SD_Threshold',[],'IEI_Mean',[], 'CV_IEI',[], 'EvRate',[],...
    'Area_um',[],'Ctr_X_um',[],'Ctr_Y_um',[],'StructNorm_Area',[]);
rec_summary = struct('Num_ROIs',[],'ToT_Num_Events',[],'Amp_Mean',[],'Amp_FoF_Mean',[],'Amp_dFoF_Mean',[],...
    'Amp_SD',[],'Amp_FoF_SD',[],'Amp_dFoF_SD',[], 'IEI_Mean',[],'EvRate_Mean',[],...
    'Area_um_Mean',[],'Dend_Area_um',[],'Bg_Area_um',[],'FoV_Area_um',[],'StructNorm_Area',[],'Events_per_Dend_Area',[]);
synchro_tbl = struct('Time_s',[], 'Eventcount',[], 'Amp_Mean',[],'Amp_FoF_Mean',[],'Amp_dFoF_Mean',[],...
    'Amp_SD',[],'Amp_FoF_SD',[],'Amp_dFoF_SD',[]);
figsize = [30 30 1200 700];
cmap_n = max(cellfun(@max, roidata.roi_bounds(:,2)));
cmap = get_colormap([1 0 0],[1 1 0],[0 1 1],cmap_n);
uint8_im = (roidata.aip./max(roidata.aip,[],'all'))*255;

%% Outputnames
warning('off','MATLAB:xlswrite:AddSheet');
filename=filename(9:end-4);
savepath = fullfile(path,filename);
if ~isdir(savepath), mkdir(savepath); end

trace_xls = strcat(savepath, '\TRACES_', filename,'.xlsx'); if isa(trace_xls,'cell'), trace_xls=trace_xls{1};end
crptrace_xls = strcat(savepath, '\CROPTRC_', filename,'.xlsx'); if isa(crptrace_xls,'cell'), crptrace_xls=crptrace_xls{1};end
binary_xls = strcat(savepath, '\BINARY_', filename,'.xlsx'); if isa(binary_xls,'cell'), binary_xls=binary_xls{1};end
eventinfo_xls = strcat(savepath, '\EVENTS_', filename,'.xlsx'); if isa(eventinfo_xls,'cell'), eventinfo_xls=eventinfo_xls{1};end
summary_xls = strcat(savepath, '\SUMMARY_', filename,'.xlsx'); if isa(summary_xls,'cell'), summary_xls=summary_xls{1};end
bg_png = strcat(savepath, '\BGMASK_', filename,'.png'); if isa(bg_png,'cell'), bg_png=bg_png{1};end
roi_ovw_png = strcat(savepath, '\ROIOVW_', filename,'.png'); if isa(roi_ovw_png,'cell'), roi_ovw_png=roi_ovw_png{1};end
dia_mtrx_xls = strcat(savepath, '\DIAMTRX_', filename,'.xlsx'); if isa(dia_mtrx_xls,'cell'), dia_mtrx_xls=dia_mtrx_xls{1};end
synchro_dist_xls = strcat(savepath, '\SYNCHRO_DIST_', filename,'.xlsx'); if isa(synchro_dist_xls,'cell'), synchro_dist_xls=synchro_dist_xls{1};end
synchro_tbl_xls = strcat(savepath, '\SYNCHRO_', filename,'.xlsx'); if isa(synchro_tbl_xls,'cell'), synchro_tbl_xls=synchro_tbl_xls{1};end
ev_struct = strcat(savepath, '\EVSTRCT_', filename,'.mat'); if isa(ev_struct,'cell'), ev_struct=ev_struct{1};end

%% Event-related
nframes = size(roidata.traces,2);
timepoints = (1:nframes)'.*ft;
tbl_onset_binary.Time_s = timepoints;
tbl_peaks_binary.Time_s = timepoints;
traces.Time_s = timepoints;
dFoF_traces.Time_s =  timepoints;
FoF_traces.Time_s =  timepoints;
syn_mtrx_Amp = NaN(nrois, nframes);
syn_mtrx_Amp_FoF = NaN(nrois, nframes);
syn_mtrx_Amp_dFoF = NaN(nrois, nframes);

if ~isempty(save_rois)
    for iRoi = 1:nrois
        tmproi = save_rois(iRoi);
        tmpdFoFtrace = roidata.dFoF_traces(tmproi,:);
        tmpFoFtrace = roidata.FoF_traces(tmproi,:);
        tmptrace = roidata.traces(tmproi,:);
        save_idx = eventdata{tmproi,7};
        nevents = numel(save_idx);
        totevents = totevents+nevents;
        save_ptr = NaN(nevents,1);
        for iE = 1:nevents, save_ptr(iE) = find(save_idx(iE)==eventdata{tmproi,3}); end
        
        %% Traces
        traces.(roinames{iRoi}) = roidata.traces(tmproi,:)';
        dFoF_traces.(roinames{iRoi}) =  roidata.dFoF_traces(tmproi,:)';
        FoF_traces.(roinames{iRoi}) =  roidata.FoF_traces(tmproi,:)';
        
        %% Cropped traces
        dFoF_croptrc = NaN(nevents,totpts);
        FoF_croptrc = NaN(nevents,totpts);
        crop_info = cell(4,nevents);
        for iE = 1:numel(save_idx)
            % Event preceding part
            if save_idx(iE)-prepts < 1
                dFoF_croptrc(iE,prepts-save_idx(iE)+2:prepts+1) = tmpdFoFtrace(1:save_idx(iE))';
                FoF_croptrc(iE,prepts-save_idx(iE)+2:prepts+1) = tmpFoFtrace(1:save_idx(iE))';
            else
                dFoF_croptrc(iE,1:prepts+1) = tmpdFoFtrace(save_idx(iE)-prepts:save_idx(iE))';
                FoF_croptrc(iE,1:prepts+1) = tmpFoFtrace(save_idx(iE)-prepts:save_idx(iE))';
            end
            % Event following part
            if save_idx(iE)+postpts > nframes
                dFoF_croptrc(iE,prepts+2:prepts+1+size(tmpdFoFtrace,2)-save_idx(iE)) = tmpdFoFtrace(save_idx(iE)+1:end);
                FoF_croptrc(iE,prepts+2:prepts+1+size(tmpdFoFtrace,2)-save_idx(iE)) = tmpFoFtrace(save_idx(iE)+1:end);
            else
                dFoF_croptrc(iE,prepts+2:end) = tmpdFoFtrace(save_idx(iE)+1:save_idx(iE)+postpts);
                FoF_croptrc(iE,prepts+2:end) = tmpFoFtrace(save_idx(iE)+1:save_idx(iE)+postpts);
            end
            crop_info{1,iE} = roinames{iRoi};
            crop_info{2,iE} = iE;
            crop_info{3,iE} = save_idx(iE);
            crop_info{4,iE} = save_idx(iE)*ft;
        end
        all_crp_dFoF = [all_crp_dFoF dFoF_croptrc'];
        all_crp_FoF = [all_crp_FoF FoF_croptrc'];
        all_crp_info = [all_crp_info crop_info];
        
        % Binary traces
        tmponset = zeros(nframes,1); tmponset(save_idx,1) = 1;
        tmppeak = zeros(nframes,1); tmppeak(eventdata{tmproi,4}(save_ptr),1) = 1;
        tbl_onset_binary.(roinames{iRoi}) = uint8(tmponset);
        tbl_peaks_binary.(roinames{iRoi}) = uint8(tmppeak);
        
        % Event info
        % Calculate IEIs
        tmpieis = diff(reshape(eventdata{tmproi,7},[nevents 1]).*ft); tmpieis = [0; tmpieis];
        for iE =1:nevents
            event_info(tblcnt).ROI = tmproi;
            event_info(tblcnt).Event = iE;
            event_info(tblcnt).Onset_frame = eventdata{tmproi,7}(iE);
            event_info(tblcnt).Onset_time = eventdata{tmproi,7}(iE)*ft;
            event_info(tblcnt).Peak_frame = eventdata{tmproi,4}(save_ptr(iE));
            event_info(tblcnt).Peak_time = eventdata{tmproi,4}(save_ptr(iE))*ft;
            event_info(tblcnt).Amp = eventdata{tmproi,5}(save_ptr(iE),3);
            event_info(tblcnt).Amp_FoF = eventdata{tmproi,5}(save_ptr(iE),2);
            event_info(tblcnt).Amp_dFoF = eventdata{tmproi,5}(save_ptr(iE),1);
            event_info(tblcnt).IEI = tmpieis(iE);
            tblcnt=tblcnt+1;
            
            syn_mtrx_Amp_dFoF(iRoi,eventdata{tmproi,7}(iE)) = eventdata{tmproi,5}(save_ptr(iE),1);
            syn_mtrx_Amp_FoF(iRoi,eventdata{tmproi,7}(iE)) = eventdata{tmproi,5}(save_ptr(iE),2);
            syn_mtrx_Amp(iRoi,eventdata{tmproi,7}(iE)) = eventdata{tmproi,5}(save_ptr(iE),3);
        end
        
        % ROI Event Summary
        roi_summary(iRoi).ROI = roinames{iRoi};
        roi_summary(iRoi).Num_Events = nevents;
        roi_summary(iRoi).Amp_Mean = mean(eventdata{tmproi,5}(save_ptr,3));
        roi_summary(iRoi).Amp_FoF_Mean = mean(eventdata{tmproi,5}(save_ptr,2));
        roi_summary(iRoi).Amp_dFoF_Mean = mean(eventdata{tmproi,5}(save_ptr,1));
        roi_summary(iRoi).Amp_SD = std(eventdata{tmproi,5}(save_ptr,3));
        roi_summary(iRoi).Amp_FoF_SD = std(eventdata{tmproi,5}(save_ptr,2));
        roi_summary(iRoi).Amp_dFoF_SD = std(eventdata{tmproi,5}(save_ptr,1));
        roi_summary(iRoi).Trc_Mean = mean(tmptrace);
        roi_summary(iRoi).Trc_SD = std(tmptrace);
        roi_summary(iRoi).Trc_dFoF_SD = std(tmpdFoFtrace);
        roi_summary(iRoi).SD_Threshold = eventdata{tmproi,9}(1);
        roi_summary(iRoi).Trc_SNR = roi_summary(iRoi).Amp_dFoF_Mean/roi_summary(iRoi).Trc_dFoF_SD;
        if isempty(roi_summary(iRoi).Amp_dFoF_Mean), roi_summary(iRoi).Trc_SNR = NaN; end
        roi_summary(iRoi).Event_confidence = eventdata{tmproi,12}(1)/sum(eventdata{tmproi,12});
        roi_summary(iRoi).IEI_Mean = mean(tmpieis);
        roi_summary(iRoi).CV_IEI = std(tmpieis)/mean(tmpieis);
        roi_summary(iRoi).EvRate = nevents/(nframes*ft);
        roi_summary(iRoi).Area_um = roidata.roi_area_um(tmproi);
        roi_summary(iRoi).StructNorm_Area = roidata.roi_area_um(tmproi)/roidata.dendritic_area_um;
        roi_summary(iRoi).Ctr_X_um = roidata.roi_centroids_um(tmproi,1);
        roi_summary(iRoi).Ctr_Y_um = roidata.roi_centroids_um(tmproi,2);
        
        % ROI snap
        f1=figure('Position',figsize,'Visible','off');
        imagesc(uint8(uint8_im));colormap('gray');axis('off');daspect([1 1 1]);
        B = roidata.roi_bounds{save_rois(iRoi),1};
        for k = 1:length(B)
            boundary = B{k};
            hold on, plot(boundary(:,2), boundary(:,1), 'Color',cmap(:,roidata.roi_bounds{save_rois(iRoi),2})', 'LineWidth', 1)
        end
        roi_png = strcat(savepath, '\ROI_',num2str(tmproi),'_', filename,'.png'); if isa(roi_png,'cell'), roi_png=roi_png{1};end
        saveas(f1,roi_png); close(f1);
    end
    
    %% Calculate ROI Distance and Synchronicity Matrix
    disp('Calculate distance matrix...');
    synchro_mtrx = ones(nrois,nrois);
    distance_mtrx = ones(nrois,nrois);
    fill_mtrx = true(nrois,nrois);
    roi_mtrx = cell(nrois,nrois);
    synchro_dist_mtrx = NaN(sum(1:nrois-1),2);
    cnt = 1;
    for iRoi1 = 1:nrois-1
        tmproi1 = save_rois(iRoi1);
        tmpbin1 = zeros(nframes,1); tmpbin1(eventdata{tmproi1,7}) = 1;
        tmpx1 = roidata.roi_centroids_um(tmproi1,1);
        tmpy1 = roidata.roi_centroids_um(tmproi1,2);
        for iRoi2 = iRoi1+1:nrois
            % Distance
            tmproi2 = save_rois(iRoi2);
            tmpx2 = roidata.roi_centroids_um(tmproi2,1);
            tmpy2 = roidata.roi_centroids_um(tmproi2,2);
            dist = norm([tmpx1 tmpy1] -[tmpx2 tmpy2]);
            distance_mtrx(iRoi1,iRoi2) = dist;
            % Synchronicity expressed as Pearsons correlation coefficient
            tmpbin2 = zeros(nframes,1); tmpbin2(eventdata{tmproi2,7}) = 1;
            av1 = mean(tmpbin1); av2 = mean(tmpbin2);
            r = sum((av1-tmpbin1).*(av2-tmpbin2)) / sqrt(sum((av1-tmpbin1).^2)*sum((av2-tmpbin2).^2));
            synchro_mtrx(iRoi1,iRoi2) =  r;
            fill_mtrx(iRoi1,iRoi2) = false;
            roi_mtrx{iRoi1,iRoi2} = strcat(roinames{iRoi1},'-',roinames{iRoi2});
            synchro_dist_mtrx(cnt,1) = r; synchro_dist_mtrx(cnt,2) = dist;
            cnt = cnt+1;
        end
    end
    trsp_mtrx = distance_mtrx';
    distance_mtrx(fill_mtrx) = trsp_mtrx(fill_mtrx);
    trsp_mtrx = synchro_mtrx';
    synchro_mtrx(fill_mtrx) = trsp_mtrx(fill_mtrx);
    trsp_mtrx = roi_mtrx';
    roi_mtrx(fill_mtrx) = trsp_mtrx(fill_mtrx);
    synchro_dist_tbl = table();
    synchro_dist_tbl.Dist_um = synchro_dist_mtrx(:,2);
    synchro_dist_tbl.PearsonR = synchro_dist_mtrx(:,1);
    
    %% Calculate Synchronicity Histogram
    for iT = 1:nframes
        evid = find(isnan(syn_mtrx_Amp(:,iT))==0);
        synchro_tbl(iT).Time_s = timepoints(iT);
        synchro_tbl(iT).Eventcount = numel(evid);
        synchro_tbl(iT).Amp_Mean = mean(syn_mtrx_Amp(evid,iT));
        synchro_tbl(iT).Amp_FoF_Mean = mean(syn_mtrx_Amp_FoF(evid,iT));
        synchro_tbl(iT).Amp_dFoF_Mean  = mean(syn_mtrx_Amp_dFoF(evid,iT));
        synchro_tbl(iT).Amp_SD = std(syn_mtrx_Amp(evid,iT));
        synchro_tbl(iT).Amp_FoF_SD = std(syn_mtrx_Amp_FoF(evid,iT));
        synchro_tbl(iT).Amp_dFoF_SD = std(syn_mtrx_Amp_dFoF(evid,iT));
    end
    
    %% Recoring Summary
    rec_summary.Num_ROIs = nrois;
    rec_summary.ToT_Num_Events = totevents;
    rec_summary.Num_Events_Mean = totevents/nrois;
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
    rec_summary.REC_time_s = nframes.*ft;
    rec_summary.Area_um_Mean = mean([roi_summary(:).Area_um]);
    rec_summary.StructNorm_Area = mean([roi_summary(:).StructNorm_Area]);
    rec_summary.Dend_Area_um = roidata.dendritic_area_um;
    rec_summary.Bg_Area_um = roidata.background_area_um;
    rec_summary.FoV_Area_um = roidata.fov_area_um;
    rec_summary.Events_per_Dend_Area = totevents/roidata.dendritic_area_um;
    
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
    
end

% Background mask and ROI overview
imsize = size(roidata.aip);
magfac = round(34000/imsize(2))*500;
rim = zeros(imsize(1), imsize(2),1);
gim = zeros(imsize(1), imsize(2),1);
bim = zeros(imsize(1), imsize(2),1);
bg_col = [0.1529 0.1882 0.6706]; % RGB
rim(roidata.background_mask==1) = bg_col(1);
gim(roidata.background_mask==1) = bg_col(2);
bim(roidata.background_mask==1) = bg_col(3);
rgbim = cat(3, rim, gim,bim);
f1=figure('Position',figsize,'Visible','off'); imshow(rgbim, 'InitialMagnification', magfac);
saveas(f1,bg_png); close(f1);

if ~isempty(save_rois)
    f1=figure('Position',figsize,'Visible','off');
    imagesc(uint8(uint8_im));colormap('gray');axis('off');daspect([1 1 1]);
    for iRoi = 1:nrois
        B = roidata.roi_bounds{save_rois(iRoi),1};
        for k = 1:length(B)
            boundary = B{k};
            hold on, plot(boundary(:,2), boundary(:,1), 'Color',cmap(:,roidata.roi_bounds{save_rois(iRoi),2})', 'LineWidth', 1)
        end
    end
    saveas(f1,roi_ovw_png); close(f1);
    
    % Save cell with events as struct
    eventinfo = struct();
    eventinfo.eventdata = eventdata;
    save(ev_struct, 'eventinfo', '-v7.3');
end

disp('Files saved');
end