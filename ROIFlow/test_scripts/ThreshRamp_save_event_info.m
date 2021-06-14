function results  = ThreshRamp_save_event_info(eventdata, roidata, path, filename)

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
all_nev = 0;
results = 0;

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
        all_nev = all_nev + nevents;
        
        %% Traces
        traces.(roinames{iRoi}) = roidata.traces(tmproi,:)';
        dFoF_traces.(roinames{iRoi}) =  roidata.dFoF_traces(tmproi,:)';
        FoF_traces.(roinames{iRoi}) =  roidata.FoF_traces(tmproi,:)';
        
    end
    results = all_nev/nrois;
%     results = all_nev;
end

end