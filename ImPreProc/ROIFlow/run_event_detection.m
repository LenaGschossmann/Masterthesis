function [eventdata, selrois, params] = run_event_detection(eventdata,roilist, dFoFtraces, FoFtraces, traces,  params)

n_frames = size(dFoFtraces,2);
peakwin = params.peak_win;
threshold_fac = params.ev_fac;
save_fac = params.ev_save_fac;
ui = inputdlg({'Threshold Factor:', 'Save Factor', 'Specific ROIs (sep: ,)'} ,'Set Values',[1 30], {num2str(threshold_fac), num2str(save_fac), ''});
if ~isempty(ui)
    threshold_fac = str2num(ui{1}); params.ev_fac = threshold_fac;
    save_fac = str2num(ui{2}); params.ev_save_fac = save_fac;
end
if isempty(ui), selrois = 1:size(dFoFtraces,1);
elseif isempty(ui{3}), selrois = 1:size(dFoFtraces,1);
else
    tmprois = split(ui{3},',');
    selrois = NaN(numel(tmprois),1);
    for iRoi = 1:numel(tmprois),selrois(iRoi)=roilist{str2num(tmprois{iRoi}),2}; end
end
nrois = numel(selrois);
rand_pts = round(n_frames*0.15);

for iRoi = 1:nrois
    tmproi = selrois(iRoi);
    tmpdFoF = dFoFtraces(tmproi,:);
    tmpFoF = FoFtraces(tmproi,:);
    tmptrc = traces(tmproi,:);
    
    %% Define threshold based on sd of dFoF values
    %     perc_thresh_val = prctile(dFoFim,perc_thresh,3);
    sd_thresh_val = std(tmpdFoF);
    onsetbinary = tmpdFoF > (sd_thresh_val*threshold_fac);
    onset_idx = find(onsetbinary == 1);
    peakbinary = false(1,n_frames); % Binarized vector indicating peaks
    
    %% Determine peak after threshold onset
    if numel(onset_idx) > 0
        numev = numel(onset_idx);
        peakidx = zeros(size(onset_idx,2)); % Vector containing the indices of peak points
        dFoF_amps = zeros(size(onset_idx,2),1);
        FoF_amps = zeros(size(onset_idx,2),1);
        amps = zeros(size(onset_idx,2),1);
        
        for iEv = 1:numev
            if onset_idx(iEv)+peakwin > n_frames, peakwin= n_frames-onset_idx(iEv); end
            [~,idxP] = max(tmpdFoF(onset_idx(iEv):onset_idx(iEv)+peakwin));
            peakidx(iEv) = onset_idx(iEv)+idxP-1;
            peakbinary(peakidx(iEv)) = true;
            
            dFoF_amps(iEv,1) =  tmpdFoF(peakidx(iEv));
            FoF_amps(iEv,1) = tmpFoF(peakidx(iEv));
            amps(iEv,1) = tmptrc(peakidx(iEv));
        end
        
        %% Mark events for revision
        % Use dFoF threshold value to mark for revision
        critical = tmpdFoF(onset_idx) < (sd_thresh_val*save_fac);
        revise_idx = onset_idx(critical);
        
        % Draw random points to quantify noise
        possiblepos = 1:numel(tmpdFoF)-peakwin;
        randpos = randi(numel(possiblepos), rand_pts, 1);
        randpos = possiblepos(randpos);
        randmax = zeros(rand_pts,1);
        for iP = 1:rand_pts % Calculate amplitude of maxima in the specified window after event onset with randomly chosen points
            randmax(iP) = max(tmpdFoF(randpos(iP):randpos(iP)+peakwin));
        end
        
        neg_dFoF_amps = abs(tmpdFoF(tmpdFoF < 0));
        
        save_idx = onset_idx;
        for iRev = 1:numel(revise_idx), save_idx(revise_idx(iRev)==save_idx) = []; end
        
        eventdata{tmproi,1} = onsetbinary';
        eventdata{tmproi,2} = peakbinary';
        eventdata{tmproi,3} = onset_idx;
        eventdata{tmproi,4} = peakidx;
        eventdata{tmproi,5} = [dFoF_amps FoF_amps amps];
        eventdata{tmproi,6} = revise_idx;
        eventdata{tmproi,7} = save_idx;
        eventdata{tmproi,8} = randmax;
        eventdata{tmproi,9} = neg_dFoF_amps';
        eventdata{tmproi,10} = [threshold_fac save_fac];
    else
        eventdata{tmproi,1} = onsetbinary';
        eventdata{tmproi,2} = peakbinary';
        eventdata{tmproi,3} = [];
        eventdata{tmproi,4} = [];
        eventdata{tmproi,5} = [];
        eventdata{tmproi,6} = [];
        eventdata{tmproi,7} = [];
        eventdata{tmproi,8} = [];
        eventdata{tmproi,9} = [];
        eventdata{tmproi,10} = [threshold_fac save_fac];
    end
end

end