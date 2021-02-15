function template = get_crosscorr_template(evdata, dFoFtraces, peakwin)
%% Define template for auto modus
% Find largest events
ev_list = [];
nrois = size(evdata,1);
temp_n = 10;

for iRoi = 1:nrois
    nev = numel(evdata{iRoi,3});
    ev_list = [ev_list; repelem(iRoi,nev)' evdata{iRoi,3}' evdata{iRoi,5}(:,1)];
end

[~,idx] = sort(ev_list(:,3), 'descend');

if numel(idx) <= temp_n, temp_n = numel(idx); end 

template = [];
for iN = 1:temp_n
    roi_id = ev_list(idx(iN),1);
    onset = ev_list(idx(iN),2);
    template = [template; dFoFtraces(roi_id,onset:onset+peakwin)];
end
template = median(template,1);

end