function tmpmtrx = mark_rois(tmpmtrx, rois)

% Declare globally shared variables
global roiINFO

if strcmp(rois,'all')
    roiidx = 1:size(roiINFO,2);
    if isempty(roiINFO(end).name), roiidx= roiidx(1:end-1); end
else
    roiidx = rois;
end
tmpmarkval = max(tmpmtrx,[],'all')*0.3;

for iRoi = 1:numel(roiidx)
    roimask = roiINFO(roiidx(iRoi)).mask;
    newmask = false(size(tmpmtrx));
    x1 = find(roimask(1,:),1,'first'); x2 = find(roimask(1,:),1,'last');
    newmask(:,x1:x2) = true;
    
    tmpmtrx(newmask) = tmpmtrx(newmask)+tmpmarkval;
end
end