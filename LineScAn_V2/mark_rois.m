function tmpmtrx = mark_rois(tmpmtrx, pr, rois)

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
    if roiINFO(roiidx(iRoi)).mode == 1
        newmask = false(size(tmpmtrx));
        x1 = find(roimask(1,:),1,'first'); x2 = find(roimask(1,:),1,'last');
        newmask(:,x1:x2) = true;
    else
        % Check if ROI still lies in plotted range
        roipr = roiINFO(roiidx(iRoi)).plotrange;
        if any(roipr ~= pr)
            if pr(1) <= roipr(1) && pr(2) <= roipr(1), newmask = false(size(roimask));
            elseif pr(1) <= roipr(1)&& pr(2) <= roipr(2), ycut = 1:pr(2)-roipr(1)+1; ypaste = roipr(1)-pr(1)+1:pr(2)-pr(1)+1;
            elseif pr(1) <= roipr(1) && pr(2) >= roipr(2), ycut = 1:roipr(2)-roipr(1)+1; ypaste = roipr(1)-pr(1)+1:roipr(2)-pr(1)+1;
            elseif pr(1) >= roipr(1) && pr(2) <= roipr(2), ycut = pr(1)-roipr(1)+1:pr(2)-roipr(1)+1; ypaste = pr(1)-pr(1)+1:pr(2)-pr(1)+1;
            elseif pr(1) >= roipr(1) && pr(2) >= roipr(2), ycut = pr(1)-roipr(1)+1:roipr(2)-roipr(1)+1; ypaste = pr(1)-pr(1)+1:roipr(2)-pr(1)+1;
            end
            newmask = false(size(tmpmtrx));
            newmask(ypaste,:) = roimask(ycut,:);
        else
            newmask = roimask;
        end
    end
    tmpmtrx(newmask) = tmpmtrx(newmask)+tmpmarkval;
end
end