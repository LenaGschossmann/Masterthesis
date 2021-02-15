function [roivals, roidFoF, yrange, avvals, alldFoF] = calc_roi_av_trace(roiidx, vals, pr)

% Declare globally shared variables
global roiINFO

roimask = roiINFO(roiidx).mask;
roipr = roiINFO(roiidx).plotrange;
% Calculate whole trace average
collapsey = sum(roimask,1);
xmin = find(collapsey,1,'first'); xmax = find(collapsey,1,'last');
avvals = mean(vals(:,xmin:xmax),2);
av = mean(avvals);

skip = 0;
if any(roipr ~= pr)
    if pr(1) <= roipr(1) && pr(2) <= roipr(1), yrange = 1;
    elseif pr(1) <= roipr(1) && pr(2) <= roipr(2), ycut = 1:pr(2)-roipr(1)+1; ypaste = roipr(1)-pr(1)+1:pr(2)-pr(1)+1;
    elseif pr(1) <= roipr(1) && pr(2) >= roipr(2), ycut = 1:roipr(2)-roipr(1)+1; ypaste = roipr(1)-pr(1)+1:roipr(2)-pr(1)+1;
    elseif pr(1) >= roipr(1) && pr(2) <= roipr(2), ycut = pr(1)-roipr(1)+1:pr(2)-roipr(1)+1; ypaste = pr(1)-pr(1)+1:pr(2)-pr(1)+1;
    elseif pr(1) >= roipr(1) && pr(2) >= roipr(2), ycut = pr(1)-roipr(1)+1:roipr(2)-roipr(1)+1; ypaste = pr(1)-pr(1)+1:roipr(2)-pr(1)+1;
    end
    newmask = false(diff(pr)+1, size(vals,2));
    newmask(ypaste,:) = roimask(ycut,:);
else
    newmask = roimask;
end

if ~skip
    % Calculate ROI average
    collapsex = sum(newmask,2);
    ymin = find(collapsex,1,'first'); ymax = find(collapsex,1,'last');
    % Av values
    roivals = avvals(pr(1)+ymin-1:pr(1)+ymax-1,1);
    % dFoF
    [~, roidFoF] = rollBase_dFoF(roivals);
    [~, alldFoF] = rollBase_dFoF(avvals);
    yrange = [ymin ymax];
else
    roivals = []; roidFoF = []; alldFoF = []; yrange = [];
end
end