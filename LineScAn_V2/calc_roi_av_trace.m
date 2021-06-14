function [avvals, alldFoF, allFoF, alldF] = calc_roi_av_trace(roiidx, vals)

% Declare globally shared variables
global roiINFO FTIME dFWIN dFSHIFT

deltawinsz = round(dFWIN/FTIME);
roimask = roiINFO(roiidx).mask;

% Calculate whole trace average
collapsey = sum(roimask,1);
xmin = find(collapsey,1,'first'); xmax = find(collapsey,1,'last');
avvals = mean(vals(:,xmin:xmax),2);

% Calculate ROI average
disp('Calculate dF/F, FoF, dF trace...');
lin_avvals = avvals';
[alldF, allFoF, alldFoF] = rollBase_dFoF(lin_avvals,deltawinsz,dFSHIFT, 'roll');
allFoF = allFoF'; alldFoF = alldFoF'; alldF = alldF';

end