function [traceidx] = create_trc(roiidx, iRoi, traceidx)

% Declare globally shared variables
global SMTHWIN roiINFO traceINFO COMPOSITE2D WINSZ PLOTRANGE

if isempty(traceINFO(1).roiID), tmpidx = 1; else, tmpidx = size(traceINFO,2)+1; end
traceINFO(tmpidx).roiID = roiINFO(roiidx(iRoi)).ID;
traceINFO(tmpidx).params{2,1} = WINSZ;
traceINFO(tmpidx).params{3,1} = SMTHWIN;
traceINFO(tmpidx).save = true;
traceINFO(tmpidx).currmode = '  Average';
traceINFO(tmpidx).showtot = 0;

% Call Analysis functions
disp('Bin px traces...');
averaged = average_linescan(COMPOSITE2D, WINSZ);
markedroi = averaged(PLOTRANGE(1):PLOTRANGE(2),:);
scmarkedroi = mark_rois(markedroi, roiidx(iRoi));
%         scmarkedroi = scale_data(markedroi, climraw);
disp('Smooth trace...');
smoothed = smooth_data(averaged, SMTHWIN);
[allvals, alldFoF, allFoF, alldF] = calc_roi_av_trace(roiidx(iRoi), smoothed);
traceidx = [traceidx tmpidx];
traceINFO(tmpidx).roi_av = {allvals};
traceINFO(tmpidx).dFoF_roi_av = {alldFoF};
traceINFO(tmpidx).FoF_roi_av = {allFoF};
traceINFO(tmpidx).dF_roi_av = {alldF};
traceINFO(tmpidx).plotmarked = {scmarkedroi};
end