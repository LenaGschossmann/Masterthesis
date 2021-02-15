function [traceidx] = create_trc(figidx, roiidx, iRoi, traceidx, pr, wsz)

% Declare globally shared variables
global THRESHOLD PEAKTHRESHOLD SMTHWIN figINFO roiINFO traceINFO IMHW FTIME COMPOSITE2D

if isempty(traceINFO(1).figID), tmpidx = 1; else, tmpidx = size(traceINFO,2)+1; end
traceINFO(tmpidx).figID = figINFO(figidx).IDs;
traceINFO(tmpidx).roiID = roiINFO(roiidx(iRoi)).ID;
traceINFO(tmpidx).fig_params{1,1} = pr;
traceINFO(tmpidx).fig_params{2,1} = wsz;
traceINFO(tmpidx).fig_params{3,1} = SMTHWIN;
traceINFO(tmpidx).fig_params{4,1} = THRESHOLD;
traceINFO(tmpidx).fig_params{4,1} = PEAKTHRESHOLD;
traceINFO(tmpidx).save = true;
traceINFO(tmpidx).currmode = '  Average';
traceINFO(tmpidx).showtot = 0;
% Call Analysis functions
averaged = average_linescan(COMPOSITE2D, wsz);
markedroi = averaged(pr(1):pr(2),:);
scmarkedroi = mark_rois(markedroi,pr, roiidx(iRoi));
%         scmarkedroi = scale_data(markedroi, climraw);
[tmpvals,tmpdFoF, yrange, allvals, alldFoF] = calc_roi_av_trace(roiidx(iRoi), averaged, pr);
if isempty(yrange)
    traceINFO(tmpidx).binned_roi_av = {NaN};
    traceINFO(tmpidx).dFoF_roi_av = {NaN};
    traceINFO(tmpidx).timestamp = {NaN};
    traceINFO(tmpidx).plotmarked = {NaN};
    msgbox(strcat('ROI # ', num2str(iRoi), ' lies outside the plotted range!'));
else
    traceidx = [traceidx tmpidx];
    traceINFO(tmpidx).binned_roi_av = {tmpvals};
    traceINFO(tmpidx).dFoF_roi_av = {tmpdFoF};
    traceINFO(tmpidx).plotmarked = {scmarkedroi};
    traceINFO(tmpidx).timestamp = {((pr(1)+yrange(1)-1)*FTIME:FTIME:(pr(1)+yrange(2)-1)*FTIME)'};
    % Extract trace parameters
    [events, smoothed] = get_trc_params(allvals, [pr(1)+yrange(1)-1 pr(1)+yrange(2)-1], [], []);
    traceINFO(tmpidx).events = events;
    traceINFO(tmpidx).smoothed = {smoothed};
end
if roiINFO(roiidx(iRoi)).mode == 1
    traceINFO(tmpidx).tot_binned_roi_av = {allvals};
    traceINFO(tmpidx).tot_dFoF_roi_av = {alldFoF};
    traceINFO(tmpidx).tot_timestamp = {(1*FTIME:FTIME:IMHW(1)*FTIME)'};
    % Extract trace parameters
    [~, smoothed] = get_trc_params(allvals,[],[],[]);
%     traceINFO(tmpidx).tot_events = events;
    traceINFO(tmpidx).tot_smoothed = {smoothed};
end
end