function [traceidx, roiidx] = access_trc()

global roiINFO traceINFO SMTHWIN WINSZ

% ROI Parameters
[~,roiidx] = find([roiINFO(:).selected] == 1);
numrois = numel(roiidx);

traceidx = [];
for iRoi = 1:numrois
    trcexists = false;
    existbin = [traceINFO(:).roiID] == roiINFO(roiidx(iRoi)).ID;
    if any(existbin)
        existidx = find(existbin);
        iEx = 1;
        while iEx <= numel(existidx)
            traceINFO(existidx(iEx)).save = 0;
            if traceINFO(existidx(iEx)).params{2,1} == WINSZ &&...
                    traceINFO(existidx(iEx)).params{3,1} == SMTHWIN
                trcexists = true;
                traceidx = [traceidx existidx(iEx)];
                break;
            end
            iEx = iEx+1;
        end
    end
    if ~trcexists
        disp('Create new trace...');
        traceidx = create_trc(roiidx, iRoi, traceidx);
    end
end
end