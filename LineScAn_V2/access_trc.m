function [traceidx, roiidx, figidx] = access_trc()

global figINFO roiINFO traceINFO FTIME SMTHWIN

% Image Processing Parameters
figures = get(groot,'Children');
findfig = strncmp({figures(:).Name}, 'Fig',3);
currfig = figures(find(findfig,1, 'first'));
if ~isempty(currfig)
    figid = currfig.UserData;
else
    warning('As no figure is open, the parameters of the last created figure are used');
    figid = max([figINFO(:).IDs]);
end
[~,figidx] = find([figINFO(:).IDs] == figid);

% Extract parameters of respective figure window
pr = figINFO(figidx).plotrange;
wsz = figINFO(figidx).avwinsize;

% ROI Parameters
[~,roiidx] = find([roiINFO(:).selected] == 1);
numrois = numel(roiidx);

traceidx = [];
for iRoi = 1:numrois
    trcexists = false;
    existidx1 = [traceINFO(:).figID] == figid;
    existidx2 = [traceINFO(:).roiID] == roiINFO(roiidx(iRoi)).ID;
    if any(existidx1 & existidx2)
        existidx = find(existidx1 & existidx2);
        iEx = 1;
        while iEx <= numel(existidx)
            traceINFO(existidx(iEx)).save = 0;
            if all(traceINFO(existidx(iEx)).fig_params{1,1} == pr) &&...
                    traceINFO(existidx(iEx)).fig_params{2,1} == wsz &&...
                    traceINFO(existidx(iEx)).fig_params{3,1} == SMTHWIN
                trcexists = true;
                traceidx = [traceidx existidx(iEx)];
                break;
            end
            iEx = iEx+1;
        end
    end
    if ~trcexists
        disp('Create new trace...');
        traceidx = create_trc(figidx, roiidx, iRoi, traceidx, pr, wsz);
    end
end
end