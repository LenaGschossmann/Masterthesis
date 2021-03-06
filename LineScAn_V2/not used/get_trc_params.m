function [events, smoothed] = get_trc_params(alldata, range, av, sd)
% Obtain event-related parameters of spike-train data
% Input:
% alldata: Vector of data representing neural activity
% range: [start end] Measurement points taken into account for extraction
% of event-related parameters(if empty, all available points are used); for
% calculation of average and sd, always all timepoints are used
% av: Average of input trace (if not yet calculated, leave empty ([]))
% sd: SD of input trace (if not yet calculated, leave empty ([]))
% Output:
% events: (struct) Information regarding extracted parameters
% smoothed: Vector of smoothed data

% Declare globally shared variables
global PLOTPEAKS RISETIME THRESHOLD PEAKTHRESHOLD SMTHWIN FTIME

% Define range of data from which parameters are obtained
if isempty(range), range = 1:numel(alldata); else, range = range(1):range(2); end
% If average and sd are not given, calculate them based on data of all timepoints
if isempty(av) || isempty(sd)
    % Smooth data
    smthdata = smooth_data(alldata, SMTHWIN);
    smoothed = smthdata(range);
    % Threshold crossing: smoothed data and 1st derivative
    av = mean(smthdata);
    sd = std(smthdata);
else
    smoothed = smooth_data(alldata(range), SMTHWIN);
end

% Define threshold for detection of event start and peak
trcthreshold = av+THRESHOLD*sd;
trcpkthreshold = av + PEAKTHRESHOLD*sd;

aboveT = smoothed >= trcthreshold;                                      % Binarized vector indicating all points above threshold
crossing = false(size(aboveT));                                         % Binarized vector indicating points of threshold crossing

% Determining points of threshold crossing (only if 2 consecutive
% points above threshold)
iPos = 2;
while iPos < numel(smoothed)-2
    if aboveT(iPos-1) == 0 && all(aboveT(iPos:iPos+2)) == 1,crossing(iPos) = true; end
    iPos = iPos+1;
end
crossidx = find(crossing == 1);                                         % Vector containing the indices of threshold crossing points

% Determine peak after threshold crossing
peaks = false(size(smoothed));                                          % Binarized vector indicating peaks
peakidx = zeros(size(crossidx));                                        % Vector containing the indices of peak points
amps = zeros(size(crossidx));                                           % Vector containing the event amplitudes (peak - average)
delidx = [];
riseThresh = ceil(RISETIME/FTIME);
for iCross = 1:numel(crossidx)
    iPos = crossidx(iCross);
    test = false(size(aboveT));
    while iPos <= numel(aboveT) && aboveT(iPos) == 1
        test(iPos) = true;
        iPos = iPos+1;
    end
    [valP,idxP] = max(smoothed(test));
    peakidx(iCross) = crossidx(iCross)+idxP-1;
    
    % Calculate delta F based on ROI average
    amps(iCross) = valP - av;
    
    % Exclusion criterion 2 - amplitude: If peak is too small, discard
    % event
    if valP <= trcpkthreshold || (peakidx(iCross)-crossidx(iCross)) > riseThresh
        delidx = [delidx iCross];
    end
end
crossing(crossidx(delidx)) = false;
crossidx(delidx) = []; peakidx(delidx) = []; amps(delidx) = [];
peaks(peakidx) = true;

% Calculate event-related parameters
if PLOTPEAKS, ieis = diff(peakidx).*FTIME; evkey = 'peak';
else, ieis = diff(crossidx).*FTIME; evkey = 'crossing';
end
cviei = std(ieis)/mean(ieis);
%     eventrate = numel(crossidx)/(numel(smoothed)*FTIME); % Eventrate calculated from number of events [Hz]
eventrate = 1/median(ieis); % Eventrate calculated from intereventintervals [Hz]
% end

% Save information in struct
events=struct();
events.average = av;
events.sd = sd;
events.threshold = trcthreshold;
events.peakthreshold = trcpkthreshold;
events.eventtype = evkey;
events.eventrate = eventrate;
events.aviei = mean(ieis);
events.avamp = mean(amps);
events.cviei = cviei;
events.suprathreshold = aboveT;
events.crossings = crossing;
events.peaks = peaks;
events.crossidx = crossidx;
events.peakidx = peakidx;
events.amps = amps;
events.ieis = ieis;
end