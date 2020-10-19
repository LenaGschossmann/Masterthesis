function [events, tmpdftrace] = detect_events(worktraces, seltrace)
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
% tmpdftrace: Vector of tmpdftrace data

%% Declare globally shared variables
global PLOTPEAKS THRESHOLDTYPE THRESHOLD PEAKTHRESHOLD SMTHWIN FTIME...
    RANGE TOTP evINFO BASEPOINTS PREPOINTS POSTPOINTS

%% Initialize variables
if strcmp(THRESHOLDTYPE, 'dF'), mode = 1; else, mode = 0; end
ptspre = PREPOINTS; ptspost = POSTPOINTS;
if ptspre+ptspost+1 > diff(RANGE), RANGE = 1:TOTP; msgbox('The range was too small and set back!'); end

peakwindow = 5; % Frames after threshold crossing in which peak is detected
dftraces = diff(worktraces,1,1);
avtraces = mean(worktraces,1);
sdtraces = std(worktraces,1);

for iTr = 1:size(worktraces,2)
    tmpdftrace = dftraces(:,iTr);
    tmptrace = worktraces(:,iTr);
    if mode %% Detect events based on delta F
        crossing = tmpdftrace > THRESHOLD;
        crossidx = find(crossing == 1);
        ii = 1;
        while ii <= numel(crossidx)-1
            if crossidx(ii+1)-crossidx(ii) < 5
                crossing(crossidx(ii)) = 0;
                crossidx(ii+1) = [];
            else
                ii = ii+1;
            end
        end
        trcthreshold = THRESHOLD;
        trcpkthreshold = NaN;
        crossing(end+1) = 0;
        
    else %% Detect events based on sd
        % Define threshold for detection of event start and peak
        trcthreshold = avtraces(iTr)+THRESHOLD*sdtraces(iTr);
        trcpkthreshold = avtraces(iTr)+ PEAKTHRESHOLD*sdtraces(iTr);
        
        aboveT = tmptrace >= trcthreshold;                                  % Binarized vector indicating all points above threshold
        crossing = false(size(aboveT));                                     % Binarized vector indicating points of threshold crossing
        % Determining points of threshold crossing (only if 2 consecutive points above threshold)
        iPos = 2;
        while iPos < numel(tmptrace)-2
            if aboveT(iPos-1) == 0 && all(aboveT(iPos:iPos+2)) == 1,crossing(iPos) = true; end
            iPos = iPos+1;
        end
        crossidx = find(crossing == 1);                                     % Vector containing the indices of threshold crossing points
    end
    
    %% Determine peak after threshold crossing
    if numel(crossidx) > 0
        numev = numel(crossidx);
        peaks = false(size(tmptrace));                                      % Binarized vector indicating peaks
        peakidx = zeros(size(crossidx));                                    % Vector containing the indices of peak points
        amps = zeros(size(crossidx));
        dFoF_amps = zeros(size(crossidx));
        baselinevals = zeros(size(crossidx));
        risetime = zeros(size(crossidx));
        decaytime = zeros(size(crossidx));
        for iEv = 1:numev
            if crossidx(iEv)+peakwindow > numel(tmptrace), peakwindow= numel(tmptrace)-crossidx(iEv); end
            [valP,idxP] = max(tmptrace(crossidx(iEv):crossidx(iEv)+peakwindow));
            if strcmp(THRESHOLDTYPE,'dF') || (strcmp(THRESHOLDTYPE,'sd') && valP >= avtraces(iTr)+PEAKTHRESHOLD*sdtraces(iTr))
                peakidx(iEv) = crossidx(iEv)+idxP-1;
                peaks(peakidx(iEv)) = true;
                
                % Calculate baseline value
                if crossidx(iEv)-BASEPOINTS < 1, startval = 1; else, startval = crossidx(iEv)-BASEPOINTS; end
                baselinevals(iEv) = mean(tmptrace(startval:crossidx(iEv)));
                
                amps(iEv) =  valP-baselinevals(iEv);
                dFoF_amps(iEv) = amps(iEv)/baselinevals(iEv);
            end
            risetime(iEv) = (peakidx(iEv)-crossidx(iEv))*FTIME;
            % Decay
            iP = peakidx(iEv);
            while decaytime(iEv) == 0 && iP <= size(tmptrace,iTr)
                if tmptrace(iP) <= baselinevals(iEv), decaytime(iEv) = (iP-peakidx(iEv))*FTIME; end
                iP = iP+1;
            end
        end
        
        %% Calculate event-related parameters
        if PLOTPEAKS, ieis = diff(peakidx).*FTIME;
        else, ieis = diff(crossidx).*FTIME;
        end
        cviei = std(ieis)/mean(ieis);
        %     eventrate = numel(crossidx)/(numel(tmpdftrace)*FTIME); % Eventrate calculated from number of events [Hz]
        eventrate = 1/median(ieis); % Eventrate calculated from intereventintervals [Hz]
        
    else
        crossidx = NaN;
        peaks = zeros(size(crossing));
        peakidx = NaN;
        amps = NaN;
        dFoF_amps = NaN;
        ieis = NaN;
        cviei = NaN;
        eventrate = NaN;
        baselinevals = NaN;
        risetime = NaN;
        decaytime = NaN;
    end
    
    % Save information in struct
    x = seltrace(iTr);
    evINFO(x).rangein = RANGE(1);
    evINFO(x).rangeout = RANGE(2);
    evINFO(x).ftime = FTIME;
    evINFO(x).smthwin = SMTHWIN;
    evINFO(x).baseframes = BASEPOINTS;
    evINFO(x).baselinevalues = baselinevals;
    evINFO(x).average = avtraces(iTr);
    evINFO(x).sd = sdtraces(iTr);
    evINFO(x).thresholdtype = THRESHOLDTYPE; 
    evINFO(x).threshold = trcthreshold;
    evINFO(x).peakthreshold = trcpkthreshold;
    evINFO(x).binarycross = crossing;
    evINFO(x).binarypeaks = peaks;
    evINFO(x).crossidx = crossidx;
    evINFO(x).peakidx = peakidx;
    evINFO(x).amps = amps;
    evINFO(x).amps_dFoF = dFoF_amps;
    evINFO(x).avamp = mean(amps);
    evINFO(x).avamp_dFoF = mean(dFoF_amps);
    evINFO(x).ieis = ieis;
    evINFO(x).aviei = mean(ieis);
    evINFO(x).cviei = cviei;
    evINFO(x).eventratetype = 'IEI based';
    evINFO(x).eventrate = eventrate;
    evINFO(x).risetimes = risetime;
    evINFO(x).decaytimes = decaytime;
end
end