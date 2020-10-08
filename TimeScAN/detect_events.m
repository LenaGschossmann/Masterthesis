function [events, tmptrace] = detect_events(worktraces, seltrace)
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
% tmptrace: Vector of tmptrace data

%% Declare globally shared variables
global PLOTPEAKS RISETIME THRESHOLDTYPE THRESHOLD PEAKTHRESHOLD SMTHWIN FTIME...
    TRACEID TRACEDATA ROILIST RANGE TOTP ISBGSUB WHREVFIG evINFO BASEPOINTS PREFRAMES POSTFRAMES

%% Prepare revision display


%% Prepare variables
if strcmp(THRESHOLDTYPE, 'dF'), mode = 1; else, mode = 0; end
ptspre = ceil(0.5/FTIME); ptspost = ceil(1.5/FTIME);

if ptspre+ptspost+1 > diff(RANGE), RANGE = 1:TOTP; msgbox('The range was too small and set back!'); end

dftraces = diff(worktraces,1,1);
avtraces = mean(worktraces,1);
sdtraces = std(worktraces,1);
evcutouts = cell(1,numel(seltrace));

for iTr = 1:size(worktraces,2)
    tmptrace = dftraces(:,iTr);
    if mode %% Detect events based on delta F
        crossing = tmptrace > THRESHOLD;
        crossidx = find(crossing == 1);
        trcthreshold = THRESHOLD;
        trcpkthreshold = NaN;
        
    else %% Detect events based on sd
        % Define threshold for detection of event start and peak
        trcthreshold = avtraces(iTr)+THRESHOLD*sdtraces(iTr);
        trcpkthreshold = avtraces(iTr)+ PEAKTHRESHOLD*sdtraces(iTr);
        
        % Exclusion criterion 1 - SD value: if sd value too low, presumably smoothing was too
        % strong and likelihood of false positives increases a lot
        % if sd/av < 0.05
        %     aboveT = false(size(tmptrace));
        %     crossing = aboveT; peaks = aboveT;
        %     crossidx = NaN; peakidx = NaN; amps = NaN; ieis = NaN; eventrate = NaN; cviei = NaN; evkey = 'NaN';
        % else
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
        for iEv = 1:numev
            evcutouts{iEv, iTr} = tmptrace(crossidx(iEv)-ptspre:crossidx(iEv)-ptspost);
            peaks = false(size(tmptrace));                                          % Binarized vector indicating peaks
            peakidx = zeros(size(crossidx));                                        % Vector containing the indices of peak points
            amps = zeros(size(crossidx));                                           % Vector containing the event amplitudes (peak - average)
            delidx = [];
            riseThresh = ceil(RISETIME/FTIME);
            for iCross = 1:numel(crossidx)
                iPos = crossidx(iCross);
                %             test = false(size(tmptrace));
                %             while iPos <= numel(tmptrace) && aboveT(iPos) == 1
                %                 test(iPos) = true;
                %                 iPos = iPos+1;
                %             end
                add = 5;
                if crossidx(iCross)+add > numel(tmptrace), add= numel(tmptrace)-crossidx(iCross); end
                [valP,idxP] = max(tmptrace(crossidx(iCross):crossidx(iCross)+add));
                peakidx(iCross) = crossidx(iCross)+idxP-1;
                
                % Calculate delta F based on ROI average
                amps(iCross) = valP - mean(tmptrace(crossidx(iCross)-BASEPOINTS:crossidx(iCross)));
                
                %             % Exclusion criterion 2 - amplitude: If peak is too small, discard
                %             % event
                %             if strcmp(THRESHOLDTYPE, 'sd')
                %                 if valP <= trcpkthreshold || (peakidx(iCross)-crossidx(iCross)) > riseThresh
                %                     delidx = [delidx iCross];
                %                 end
                %             end
            end
            crossing(crossidx(delidx)) = false;
            crossidx(delidx) = []; peakidx(delidx) = []; amps(delidx) = [];
            peaks(peakidx) = true;
        end
        
        %% Calculate event-related parameters
        if PLOTPEAKS, ieis = diff(peakidx).*FTIME; evkey = 'peak';
        else, ieis = diff(crossidx).*FTIME; evkey = 'crossing';
        end
        cviei = std(ieis)/mean(ieis);
        %     eventrate = numel(crossidx)/(numel(tmptrace)*FTIME); % Eventrate calculated from number of events [Hz]
        eventrate = 1/median(ieis); % Eventrate calculated from intereventintervals [Hz]
        % end
    else
        crossidx = NaN;
        peaks = zeros(size(crossing));
        peakidx = NaN;
        amps = NaN;
        ieis = NaN;
        cviei = NaN;
        eventrate = NaN;
    end
    
    
    % Save information in struct
    x = seltrace(iTr);
    evINFO(x).rangein = RANGE(1);
    evINFO(x).rangeout = RANGE(2);
    evINFO(x).ftime = FTIME;
    evINFO(x).smthwin = SMTHWIN;
    evINFO(x).basepoints = BASEPOINTS;
    evINFO(x).preframes = PREFRAMES;
    evINFO(x).postframes = POSTFRAMES;
    evINFO(x).average = avtraces(iTr);
    evINFO(x).sd = sdtraces(iTr);
    evINFO(x).threshold = trcthreshold;
    evINFO(x).peakthreshold = trcpkthreshold;
    evINFO(x).binarycross = crossing;
    evINFO(x).binarypeaks = peaks;
    evINFO(x).crossidx = crossidx;
    evINFO(x).peakidx = peakidx;
    evINFO(x).amps = amps;
    evINFO(x).avamp = mean(amps);
    evINFO(x).eventtype = THRESHOLDTYPE;
    evINFO(x).ieis = ieis;
    evINFO(x).aviei = mean(ieis);
    evINFO(x).cviei = cviei;
    evINFO(x).eventratetype = 'IEI based';
    evINFO(x).eventrate = eventrate;
end
end