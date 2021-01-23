function [events, tmpdftrace] = detect_events(worktraces, seltrace)
% Obtain event-related parameters of spike-train data
% Input:
% Output:
% events: (struct) Information regarding extracted parameters
% tmpdftrace: Vector of tmpdftrace data

%% Declare globally shared variables
global PLOTPEAKS THRESHOLDTYPE THRESHOLD PEAKTHRESHOLD SMTHWIN FTIME...
    RANGE TOTP evINFO BASEPOINTS PREPOINTS POSTPOINTS PEAKWIN RANDNOISEPTS...
    RANDMAXTHRESH

%% Initialize variables
if strcmp(THRESHOLDTYPE, 'dF'), mode = 1; else, mode = 0; end
ptspre = PREPOINTS; ptspost = POSTPOINTS;
if ptspre+ptspost+1 > diff(RANGE), RANGE = 1:TOTP; msgbox('The range was too small and set back!'); end

dftraces = diff(worktraces,1,1);
avtraces = mean(worktraces,1);
sdtraces = std(worktraces,1);

for iTr = 1:size(worktraces,2)
    tmpdftrace = dftraces(:,iTr);
    tmptrace = worktraces(:,iTr);
    if mode %% Detect events based on delta F
        onset = tmpdftrace > THRESHOLD;
        onsetidx = find(onset == 1);
        trcthreshold = THRESHOLD;
        trcpkthreshold = NaN;
        onset(end+1) = 0;
        
    else %% Detect events based on sd
        % Define threshold for detection of event start and peak
        trcthreshold = avtraces(iTr)+THRESHOLD*sdtraces(iTr);
        trcpkthreshold = avtraces(iTr)+ PEAKTHRESHOLD*sdtraces(iTr);
        
        aboveT = tmptrace >= trcthreshold;                                  % Binarized vector indicating all points above threshold
        onset = false(size(aboveT));                                     % Binarized vector indicating points of threshold onset
        % Determining points of threshold onset (only if 2 consecutive points above threshold)
        iPos = 2;
        while iPos < numel(tmptrace)-2
            if aboveT(iPos-1) == 0 && all(aboveT(iPos:iPos+2)) == 1,onset(iPos) = true; end
            iPos = iPos+1;
        end
        onsetidx = find(onset == 1);                                     % Vector containing the indices of threshold onset points
    end
    
    %% Determine peak after threshold onset
    if numel(onsetidx) > 0
        numev = numel(onsetidx);
        peaks = false(size(tmptrace));                                      % Binarized vector indicating peaks
        peakidx = zeros(size(onsetidx));                                    % Vector containing the indices of peak points
        amps = zeros(size(onsetidx));
        dFoF_amps = zeros(size(onsetidx));
        baselinevals = zeros(size(onsetidx));
        risetime = zeros(size(onsetidx));
        decaytime = zeros(size(onsetidx));
        for iEv = 1:numev
            if onsetidx(iEv)+PEAKWIN > numel(tmptrace), PEAKWIN= numel(tmptrace)-onsetidx(iEv); end
            [valP,idxP] = max(tmptrace(onsetidx(iEv):onsetidx(iEv)+PEAKWIN));
            if strcmp(THRESHOLDTYPE,'dF') || (strcmp(THRESHOLDTYPE,'sd') && valP >= avtraces(iTr)+PEAKTHRESHOLD*sdtraces(iTr))
                peakidx(iEv) = onsetidx(iEv)+idxP-1;
                peaks(peakidx(iEv)) = true;
                
                % Calculate baseline value
                if onsetidx(iEv)-BASEPOINTS < 1, startval = 1; else, startval = onsetidx(iEv)-BASEPOINTS; end
                baselinevals(iEv) = mean(tmptrace(startval:onsetidx(iEv)));
                
                amps(iEv) =  valP-baselinevals(iEv);
                dFoF_amps(iEv) = amps(iEv)/baselinevals(iEv);
            end
            risetime(iEv) = (peakidx(iEv)-onsetidx(iEv))*FTIME;
            % Decay
            iP = peakidx(iEv);
            while decaytime(iEv) == 0 && iP <= size(tmptrace,iTr)
                if tmptrace(iP) <= baselinevals(iEv), decaytime(iEv) = (iP-peakidx(iEv))*FTIME; end
                iP = iP+1;
            end
        end
        
        %% Discard event if amplitude lies within the random noise amplitude distribution
        possiblepos = 1:numel(tmptrace)-PEAKWIN;
        delpos = [];
        for iEv = 1:numev, delpos = [delpos peakidx(iEv)-PEAKWIN:peakidx(iEv)]; end
        possiblepos(delpos) = [];
        randpos = randi(numel(possiblepos), RANDNOISEPTS, 1);
        randpos = possiblepos(randpos);
        randmax = zeros(RANDNOISEPTS,1);
        for iP = 1:RANDNOISEPTS % Calculate amplitude of maxima in the specified window after event onset with randomly chosen points
            tmpmax = max(tmptrace(randpos(iP):randpos(iP)+PEAKWIN));
            randamp(iP) = tmpmax - get_base(tmptrace, randpos(iP));
        end
        randthresh = RANDMAXTHRESH*max(randamp);
        del = [];
        for iEv = 1:numev
            if amps(iEv) <= randthresh, del = [del iEv]; end % Exclude an event if its amplitude lies within the boundaries of noise amplitude
        end
        if ~isempty(del)
            onsetidx(del) = []; peakidx(del) = [];
            amps(del) = []; dFoF_amps(del) = []; baselinevals(del) = [];
%             risetime(del) = []; decaytime(del) = [];
            onset(del) = false; peaks(del) = false;
        end
    end
    
    %% Calculate event-related parameters
    if numel(onsetidx) > 1
        if PLOTPEAKS, ieis = diff(peakidx).*FTIME;
        else, ieis = diff(onsetidx).*FTIME;
        end
        cviei = std(ieis)/mean(ieis);
        %     eventrate = numel(onsetidx)/(numel(tmpdftrace)*FTIME); % Eventrate calculated from number of events [Hz]
        eventrate = 1/median(ieis); % Eventrate calculated from intereventintervals [Hz]
    else
        ieis= NaN;
        cviei= NaN;
        eventrate = NaN;
        onsetidx = NaN;
        peaks = zeros(size(onset));
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
    evINFO(x).binaryonset = onset;
    evINFO(x).binarypeaks = peaks;
    evINFO(x).onsetidx = onsetidx;
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
    evINFO(x).randpoints = randpos;
    evINFO(x).rand_amp= randamp;
    evINFO(x).rand_thresh= randthresh;
end

%% Local functions
    function base = get_base(trc, pos)
        if pos > BASEPOINTS
            base = mean(trc(pos-BASEPOINTS:pos));
        else
            base = mean(trc(1:BASEPOINTS));
        end
    end

end