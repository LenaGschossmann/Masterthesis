function [infostruct] = create_evInfo()

global SMTHWIN TRACEID TRACEDATA ROILIST RANGE TOTP ISBGSUB FTIME evINFO...
    BASEPOINTS PREFRAMES POSTFRAMES

infostruct = struct();
for iTr = 1:numel(ROILIST)
    infostruct(iTr).id = ROILIST{iTr};
    infostruct(iTr).bg_sub = ISBGSUB;
    infostruct(iTr).rangein = NaN;
    infostruct(iTr).rangeout = NaN;
    infostruct(iTr).ftime = FTIME;
    infostruct(iTr).smthwin = NaN;
    infostruct(iTr).basepoints = NaN;
    infostruct(iTr).preframes = NaN;
    infostruct(iTr).postframes = NaN;
    infostruct(iTr).average = NaN;
    infostruct(iTr).sd = NaN;
    infostruct(iTr).threshold = NaN;
    infostruct(iTr).peakthreshold = NaN;
    infostruct(iTr).binarycross = NaN;
    infostruct(iTr).binarypeaks = NaN;
    infostruct(iTr).crossidx = NaN;
    infostruct(iTr).peakidx = NaN;
    infostruct(iTr).amps = NaN;
    infostruct(iTr).avamp = NaN;
    infostruct(iTr).eventtype = NaN;
    infostruct(iTr).ieis = NaN;
    infostruct(iTr).aviei = NaN;
    infostruct(iTr).cviei = NaN;
    infostruct(iTr).eventratetype = 'IEI based';
    infostruct(iTr).eventrate = NaN;
end
end