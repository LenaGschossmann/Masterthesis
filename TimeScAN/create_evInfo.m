function [infostruct] = create_evInfo()

global ROILIST ISBGSUB FTIME

infostruct = struct();
for iTr = 1:numel(ROILIST)
    infostruct(iTr).id = ROILIST{iTr};
    infostruct(iTr).bg_sub = ISBGSUB;
    infostruct(iTr).rangein = NaN;
    infostruct(iTr).rangeout = NaN;
    infostruct(iTr).ftime = FTIME;
    infostruct(iTr).smthwin = NaN;
    infostruct(iTr).baseframes = NaN;
    infostruct(iTr).preframes = NaN;
    infostruct(iTr).postframes = NaN;
    infostruct(iTr).average = NaN;
    infostruct(iTr).sd = NaN;
    infostruct(iTr).threshold = NaN;
    infostruct(iTr).thresholdtype = NaN;
    infostruct(iTr).peakthreshold = NaN;
    infostruct(iTr).binaryonset = NaN;
    infostruct(iTr).binarypeaks = NaN;
    infostruct(iTr).onsetidx = NaN;
    infostruct(iTr).peakidx = NaN;
    infostruct(iTr).amps = NaN;
    infostruct(iTr).amps_dFoF = NaN;
    infostruct(iTr).avamp = NaN;
    infostruct(iTr).avamp_dFoF = NaN;
    infostruct(iTr).ieis = NaN;
    infostruct(iTr).aviei = NaN;
    infostruct(iTr).cviei = NaN;
    infostruct(iTr).eventratetype = 'IEI based';
    infostruct(iTr).eventrate = NaN;
    infostruct(iTr).baselinevalues = NaN;
    infostruct(iTr).risetimes = NaN;
    infostruct(iTr).decaytimes = NaN;
    infostruct(iTr).accepted = 'pending'; % values: pending, accepted
end
end