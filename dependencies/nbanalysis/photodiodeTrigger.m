function [STIM] = photodiodeTrigger(filename,STIM,refresh)


[BNC, BNCLabels,fs] = getBNCData({'ainp1'},strcat(filename,'.ns6'));

refreshdelay = fs./refresh +200;
STIM.onsets_p = [];
for tr = 1:length(STIM.trial)
    trstart      = floor((STIM.trstart(tr) + (500.*(fs/1000)))); % only consider samples > 50 ms in (there's an initial flash)
    trend        = floor(STIM.trend(tr));
    refwin       = trstart:trend;
    sig          = BNC(refwin,1);
    thr          = mean(sig) -  2*std(sig,0,1);
   
    if ~isfield('temporal_freq',STIM)
        detected  = refwin(find(sig<thr,1,'first')); % onset of first cycle
    elseif STIM.temporal_freq(tr) <= 1
         detected  = refwin(find(sig<thr,1,'first')); % onset of first cycle
    else
        allp   =  find(sig<thr);
        pts    = refwin(allp(diff(allp)>1));
        detected = [pts(1) pts((find(diff(pts)>(refreshdelay)))+1)];
    end
    
    if isempty(detected)
        
        STIM.onsets_p(tr) = nan; 
    else
    
    STIM.onsets_p(tr) = detected(1);
    end
    
end

if size(STIM.onsets_p,2) > size(STIM.onsets_p,1)
    STIM.onsets_p = STIM.onsets_p';
end