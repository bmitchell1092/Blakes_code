function [spktr] = triggerPointData(PRE,POST,onsets,pins,NEV,SPK)

if nargin < 6
    
    spktr = zeros(length([PRE:POST]),length(pins),length(onsets));
    
    for ch = 1:length(pins)
        
        SPK    = double(NEV.Data.Spikes.TimeStamp(find(NEV.Data.Spikes.Electrode == pins(ch))));
        for tr = 1:length(onsets)
            refwin = onsets(tr) + PRE : onsets(tr) + POST;
            id                   = find(ismember(refwin,SPK(SPK>=refwin(1) & SPK<=refwin(end))));
            spktr(id,ch,tr)      = 1;
        end
    end
    
else
    
   
    spktr = zeros(length([PRE:POST]),length(onsets));
    
    for tr = 1:length(onsets)
        refwin = onsets(tr) + PRE : onsets(tr) + POST;
        id                   = find(ismember(refwin,SPK(SPK>=refwin(1) & SPK<=refwin(end))));
        spktr(id,tr)         = 1;
    end
    

end

