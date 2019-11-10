function [pEvC, pEvT] = parsEventCodesML(EventCodes,EventTimes,neighteen)
% Sept 2014, MAC
% output same structure as BHV.CodeNumbers & BHV.CodeTimes
% Note however, that BHV.CodeTimes is from each trail start while pEvT is
% from start of BR recording (more useful fro trigering data)


% these first 2 if statements allow function to be run as:
% parsEventCodesML(NEV.Data.SerialDigitalIO.UnparsedData,NEV.Data.SerialDigitalIO.TimeStampSec)
if nargin == 2
    neighteen = 3; 
end

if ~any(EventCodes == 9)
   EventCodes = EventCodes - 128;
end
if  EventTimes(1)<1
    EventTimes = EventTimes * 1000; %ms, to match 1kHz
end

% get trial start index
stind = find(EventCodes == 9);
d = diff(diff(stind)); d = [NaN; d; NaN];
stind = stind(d ==0);
ntr   = length(stind);
if EventCodes(end) ~= 18 % data collection was stoped within a trial
    ntr = ntr-1;
end

% setup output vars
pEvC = cell(1,ntr);
pEvT = cell(1,ntr);


for tr = 1:ntr
    ind = stind(tr)-1;
    ct = 0;
    while EventCodes(ind) ~= 18
        ct = ct + 1;
        
        pEvC{tr}(ct,1) = EventCodes(ind);
        pEvT{tr}(ct,1) = EventTimes(ind);
        
        ind = ind+1;
    end
    % get final 18s
    for ending = 1:neighteen
        pEvC{tr}(ct+ending,1) = EventCodes(ind+ending-1);
        pEvT{tr}(ct+ending,1) = EventTimes(ind+ending-1);
    end

end