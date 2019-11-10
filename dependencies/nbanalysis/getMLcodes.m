function [evcodes evtimes evtrials] = getMLcodes(NEV,findcodes)

% right now findcodes works only with stim onsets/offsets 
% get event codes from NEV, then clear

EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventTimes = floor(NEV.Data.SerialDigitalIO.TimeStampSec .* 1000); %ms, to match 1kHz

EventCodes = EventCodes(find(EventCodes==9,1,'first'):end); %1st nine
EventCodes = EventCodes(1:find(EventCodes==18,1,'last'));   %last eighteen, truncate

idxvec     = 1:numel(EventCodes);
nine       = idxvec(EventCodes==9);
eighteen   = idxvec(EventCodes==18);

if mod(length(nine),3) ~=0 || mod(length(eighteen),3) ~=0
    error('error in event codes\');
else
    nine = nine(1:3:end);
    eighteen = eighteen(1:3:end);
end

if length(eighteen) ~= length(nine) && eighteen(end) < nine(end)
    nine(end) = [];
    elseif length(eighteen) ~= length(nine) && eighteen(end) < nine(end)
end

% find codes and sample times
idxtr = 0;
evcodes = []; 
evtimes = []; 
evtrials = []; 
for tr = 1:length(nine)
    idxtr   = idxtr + 1;
    trcodes = double(EventCodes(nine(tr):eighteen(tr))); 
    trtimes = double(EventTimes(nine(tr):eighteen(tr)));
    existid = find(ismember(trcodes,findcodes)); 
    existcodes = trcodes(find(ismember(trcodes,findcodes))); 
    if mod(length(existcodes),2) ~= 0
        existcodes(end) = []; % keep only stim onset/offset pairs
        existid(end) = []; 
    end
    if ~isempty(existcodes)
    evcodes = [evcodes; existcodes];
    evtimes = [evtimes trtimes(existid)]; 
    evtrials = [evtrials; repmat(tr,length(existcodes),1)]; 
    end
    clear existcodes trcodes trtimes
end




