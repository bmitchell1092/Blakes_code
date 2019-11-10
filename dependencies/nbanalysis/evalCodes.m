
function [codesbtw tcodesbtw trackcond nine eighteen succtr] = evalCodes(EventCodes,EventTimes,succ) 
% find trial start and end times 

if nargin == 2
    succ = 1; 
end

EventCodes = EventCodes(find(EventCodes==9,1,'first'):end); %1st nine
EventCodes = EventCodes(1:find(EventCodes==18,1,'last'));   %last nine, truncate


idxvec     = 1:numel(EventCodes);
nine       = idxvec(EventCodes==9);
eighteen   = idxvec(EventCodes==18);

if mod(length(nine),3) ~=0
    error('error in event codes\');
else
    nine = nine(1:3:end);
 
end

if mod(length(eighteen),3) ~=0
    error('error in event codes\');
else
    eighteen = eighteen(1:3:end);
      
end

if length(eighteen) ~= length(nine) & eighteen(end) < nine(end)
    nine(end) = []; 
    
elseif length(eighteen) ~= length(nine) & eighteen(end) < nine(end)
end
      

% make and sort correct trials by condition
mxevent   = 116; % last ML event code

idxtr = 0;
trackcond = []; 
succtr = []; 
for tr = 1:length(nine)
    
    h_codesbtw = EventCodes(nine(tr):eighteen(tr));
    
    % successful trial?
    succtrls = find(h_codesbtw == 96);
    
    if succ == 1 & isempty(succtrls)
        continue
    else
    
    idxtr = idxtr + 1; %track number of successful trs for indexing
    codesbtw{idxtr}  =  double(EventCodes(nine(tr):eighteen(tr))); %hang onto codes for trial
    
    tcodesbtw{idxtr} =  double(EventTimes(nine(tr):eighteen(tr))); %hang onto codes for trial
    
    succtr = [succtr tr]; 
    end
    clear h_codesbtw h_c 
end


