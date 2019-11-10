function [spktr,wavespktr] = trigBinaryData(SPK,PRE,POST,triggerpts,spkwaves)

wavespktr = []; 

if nargin == 5
    wavespktr = zeros(61,length(PRE:POST),length(triggerpts));
    for tr = 1:length(triggerpts)
        refwin               = floor(double(triggerpts(tr) + PRE :  triggerpts(tr) + POST));
        wavespktr(:,ismember(refwin,SPK),tr)         = spkwaves(:,ismember(SPK,refwin));
        clear refwin
    end
end

% whole trial together:
spktr = zeros(length(PRE:POST),length(triggerpts));

for tr = 1:length(triggerpts)
        refwin                           = floor(double(triggerpts(tr) + PRE :  triggerpts(tr) + POST)); 
        spktr(ismember(refwin,SPK),tr)   = 1;
        clear refwin
end




