function [UNIT] = removeTrialsfromStruct(UNIT,rmtr)

fds     = fields(UNIT); 
start_n = size(UNIT.onsets,1);

for fd = 1:length(fds)
    if ismember({'sdfcyc','cycn','trn'},fds{fd}), continue, end
    if size(UNIT.(fds{fd}),1) == start_n
        UNIT.(fds{fd})(rmtr) = [];
    elseif size(UNIT.(fds{fd}),2) == start_n
        UNIT.(fds{fd})(:,rmtr) = [];
    elseif size(UNIT.(fds{fd}),3) == start_n
         UNIT.(fds{fd})(:,:,rmtr) = [];
    end
end

if any(ismember({'sdfcyc','cycn','trn'},fds))
    rmid = ismember(UNIT.trn,rmtr); 
    UNIT.trn(rmid) = []; 
    UNIT.sdfcyc(:,rmid) = []; 
    UNIT.cycn(:,rmid) = []; 
end
