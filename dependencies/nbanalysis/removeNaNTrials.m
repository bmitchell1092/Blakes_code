function [r_spktr,r_spkcyc,r_spktintr,r_trls,r_spktincyc,r_cycn,r_cyctn,r_STIM] = removeNaNTrials(spktr,spkcyc,STIM)

% remove trials that contain NaNs for plotting (should be unrewarded trials if
% flagRewardedTrialsOnly == 1 or badobs trials)

ct          = 0;
flds        = fields(STIM);
r_spktr     = []; 
r_spkcyc    = []; 
r_spktintr  = []; % indices of spike times relative to trial window
r_trls      = []; % trial n for each spike in r_spktintr
r_spktincyc = []; % indices of spike times relative to cycle window
r_cycn      = [];
r_cyctn     = [];
r_STIM      = [];

for i = 1:numel(flds)
    r_STIM.(flds{i}) = [];
end

if ~isempty(spkcyc)
    
    %photodiode version:
    for tr = 1:size(spktr,2)
        
        if any(isnan(spktr(:,tr)))
            continue
            
        else
            
            ct                   = ct + 1;
            r_spktr(:,ct)        = spktr(:,tr);
            r_spkcyc(:,:,ct)     = spkcyc(:,:,tr);
            
            r_spktintr           = [r_spktintr find(spktr(:,tr))'];
            r_trls               = [r_trls repmat(ct,1,length(find(spktr(:,tr))))];
            
            for c = 1:size(spkcyc,2)
                r_spktincyc      = [r_spktincyc find(spkcyc(:,c,tr))'];
                r_cycn           = [r_cycn repmat(c,1,length(find(spkcyc(:,c,tr))))];
                r_cyctn          = [r_cyctn repmat(ct,1,length(find(spkcyc(:,c,tr))))];
            end
            
            for i = 1:numel(flds)
                r_STIM.(flds{i}) = [r_STIM.(flds{i}) STIM.(flds{i})(tr)];
            end
        end
        
    end
    
else
    
    % event code version
    for tr = 1:size(spktr,2)
        
        if any(isnan(spktr(:,tr)))
            continue
            
        else
            
            ct               = ct + 1;
            r_spktr(:,ct)    = spktr(:,tr);
            r_spktintr       = [r_spktintr find(spktr(:,tr))'];
            r_trls           = [r_trls repmat(ct,1,length(find(spktr(:,tr))))];
            
     
            for i = 1:numel(flds)
                r_STIM.(flds{i}) = [r_STIM.(flds{i}) STIM.(flds{i})(tr)];
            end
            
        end
        
    end
    
    
end


end