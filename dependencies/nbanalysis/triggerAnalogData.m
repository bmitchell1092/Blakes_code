function [trneural,cycneural] = triggerAnalogData(neural,pre,post,r_spkTPs,r_all_cyc,fs,per,varargin)

% trigger analog data

for tr = 1:length(r_spkTPs)
    
    if size(neural,2) > 1
        refwin           = floor([r_spkTPs(tr)./(fs./1000) + pre : r_spkTPs(tr)./(fs./1000) + post]);
        trneural(:,:,tr) = neural(refwin,:);
    else
        refwin           = floor([r_spkTPs(tr)./(fs./1000) + pre : r_spkTPs(tr)./(fs./1000) + post]);
        trneural(:,tr)   = neural(refwin,:);
    end
    
end


% trigger by cycle for drifting gratings 

if ~isempty(r_all_cyc)
    sizes = cellfun(@(C) size(C,2), r_all_cyc);
    ncyc  = min(sizes(sizes>0));
    
    for tr = 1:length(r_spkTPs)
        
        for c = 1:ncyc
            if size(neural,2) > 1
                refwin              = floor([r_all_cyc{tr}(c)./(fs/1000) : r_all_cyc{tr}(c)./(fs/1000) + (per./(fs/1000))]);
                cycneural(:,:,c,tr) = neural(refwin,:);
            else
                refwin              = floor([r_all_cyc{tr}(c)./(fs/1000) : r_all_cyc{tr}(c)./(fs/1000) + (per./(fs/1000))]);
                cycneural(:,c,tr)   = neural(refwin,:);
            end
        end
        
    end
end