function [spktr,spkcyc,per] = triggerDRFTdata(SPK,STIM,spkTPs,all_cyc,fs,pre,post,photo_yes)

per    = []; 
spkcyc = [];
spktr  = []; 
if photo_yes
    
    tf      =  STIM.temporal_freq(1);
    per    =  floor(1000./tf).*(fs./1000);
    
    % separate by cycle:
    spkcyc = nan(per + 1,max(cellfun(@length,all_cyc)),size(all_cyc,2));
    for tr = 1:size(all_cyc,2)
        if ~isempty(all_cyc{tr})
            for c = 1:length(all_cyc{tr})
                clear refwin
                spkcyc(:,c,tr)   = 0;
                refwin           = all_cyc{tr}(c) :  all_cyc{tr}(c) + per;
                id               = find(ismember(refwin,SPK(SPK>refwin(1) & SPK<refwin(end))));
                spkcyc(id,c,tr)  = 1;
                
            end
        end
    end
    
    % whole trial together:
    spktr = nan(length([pre:post]),size(all_cyc,2));
    for tr = 1:size(all_cyc,2)
        if ~isempty(all_cyc{tr})
            spktr(:,tr)          = 0;
            refwin               = all_cyc{tr}(1) + pre :  all_cyc{tr}(1) + post;
            id                   = find(ismember(refwin,SPK(SPK>refwin(1) & SPK<refwin(end))));
            spktr(id,tr)         = 1;
            clear refwin
        end
    end
    
else
    
    % if you don't want photodiode signal triggering, use just event codes
    
    spktr  = nan(length([pre:post]),length(STIM.tilt));
    spkcyc = []; 
    for tr = 1:size(spkTPs,1)
        if ~isempty(spkTPs(tr))
            spktr(:,tr)          = 0;
            refwin               = spkTPs(tr,1) + pre :  spkTPs(tr,1) + post;
            id                   = find(ismember(refwin,SPK(SPK>refwin(1) & SPK<refwin(end))));
            spktr(id,tr)         = 1;
            clear refwin
        end
    end
    
     
end

end







       