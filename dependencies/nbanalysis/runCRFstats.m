function [anovapvals,ttestpvals] = runCRFstats(unqc,nunqc,spkdata,tvec)

%anova (5x2): Compare responses at each contrast level
%and for each ocular conditions

y = []; g1 = []; g2= [];

for c = 1:length(unqc)
    for fc = 1:length(nunqc)
        if fc == 1
            
            y  = [y nanmean(spkdata{c,1}(tvec > 0,:),1)];
            g1 = [g1; repmat(fc,size(spkdata{c,1},2),1)]; % ocular condition
            g2 = [g2; repmat(c,size(spkdata{c,1},2),1)];  % contrast level
            
        else
            
            y  = [y nanmean(spkdata{c,end}(tvec > 0,:),1)];
            g1 = [g1; repmat(fc,size(spkdata{c,end},2),1)]; % ocular condition
            g2 = [g2; repmat(c,size(spkdata{c,end},2),1)];  % contrast level
            
        end
      
        
    end
end
if all(isnan(g1)) | all(isnan(g2))
    anovapvals = nan;
else
    anovapvals = anovan(y',{g1,g2});
end



% t-tests at each level (This is really to check for difference in max saturation
% point i.e., at the highest contrast levels is there a difference?)

for c = 1:length(unqc)
    
    clear p
    
    [~,p] = ttest2(nanmedian(spkdata{c,1}(tvec > UNITDAT.latency(conservative(1)),:),1),nanmedian(spkdata{c,end}(tvec > UNITDAT.latency(conservative(1)),:),1),'tail','left');
    
    if all(isnan(nanmean(spkdata{c,1}(tvec > 0,:),1))) | all(isnan(nanmean(spkdata{c,end}(tvec > 0,:),1)))
        ttestpvals(c) = nan;
    else
        ttestpvals(c) = p;
    end
    
    
end
