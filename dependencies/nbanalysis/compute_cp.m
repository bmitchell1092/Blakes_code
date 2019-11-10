function [auc] = compute_cp(pref, nonpref,num_of_intervals,num_of_thresholds)
 
pref_min     = min(pref);   
nonpref_min  = min(nonpref);
pref_max     = max(pref);
nonpref_max  = max(nonpref);   
startbin     = min(pref_min,nonpref_min);
stopbin      = max(pref_max,nonpref_max);
range        = stopbin-startbin;
bininterval  = range/num_of_intervals;
edge         = [startbin:bininterval:stopbin];

prefhist = histc(pref(:,1),edge);
[r,c]= size(prefhist);
if c>1
    fprintf('\n rotating matrice');
    prefhist = prefhist';
end
nonprefhist = histc(nonpref(:,1),edge);
[r,c]= size(nonprefhist);
if c>1
    fprintf('\n rotating matrice');
    nonprefhist = nonprefhist';
end   
pref_hist  = prefhist;
nonpref_hist  = nonprefhist;      
    
% THRESHOLD  
prefhist     = prefhist/sum(prefhist); % normalize number of trials
nonprefhist  = nonprefhist/sum(nonprefhist); % normalize number of trials
bestfit      = round(length(prefhist)/num_of_thresholds);
act_numofcrit     = round(length(prefhist)/bestfit);
if act_numofcrit ~= num_of_thresholds
    fprintf('\n sorry, but we have to use %d criteria',act_numofcrit);   
end
criterionsize = bestfit-1;
lastbin = length(prefhist);       
for criterion = 1:act_numofcrit    
    binpos = criterion+criterionsize;
    falsealarms(criterion) = sum(nonprefhist([binpos:lastbin]));
    hits(criterion)        = sum(prefhist([binpos:lastbin]));
end   
          
% INTEGRATION
auc = abs(trapz(falsealarms,hits)); 
