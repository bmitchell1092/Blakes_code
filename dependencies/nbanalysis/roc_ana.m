function [AUC, significance] = roc_ana(data,num_of_intervals)

pref    = data(:,1);  
nonpref = data(:,2);  


if nargin <2
num_of_intervals  = 12; % I?d play around/try different values to see what makes the most sense (default: same as thresholds)
end

num_of_thresholds = num_of_intervals + 1; % I'd play around/try different values to see what makes the most sense (Pasternak et al. used 12)

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
pref_hist     = prefhist;
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
    hits(criterion)               = sum(prefhist([binpos:lastbin]));
end   
    
       
% INTEGRATION
AUC = abs(trapz(falsealarms,hits)); 


numperm                       = 10000;
[~,shuffled_distr]            = new_sigtest(AUC,pref,nonpref,numperm,num_of_intervals,num_of_thresholds); 


crit_85                       = numperm*.85;
sig_85                        = AUC > shuffled_distr(crit_85);  

crit_90                       = numperm*.9;
sig_90                        = AUC > shuffled_distr(crit_90);  

crit_95                       = numperm*.95;
sig_95                        = AUC > shuffled_distr(crit_95);  

crit_99                       = numperm*.99;
sig_99                        = AUC > shuffled_distr(crit_99); 

crit_999                      = numperm*.999;
sig_999                       = AUC > shuffled_distr(crit_999); 


significance.sig_85           = sig_85;
significance.sig_90           = sig_90; 
significance.sig_95           = sig_95; 
significance.sig_99           = sig_99; 
significance.sig_999          = sig_999; 