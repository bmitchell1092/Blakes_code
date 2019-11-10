function [significance,shuffled_distr] = new_sigtest(cp,pref,nonpref,numperm,num_of_intervals,num_of_thresholds)
% 030504 - AM
           
% reshuffle distributions
rand('state',sum(100*clock)); % reset randomizer

prefsize    = length(pref);
nonprefsize = length(nonpref);
 
% first, we create a pool of all response values:
allvals = vertcat(pref,nonpref);
% now, we shuffle this guy...
nn=length(allvals);
% n-times
for permutnum = 1:numperm;
    shuffledvals=allvals;  
    for i=1:nn-1,
        j=floor(rand(1,1)*(nn-i+1)+i);
        temp=shuffledvals(i);
        shuffledvals(i)=shuffledvals(j);
        shuffledvals(j)=temp;
    end
    % then, we divide up into the two distributions...
    pref_shuffled = shuffledvals(1:prefsize);
    nonpref_shuffled = shuffledvals((prefsize+1):(prefsize+nonprefsize));
    % and compute cp again...
    cp_shuffle(permutnum) = compute_cp(pref_shuffled,...
        nonpref_shuffled,num_of_intervals,num_of_thresholds);
end % permutation loop
 
shuffled_distr = sort(cp_shuffle);
if cp >= 0.5 % see, whether it is exceeding highest vals...
    ninefivecrit = numperm*.975;
    max_critval  = shuffled_distr(ninefivecrit);
    if cp > max_critval
        significance = 1;
    elseif cp <= max_critval
        significance = 0;  
    end
elseif cp < 0.5 % smaller than smallest ?
    ninefivecrit = numperm*.975;
    ninefivecrit = numperm-ninefivecrit;
    min_critval  = shuffled_distr(ninefivecrit);
    if cp < min_critval
        significance = 1;
    elseif cp >= min_critval
        significance = 0;   
    end
end  