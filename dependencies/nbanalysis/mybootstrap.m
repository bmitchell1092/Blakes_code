function [values,percentiles] = mybootstrap(data,reps)

ntrials = length(data); 
for r = 1:reps
    
    resamp     = datasample(data,ntrials,'replace',true); 
    samples(r) = mean(resamp); 
      
end

percentiles    = [.025 .05 .95 .975]; 
[values]       = quantile(samples,percentiles); 
