function [tedges,resp,binsize] = myPSTH(triallen,PRE,thesetrs,r_spktintr,r_trls,fs,binsize,varargin)

if ~exist('binsize')
    binsize  = 5; % ms
end

eachtr = 1; % run psth on each trial 


if eachtr == 1
    
    nthesetr = find(thesetrs);
    ntrials  = length(find(thesetrs));
    lastBin  = binsize * ceil((triallen-1)*(1000/(fs*binsize)));
    edges    = 0 : binsize : lastBin;
    tedges   = edges + PRE;
    
    for tr = 1:ntrials
        
        x            = (mod(r_spktintr(ismember(r_trls,nthesetr(tr)))-1,triallen)+1)*(1000/fs);
        
        resp(:,tr)   = (histc(x,edges)) / (1*binsize) .* (1000/binsize);
    end
    
    
else
    

    ntrials  = length(find(thesetrs));
    lastBin  = binsize * ceil((triallen-1)*(1000/(fs*binsize)));
    edges    = 0 : binsize : lastBin;
    tedges   = edges + PRE;
    x        = (mod(r_spktintr(ismember(r_trls,find(thesetrs)))-1,triallen)+1)*(1000/fs);
    
    resp     = (histc(x,edges)) / (ntrials*binsize) .* (1000/binsize);
    
end

end