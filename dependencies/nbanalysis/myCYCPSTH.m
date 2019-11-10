function [tedges,resp,binsize] = myCYCPSTH(triallen,thesetrs,r_spkcyc,fs,binsize,varargin)

if ~exist('binsize')
    binsize  = 5; % ms
end

eachtr = 1; % run psth on each trial/cyc

cyctrs = reshape([r_spkcyc(:,:,thesetrs)],[size(r_spkcyc,1) length(find(thesetrs))*size(r_spkcyc,2)]);

if eachtr == 1
    
    ntrials  = size(cyctrs,2);
    lastBin  = binsize * ceil((triallen-1)*(1000/(fs*binsize)));
    edges    = 0 : binsize : lastBin;
    tedges   = edges;
    
    for tr = 1:ntrials
        
        x        = (mod((find(cyctrs(:,tr)))-1,triallen)+1)*(1000/fs);
        
    
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