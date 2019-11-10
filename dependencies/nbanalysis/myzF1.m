function [zF1, f,freqrange,frid,pvec,meanpower,stdpower,powF1] = myzF1(resp,binsize,tf)

if size(resp,2) > 1
    eachtr = 1;
else
    eachtr = 0;
end

if eachtr == 1
    
    for tr = 1:size(resp,2)
        
        spkdens         = resp(:,tr);
        NFFT            = 2^nextpow2(length(spkdens));       % Next power of 2 from length of y
        Y               = fft(spkdens,NFFT)/length(spkdens); % keep units
        pvec(:,tr)      = abs(Y(1:NFFT/2+1,:));
    end
        
        f               = (1000./binsize)/2*linspace(0,1,NFFT/2+1);
        frid            = find(abs(f-tf) == min(abs(f - tf)));
        
        freqrange       = [f>=(1/tf) & f<=(1000/binsize/2)];
        
        
      
        meanpower       = nanmean(pvec(freqrange,:),1);
        stdpower        = nanstd(pvec(freqrange,:),0,1);
        powF1           = pvec(frid,:);
        
        zF1             = (powF1 - meanpower)./(stdpower);
    
 
else
    
    spkdens   = resp;
    NFFT      = 2^nextpow2(length(spkdens)); % Next power of 2 from length of y
    Y         = fft(spkdens,NFFT)/length(spkdens); % keep units
    
    f         = (1000./binsize)/2*linspace(0,1,NFFT/2+1);
    frid      = find(abs(f-tf) == min(abs(f - tf)));
    
    freqrange = [f>=(1/tf) & f<=(1000/binsize/2)];
    
    
    pvec      = abs(Y(:,1:NFFT/2+1));
    meanpower = mean(pvec(:,freqrange));
    stdpower  = nanstd(pvec(:,freqrange),0,2);
    powF1     = pvec(:,frid);
    
    zF1       = (powF1 - meanpower)./(stdpower);
    
end


end