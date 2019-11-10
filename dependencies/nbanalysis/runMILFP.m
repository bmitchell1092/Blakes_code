
function [MI] = runMILFP(lowdat,highdat,Fs)

totchan = size(lowdat,2);
tottr   = size(lowdat,3); 


binwidth = 30;
bins     = [0:binwidth:360];


for tr = 1:size(highdat,3)
for LFchan = 1:size(lowdat,2)
    for HFchan = 1:size(highdat,2)
        
        LFphase = rad2deg(angle(hilbert(lowdat(:,LFchan,tr))));
        LFphase(LFphase < 0) = 360 + LFphase(LFphase < 0 );
     
        [~,phIDX] = histc(LFphase,bins);
        [srtph,h_srtphIDX] = sort(phIDX);
        srtphIDX = floor(h_srtphIDX);
        srtphIDX(srtphIDX>size(highdat,1)) = size(highdat,1);
        if length(find(srtphIDX == size(highdat,1) > 3))
            error('error in phase index length\');
        end
        
        %HF amplitude:
        AmpEnv    = abs(hilbert(highdat(:,HFchan,tr))); 
        srtAmpEnv = AmpEnv(srtphIDX);
        
        
        for bn = 1:length(bins) - 1
            
            t_binAmpEnv = srtAmpEnv;
            t_binAmpEnv(srtph ~= bn) = NaN;
            
            binAmpEnv(:,bn,LFchan,HFchan,tr) = nanmean(t_binAmpEnv,1);
            
            
        end
        
        clear AmpEnv srtph h_srtphIDX phIDX t_binAmpEnv rand_t_binAmpENv ;
    end
end
end


mn_binAmpEnv         = nanmean(squeeze(binAmpEnv),4);
%UniformMeanBinAmpEnv = repmat(mean(mn_binAmpEnv,1),length(bins),1);

%Calculate H

Hmax=log10(size(mn_binAmpEnv,1));

%hfr is amp channel. lfr is phase channel
for lfr=1:size(mn_binAmpEnv,2)
    for hfr=1:size(mn_binAmpEnv,3)
        
        
        sum_amps=sum(mn_binAmpEnv(:,lfr,hfr),1);
        rep_sum=repmat(sum_amps,1,size(mn_binAmpEnv,1));
        p=squeeze(mn_binAmpEnv(:,lfr,hfr))./squeeze(rep_sum)';
        
        H=-1*sum(p.*log10(p)); %changed to log10 3/14
        
        MI(:,lfr,hfr)=(Hmax-H)./Hmax;
        clear sum_amps rep_sum p H
    end
end



end