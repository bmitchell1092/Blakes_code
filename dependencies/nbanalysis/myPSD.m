function [pvec,f] = myPSD(neural)

for tr = 1:size(neural,2)

        resp            = neural(:,tr);
        NFFT            = 2^nextpow2(length(resp)); % Next power of 2 from length of y
        Y               = fft(resp,NFFT)/length(resp); % keep units
        pvec(:,tr)      = abs(Y(1:NFFT/2+1,:));
       
  
end

      f               = (1000)/2*linspace(0,1,NFFT/2+1);

        
end