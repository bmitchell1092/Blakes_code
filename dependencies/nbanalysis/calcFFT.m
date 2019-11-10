function [power,freq] = calcFFT(data)

power = [];
freq  = [];
NFFT            = 2^nextpow2(size(data,1));       % Next power of 2 from length of y
for tr = 1:size(data,2)
    
    spkdens         = data(:,tr);
    Y               = fft(spkdens,NFFT)/length(spkdens); % keep units
    power(:,tr)     = abs(Y(1:NFFT/2+1,:));
end
freq                = (1000./1)/2*linspace(0,1,NFFT/2+1);