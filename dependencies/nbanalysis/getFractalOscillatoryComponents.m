function [Frac] = getFractalOscillatoryComponents(data)


% win       = movingwin(1)*srate;
% step      = movingwin(2)*srate;
% movingwin = [3 1]; % [window size, sliding step]

% %  load ECoG data from one sensor recorded in the left occipital of one
% %  macaque in eyes-closed awake state, totally 5 mins
% %  load('ECoG_data.mat');

% separate fractal and oscillatory components using sliding window
% nwin = floor((length(data) - win)/step);
% sig  = zeros(win,nwin);
% for i = 1 : nwin
%     sig(:,i) = data(ceil((i-1)*step)+1 : ceil((i-1)*step)+win);
% end

srate     = 1000;    % sampling frequency
frange    = [1 100];
tic
Frac = amri_sig_fractal(data,srate,'detrend',1,'frange',frange);
%Frac.time = (0:step/srate:step*(nwin-1)/srate)';
toc

% fitting power-law function to the fractal power spectra
Frange = [1, 100]; % define frequency range for power-law fitting
Frac   = amri_sig_plawfit(Frac,Frange);

% % % show averaged fractal and oscillatory power spectrum
% figure;
% subplot(2,1,1);
% loglog(Frac.freq,mean(Frac.mixd,2),'b'); hold on
% loglog(Frac.freq,mean(Frac.frac,2),'r');
% subplot(2,1,2);
% plot(Frac.freq, mean(Frac.osci,2));
% hold on; 
% plot(Frac.Freq, mean(Frac.Plaw,2));
% 
% 
% figure;
% plot(Frac.freq, mean(Frac.osci(:,detrls),2)); hold on; 
% plot(Frac.freq, mean(Frac.osci(:,ndetrls),2)); hold on; 
% plot(Frac.freq, mean(Frac.osci(:,bintrls),2)); hold on; 
