%% get matrix for color characterization (for L,M,S cone stimuli)
addpath('/users/kaciedougherty/documents/code/mlanalysis')
addpath('/users/kaciedougherty/documents/code/processdata_code')
clear all 
close all 
pack

% testing code!
RGUNsamps = normrnd(620,220,[2000 1]);
GGUNsamps = normrnd(550,200,[2000 1]);
BGUNsamps = normrnd(450,200,[2000 1]);

rh   = histfit(RGUNsamps,20); 
wvr  = floor(rh(2).XData); froi_r = wvr(wvr>=380 & wvr<=780); 
denr = rh(2).YData; denoi_r = calcNormData(denr(wvr>=380 & wvr<=780)','within')';
close 

gh   = histfit(GGUNsamps); 
wvg  = floor(gh(2).XData); froi_g  = wvg(wvg>=380 & wvg<=780);
deng = gh(2).YData; denoi_g = calcNormData(deng(wvg>=380 & wvg<=780)','within')';
close

bh   = histfit(BGUNsamps); 
wvb  = floor(bh(2).XData); froi_b = wvb(wvb>=380 & wvb<=780);   
denb = bh(2).YData; denoi_b = calcNormData(denb(wvb>=380 & wvb<=780)','within')'; 
close; 

RGUN = denoi_r'; GGUN = denoi_g'; BGUN = denoi_b'; 
if size(denoi_r) ~= size(denoi_g) ~= size(denoi_b)
    RGUN = RGUN(1:min([length(denoi_r) length(denoi_g) length(denoi_b)])); 
    GGUN = GGUN(1:min([length(denoi_r) length(denoi_g) length(denoi_b)]));
    BGUN = BGUN(1:min([length(denoi_r) length(denoi_g) length(denoi_b)]));
end

%%
samp        = 4;                % sampling resoultion nm 
wavelengths = [1:length(RGUN)]; % [380:samp:780];

% get cone spectra (S) and name emission spectra (P):
[conespectraSampled]=Stockman_Sharpe_cone_fundamentals();
S = conespectraSampled(1:length(wavelengths),1:3)'...
    *samp;                      % conespectraSampled(:,1:3); clear conespectraSampled;
P = [RGUN GGUN BGUN]./1000;     % rows are wavelength, columns are guns (emission spectra)


w = [1 0 0; 0 1 0; 0 0 1];      % desired phosophor intensities, w
c = P*w;                        % spectrum of colored stimulus, c

% calculate cone coordinates:
s = S*c; 





