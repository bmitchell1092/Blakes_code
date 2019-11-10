function [prd,xprd] = runCRFFit(spkdata,unqc)

global Data
Data = [unqc.*100 spkdata]';
[Gr Gc q s b fbest] = fitCRdata;
xprd = [0:100];
for pc = 1:length(xprd)
    prd(pc) = Gr*[xprd(pc)^(q+s)]/[xprd(pc)^q + Gc^q]+b; % prediction
end

