function runTwoGrpDataSim(nsim,datagrp1,datagrp2,tlabel,varargin)


sdatagrp1     = []; 
sdatagrp2     = [];
diffresp      = []; 

actualdiff    = mean(datagrp1,2) - mean(datagrp2,2);
alphaval      = 0.05; 

for ns = 1:nsim
    
drawgrp1      = datasample(datagrp1,size(datagrp1,2),2);  
drawgrp2      = datasample(datagrp2,size(datagrp1,2),2);  

sdatagrp1     = [sdatagrp1 mean(drawgrp1,2)];
sdatagrp2     = [sdatagrp2 mean(drawgrp2,2)];

diffresp      = [diffresp mean(drawgrp1,2)- mean(drawgrp2,2)]; 
 
end

grp1diffdist  = [];
grp2diffdist  = [];

for ns = 1:floor(nsim./2)
    
    grp1pts       = datasample(sdatagrp1,2,2);
    grp1diffdist = [grp1diffdist diff([grp1pts])];
    
    grp2pts      = datasample(sdatagrp2,2,2);
    grp2diffdist = [grp2diffdist diff([grp2pts])];
end
combinednull = [grp1diffdist grp2diffdist]; 


% sdatagrp1  = datasample(datagrp1,nsim,2);  
% sdatagrp2  = datasample(datagrp2,nsim,2);  
% diffresp   = [sdatagrp1 -sdatagrp2]; 
%  
% 
% grp1diffdist = [];
% grp2diffdist = [];
% 
% for ns = 1:floor(nsim./2)
%     
%     grp1pts       = datasample(sdatagrp1,2,2);
%     grp1diffdist = [grp1diffdist diff([grp1pts])];
%     
%     grp2pts      = datasample(sdatagrp2,2,2);
%     grp2diffdist = [grp2diffdist diff([grp2pts])];
% end
% combinednull = [grp1diffdist grp2diffdist]; 

figure, set(gcf,'Color','w'); 
subplot(1,2,1) 
hist(sdatagrp1); 
h = findobj(gca,'Type','patch'); 
h.FaceColor = [.2 .5 1]; 
hold on; 
hist(sdatagrp2)
set(gca,'Box','off','TickDir','out'); 
xlabel('responses'); ylabel('counts'); 
legend(gca,{'monocular','binocular'}); 
if exist('tlabel')
    title(gca,tlabel);
end

subplot(1,2,2) 
hist(combinednull); 
h = findobj(gca,'Type','patch'); 
h.FaceColor = [.2 .5 1]; 
hold on; 
hist(diffresp); 
hold on; 
set(gca,'Box','off','TickDir','out'); 
xlabel('monoc-binoc resp'); ylabel('counts')
labels = {'null diff','simulated diffresp'}; 
legend(gca,labels);
if exist('tlabel')
    title(gca,tlabel);
end

[npool,meanpool,stdpool] = pooledmeanstd(length(sdatagrp2),mean(sdatagrp2),std(sdatagrp2,0,2),length([sdatagrp1]),mean([sdatagrp1]),std([sdatagrp1],0,2));

effectsize   = (mean(sdatagrp1) - mean(sdatagrp2))./stdpool;
power        = 0.6; 
n            = sampsizepwr('t',[mean(sdatagrp2) std(sdatagrp2,0,2)],mean(sdatagrp1),power); 

combinednull = sort(combinednull,'ascend'); 
cutoff       = combinednull(length(combinednull).*(1-alphaval)); 
v  = vline(cutoff); set(v,'Color','r','LineWidth',2);
hold on;
v2 = vline(actualdiff); set(v2,'Color','k','LineWidth',2,'LineStyle','-');
if exist('tlabel')
    fprintf('for %s\n',tlabel)
end

fprintf('\n effect size : %f\n',effectsize);
fprintf('\n n trials    : %f\n',n);

