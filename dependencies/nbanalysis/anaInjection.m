clear 
addpath(genpath('/Users/kaciedougherty/Documents/Code/fNPMK'))
addpath(genpath('/users/kaciedougherty/documents/code/nbanalysis/'))

preBRdatafile   = '190611_B_rfori002';
postBRdatafile  = '190611_B_rfori003';
drname          = ['/volumes/bigdrobo2/drobo2/data/rig022/' preBRdatafile(1:8) '/'];

el              = 'eA'; %'eA';
pre             = 50;
post            = 250;
chans           = 9:32; 

orange    = [238 106 73]./255;
green     = [146 163 69]./255;
blue      = [28 116 165]./255;

%% get baseline and post injection data 

clear BRdatafile; 
BRdatafile                                             = [drname preBRdatafile];
[preLFP,preSDF,presdf,prelfp,preCSD,~,preSTIM]         = getInjectionData(BRdatafile,el,chans,pre,post); 

clear BRdatafile; 
BRdatafile                                             = [drname postBRdatafile];
[postLFP,postSDF,postsdf,postlfp,postCSD,TTL,postSTIM] = getInjectionData(BRdatafile,el,chans,pre,post); 

%% remove trials when spiking comes back after injection 

rmtrials   = 268:517; 
h_postSDF  = postSDF; 
h_postCSD  = postCSD; 
h_postSTIM = postSTIM; 
h_postLFP  = postLFP ;

postSDF(:,:,rmtrials) = []; 
postCSD(:,:,rmtrials) = []; 
postLFP(:,:,rmtrials) = []; 
fds = fields(postSTIM); 
for f = 1:length(fds)
postSTIM.(fds{f})(rmtrials) = []; 
end
    
%% check for a visual response in the LFP and SDFs

pvalSDF      = nan(length(chans),1);
pval_postSDF = nan(length(chans),1);
for ch = 1:size(preSDF,2)  
    [~,pvalSDF(ch)] = ttest(squeeze(nanmean(preSDF(1:pre,ch,:),1)),squeeze(nanmean(preSDF(pre:pre+post,ch,:),1)));  
    if pvalSDF(ch) < 0.05
      [~,pval_postSDF(ch)] = ttest(squeeze(nanmean(postSDF(1:pre,ch,:),1)),squeeze(nanmean(postSDF(pre:pre+post,ch,:),1)));  
    end
end

pvalLFP      = nan(length(chans),1);
pval_postLFP = nan(length(chans),1);
for ch = 1:size(preLFP,2)
    [~,pvalLFP(ch)] = ttest(squeeze(nanmean(preLFP(1:pre,ch,:),1)),squeeze(nanmean(preLFP(pre:pre+post,ch,:),1))); 
    if pvalLFP(ch) < 0.05
      [~,pval_postLFP(ch)] = ttest(squeeze(nanmean(postLFP(1:pre,ch,:),1)),squeeze(nanmean(postLFP(pre:pre+post,ch,:),1)));  
    end
end


pval_clfp      = nan(length(chans),1);
pval_csdf = nan(length(chans),1);
for i = 1:length(chans)
    [~,pval_clfp(i)] = ttest2(abs(movmean(postlfp(:,i),400,1)),abs(movmean(prelfp(:,i),400,1)),'tail','left');
    [~,pval_csdf(i)] = ttest2(abs(movmean(postsdf(:,i),400,1)),abs(movmean(presdf(:,i),400,1)),'tail','left');
end


%% normalize data relative to pre-stimulus baseline (JUST SUBTRACTION)
clear base
base         = repmat(nanmean(nanmean(preSDF(1:abs(pre),:,:),3),1),[size(preSDF,1) 1 size(preSDF,3)]); 
sub_preSDF    = (preSDF - base); 
clear base
base         = repmat(nanmean(nanmean(postSDF(1:abs(pre),:,:),3),1),[size(postSDF,1) 1 size(postSDF,3)]); 
sub_postSDF   = (postSDF - base); 


clear base
base         = repmat(nanmean(nanmean(preCSD(1:abs(pre),:,:),3),1),[size(preCSD,1) 1 size(preCSD,3)]); 
sub_preCSD    = (preCSD - base);  
clear base
base         = repmat(nanmean(nanmean(postCSD(1:abs(pre),:,:),3),1),[size(postCSD,1) 1 size(postCSD,3)]); 
sub_postCSD   = (postCSD - base); 


clear base
base         = repmat(nanmean(nanmean(preLFP(1:abs(pre),:,:),3),1),[size(preLFP,1) 1 size(preLFP,3)]); 
sub_preLFP    = (preLFP - base); 
clear base
base         = repmat(nanmean(nanmean(postLFP(1:abs(pre),:,:),3),1),[size(postLFP,1) 1 size(postLFP,3)]); 
sub_postLFP   = (postLFP - base); 

%% normalize data relative to pre-stimulus baseline 
clear base
base         = repmat(nanmean(nanmean(preSDF(1:abs(pre),:,:),3),1),[size(preSDF,1) 1 size(preSDF,3)]); 
pc_preSDF    = (preSDF - base) ./ base; 
clear base
base         = repmat(nanmean(nanmean(postSDF(1:abs(pre),:,:),3),1),[size(postSDF,1) 1 size(postSDF,3)]); 
pc_postSDF   = (postSDF - base) ./ base; 
clear base
base         = repmat(nanmean(nanmean(preCSD(1:abs(pre),:,:),3),1),[size(preCSD,1) 1 size(preCSD,3)]); 
pc_preCSD    = (preCSD - base) ./base;  
clear base
base         = repmat(nanmean(nanmean(postCSD(1:abs(pre),:,:),3),1),[size(postCSD,1) 1 size(postCSD,3)]); 
pc_postCSD   = (postCSD - base)./base; 
clear base
base         = repmat(nanmean(nanmean(preLFP(1:abs(pre),:,:),3),1),[size(preLFP,1) 1 size(preLFP,3)]); 
pc_preLFP    = (preLFP - base) ./ base; 
clear base
base         = repmat(nanmean(nanmean(postLFP(1:abs(pre),:,:),3),1),[size(postLFP,1) 1 size(postLFP,3)]); 
pc_postLFP   = (postLFP - base) ./ base; 

%% normalize data relative to PRE-INJECTION pre-stimulus baseline 
clear base
base         = repmat(nanmean(nanmean(preSDF(1:abs(pre),:,:),3),1),[size(preSDF,1) 1 size(preSDF,3)]); 
bc_preSDF    = (preSDF - base) ./ base; 
base         = repmat(nanmean(base,3),[1 1 size(postSDF,3)]); 
bc_postSDF   = (postSDF - base) ./ base; 

clear base
base         = repmat(nanmean(nanmean(preCSD(1:abs(pre),:,:),3),1),[size(preCSD,1) 1 size(preCSD,3)]); 
bc_preCSD    = (preCSD - base) ./base;  
base         = repmat(nanmean(base,3),[1 1 size(postCSD,3)]); 
bc_postCSD   = (postCSD - base)./base; 

clear base
base         = repmat(nanmean(nanmean(preLFP(1:abs(pre),:,:),3),1),[size(preLFP,1) 1 size(preLFP,3)]); 
bc_preLFP    = (preLFP - base) ./ base; 
base         = repmat(nanmean(base,3),[1 1 size(postLFP,3)]); 
bc_postLFP   = (postLFP - base) ./ base; 

%% response before and after injection (MUA, CSD, LFP)
fds = {'preSDF','postSDF','preCSD','postCSD','preLFP','postLFP','pvalSDF','pvalLFP','pval_postSDF','pval_postLFP','tags',...
    'TTL','cSDF','cLFP','chans','pval_csdf','pval_clfp'};
for f = 1:length(fds)
    DATA.(fds{f}) = [];
end

plot_types   = {'CSD','SDF','LFP'};
DATA.prefile = preBRdatafile; 
DATA.pstfile = postBRdatafile; 
DATA.tags    = {'nA/mm^3','spks/s','microV'}; 
DATA.chans   = chans;
DATA.preSDF  = sub_preSDF;
DATA.postSDF = sub_postSDF;
DATA.preCSD  = sub_preCSD; 
DATA.postCSD = sub_postCSD;  
DATA.preLFP  = sub_preLFP; 
DATA.postLFP = sub_postLFP; 

DATA.pvalSDF      = pvalSDF; 
DATA.pvalLFP      = pvalLFP; 
DATA.pval_postSDF = pval_postSDF; 
DATA.pval_postLFP = pval_postLFP; 

plotInjectionData(DATA,pre,post,plot_types)

%%

fds = {'preSDF','postSDF','preCSD','postCSD','preLFP','postLFP','pvalSDF','pvalLFP','pval_postSDF','pval_postLFP','tags',...
    'TTL','cSDF','cLFP','chans','pval_csdf','pval_clfp'};
for f = 1:length(fds)
    DATA.(fds{f}) = [];
end

plot_types   = {'CSD','SDF','LFP'};
DATA.prefile = preBRdatafile; 
DATA.pstfile = postBRdatafile; 
DATA.tags    = {'nA/mm^3','spks/s','microV'}; 
DATA.chans   = chans;
DATA.preSDF  = preSDF;
DATA.postSDF = postSDF;
DATA.preCSD  = preCSD; 
DATA.postCSD = postCSD;  
DATA.preLFP  = preLFP; 
DATA.postLFP = postLFP; 

DATA.pvalSDF      = pvalSDF; 
DATA.pvalLFP      = pvalLFP; 
DATA.pval_postSDF = pval_postSDF; 
DATA.pval_postLFP = pval_postLFP; 

plotInjectionData(DATA,pre,post,plot_types)

%%  baseline corrected response before and after injection (MUA, CSD, LFP)
fds = {'preSDF','postSDF','preCSD','postCSD','preLFP','postLFP','pvalSDF','pvalLFP','pval_postSDF','pval_postLFP','tags',...
    'TTL','cSDF','cLFP','chans','pval_csdf','pval_clfp'};
for f = 1:length(fds)
    DATA.(fds{f}) = [];
end

plot_types   = {'CSD','SDF','LFP'};
DATA.tags    = {'% change/100 CSD','% change/100 MUA','% change/100 LFP'}; 
DATA.prefile = preBRdatafile; 
DATA.pstfile = postBRdatafile; 
DATA.chans   = chans;
DATA.preSDF  = pc_preSDF;
DATA.postSDF = pc_postSDF; 
DATA.preCSD  = pc_preCSD; 
DATA.postCSD = pc_postCSD;   
DATA.preLFP  = pc_preLFP;  
DATA.postLFP = pc_postLFP;  

DATA.pvalSDF      = pvalSDF; 
DATA.pvalLFP      = pvalLFP; 
DATA.pval_postSDF = pval_postSDF; 
DATA.pval_postLFP = pval_postLFP; 

plotInjectionData(DATA,pre,post,plot_types)

%%  
fds = {'preSDF','postSDF','preCSD','postCSD','preLFP','postLFP','pvalSDF','pvalLFP','pval_postSDF','pval_postLFP','tags',...
    'TTL','cSDF','cLFP','chans','pval_csdf','pval_clfp'};
for f = 1:length(fds)
    DATA.(fds{f}) = [];
end

plot_types   = {'CSD','SDF','LFP'};
DATA.tags    = {'pre-- % change/100 CSD','pre-- % change/100 MUA','pre-- % change/100 LFP'}; 
DATA.prefile = preBRdatafile; 
DATA.pstfile = postBRdatafile; 
DATA.chans   = chans;
DATA.preSDF  = bc_preSDF;
DATA.postSDF = bc_postSDF; 
DATA.preCSD  = bc_preCSD; 
DATA.postCSD = bc_postCSD;   
DATA.preLFP  = bc_preLFP;  
DATA.postLFP = bc_postLFP;  

DATA.pvalSDF      = pvalSDF; 
DATA.pvalLFP      = pvalLFP; 
DATA.pval_postSDF = pval_postSDF; 
DATA.pval_postLFP = pval_postLFP; 

plotInjectionData(DATA,pre,post,plot_types)

%% plot whole file irrespective of stimulus onsets 
fds = {'preSDF','postSDF','preCSD','postCSD','preLFP','postLFP','pvalSDF','pvalLFP','pval_postSDF','pval_postLFP','tags',...
    'TTL','cSDF','cLFP','chans','pval_csdf','pval_clfp'};

for f = 1:length(fds)
    DATA.(fds{f}) = [];
end

plot_types   = {'cSDF','cLFP','cTTL'};
DATA.tags    = {'spks/s','abs(microV)'}; 
DATA.cLFP    = abs(movmean(postlfp,400,1)); 
DATA.cSDF    = movmean(postsdf,400,1); 
DATA.TTL     = TTL;
DATA.chans   = chans;
DATA.pstfile = postBRdatafile; 

DATA.pval_clfp = pval_clfp; 
DATA.pval_csdf = pval_csdf; 

plotInjectionData(DATA,[],[],plot_types)


%% L4 activity!!

ch = find(chans == 26); 
figure,set(gcf,'color','w','position',[1500 500 800 400]); 
csd_group_drug = []; csd_resp = []; csd_group_eye = [];
sdf_group_drug = []; sdf_resp = []; sdf_group_eye = [];
for i = 1:4
    switch i
        case 1
            REbc = preCSD(:,ch,preSTIM.eye == 2) - nanmean(preCSD(1:pre,ch,preSTIM.eye == 2),1);
            LEbc = preCSD(:,ch,preSTIM.eye == 3) - nanmean(preCSD(1:pre,ch,preSTIM.eye == 3),1);
            BIbc = preCSD(:,ch,preSTIM.eye == 1) - nanmean(preCSD(1:pre,ch,preSTIM.eye == 1),1);
            ylab = 'naA/mm^3';
        case 2
            REbc = postCSD(:,ch,postSTIM.eye == 2) - nanmean(postCSD(1:pre,ch,postSTIM.eye == 2),1);
            LEbc = postCSD(:,ch,postSTIM.eye == 3) - nanmean(postCSD(1:pre,ch,postSTIM.eye == 3),1);
            BIbc = postCSD(:,ch,postSTIM.eye == 1) - nanmean(postCSD(1:pre,ch,postSTIM.eye == 1),1);
            ylab = 'naA/mm^3';
        case 3
            REbc = preSDF(:,ch,preSTIM.eye == 2) - nanmean(preSDF(1:pre,ch,preSTIM.eye == 2),1);
            LEbc = preSDF(:,ch,preSTIM.eye == 3) - nanmean(preSDF(1:pre,ch,preSTIM.eye == 3),1);
            BIbc = preSDF(:,ch,preSTIM.eye == 1) - nanmean(preSDF(1:pre,ch,preSTIM.eye == 1),1);
            ylab = 'spks/s';
        case 4
            REbc = postSDF(:,ch,postSTIM.eye == 2) - nanmean(postSDF(1:pre,ch,postSTIM.eye == 2),1);
            LEbc = postSDF(:,ch,postSTIM.eye == 3) - nanmean(postSDF(1:pre,ch,postSTIM.eye == 3),1);
            BIbc = postSDF(:,ch,postSTIM.eye == 1) - nanmean(postSDF(1:pre,ch,postSTIM.eye == 1),1);
            ylab = 'spks/s';
    end
    h(i) = subplot(1,4,i);
    plot(-pre:post,nanmean(REbc,3),'color',[1 0 0],'linewidth',2); hold on;
    plot(-pre:post,nanmean(LEbc,3),'color',[.5 0 .5],'linewidth',2); hold on;
    plot(-pre:post,nanmean(BIbc,3),'color',[0 0 0],'linewidth',2); hold on;
    set(gca,'tickdir','out','box','off','linewidth',2);
    xlim([-pre post]); ylabel(ylab); xlabel('t(ms)');
    
    refwin             = ismember(-pre:post,80:160);
    
    if i == 1 || i == 2
        csd_resp               = [csd_resp; squeeze(nanmean(REbc(refwin,:,:),1)); squeeze(nanmean(LEbc(refwin,:,:),1)); squeeze(nanmean(BIbc(refwin,:,:),1))]; %#ok<AGROW>
        
        if i == 1
            csd_group_drug = [csd_group_drug; zeros(size(REbc,3),1); zeros(size(LEbc,3),1); zeros(size(BIbc,3),1)]; %#ok<AGROW>
            csd_group_eye      = [csd_group_eye; repmat(2,size(REbc,3),1); repmat(3,size(LEbc,3),1); repmat(1,size(BIbc,3),1)];  %#ok<AGROW>
        else
            csd_group_drug = [csd_group_drug; ones(size(REbc,3),1); ones(size(LEbc,3),1); ones(size(BIbc,3),1)]; %#ok<AGROW>
            csd_group_eye      = [csd_group_eye; repmat(5,size(REbc,3),1); repmat(6,size(LEbc,3),1); repmat(4,size(BIbc,3),1)];  %#ok<AGROW>
        end
        
    end
    
    if i == 3 || i == 4
        sdf_resp               = [sdf_resp; squeeze(nanmean(REbc(refwin,:,:),1)); squeeze(nanmean(LEbc(refwin,:,:),1)); squeeze(nanmean(BIbc(refwin,:,:),1))]; %#ok<AGROW>
        
        if i == 3
            sdf_group_drug = [sdf_group_drug; zeros(size(REbc,3),1); zeros(size(LEbc,3),1); zeros(size(BIbc,3),1)]; %#ok<AGROW>
            sdf_group_eye      = [sdf_group_eye; repmat(2,size(REbc,3),1); repmat(3,size(LEbc,3),1); repmat(1,size(BIbc,3),1)];  %#ok<AGROW>
        else
            sdf_group_drug = [sdf_group_drug; ones(size(REbc,3),1); ones(size(LEbc,3),1); ones(size(BIbc,3),1)]; %#ok<AGROW>
            sdf_group_eye      = [sdf_group_eye; repmat(5,size(REbc,3),1); repmat(6,size(LEbc,3),1); repmat(4,size(BIbc,3),1)];  %#ok<AGROW>
        end
        
    end
    
end

preylm = get(h(1),'Ylim'); set(h(1),'YLim',preylm); set(h(2),'YLim',preylm); 
preylm = get(h(3),'Ylim'); set(h(3),'YLim',preylm); set(h(4),'YLim',preylm); 
for i = 1:4
    subplot(1,4,i)
    vline(0);
end

%%
% CSD
idx1   = csd_group_eye == 1 & csd_group_drug == 0; 
idx2   = csd_group_eye == 2 & csd_group_drug == 0; 
idx3   = csd_group_eye == 3 & csd_group_drug == 0; 

idx4   = csd_group_eye == 4 & csd_group_drug == 1; 
idx5   = csd_group_eye == 5 & csd_group_drug == 1; 
idx6   = csd_group_eye == 6 & csd_group_drug == 1; 
Y      = [nanmean(csd_resp(idx1)) nanmean(csd_resp(idx2)) nanmean(csd_resp(idx3)) ...
    nanmean(csd_resp(idx4)) nanmean(csd_resp(idx5)) nanmean(csd_resp(idx6))];
X      = 1:6;
E      =  [nanstd(csd_resp(idx1))/sqrt(length(idx1)) ...
    nanstd(csd_resp(idx2))/sqrt(length(idx2)) ....
    nanstd(csd_resp(idx3))/sqrt(length(idx3))  ...
    nanstd(csd_resp(idx4))/sqrt(length(idx4))  ...
    nanstd(csd_resp(idx5))/sqrt(length(idx5))  ...
    nanstd(csd_resp(idx6))/sqrt(length(idx6))];


P         = nan(length(X),length(X));
[~,pval]  = ttest2(csd_resp(csd_group_eye == 2),csd_resp(csd_group_eye == 5)); 
P(2,5)    = pval;
P(5,2)    = pval;
[~,pval]  = ttest2(csd_resp(csd_group_eye == 3),csd_resp(csd_group_eye == 6)); 
P(3,6)    = pval;
P(6,3)    = pval;
[~,pval]  = ttest2(csd_resp(csd_group_eye == 1),csd_resp(csd_group_eye == 4)); 
P(1,4)    = pval;
P(4,1)    = pval;




figure,set(gcf,'color','w','position',[1500 500 800 400]); 
superbar(X,abs(Y),'E',E,'P',P,'BarFaceColor',[0 0 0;1 0 0; 0.5 0 0.5;1 1 1; 1 1 1;1 1 1], ...
    'BarEdgeColor',[0 0 0;1 0 0; 0.5 0 0.5; 0 0 0;1 0 0; 0.5 0 0.5],'BarLineWidth',3,...
    'PLineBackingColor','none','PStarShowNS',false,'PStarColor',[1 0 1],'PLineColor',[1 0 1],'PStarFontSize',30); hold on;

set(gca,'tickdir','out','box','off');
set(gca,'linewidth',2);

P         = nan(length(X),length(X));
[~,pval]  = ttest2(csd_resp(csd_group_eye == 1),csd_resp(csd_group_eye == 2)); 
P(1,2)    = pval;
P(2,1)    = pval;
[~,pval]  = ttest2(csd_resp(csd_group_eye == 1),csd_resp(csd_group_eye == 3)); 
P(1,3)    = pval;
P(3,1)    = pval;
[~,pval]  = ttest2(csd_resp(csd_group_eye == 2),csd_resp(csd_group_eye == 3)); 
P(2,3)    = pval;
P(3,2)    = pval;

superbar(X,abs(Y),'E',E,'P',P,'BarFaceColor',[0 0 0;1 0 0; 0.5 0 0.5;1 1 1; 1 1 1;1 1 1], ...
    'BarEdgeColor',[0 0 0;1 0 0; 0.5 0 0.5; 0 0 0;1 0 0; 0.5 0 0.5],'BarLineWidth',3,...
    'PLineBackingColor','none','PStarColor',[0 0 1],'PLineColor',[0 0 1],'PStarShowNS',false,'PStarFontSize',30); hold on;


P         = nan(length(X),length(X));
[~,pval]  = ttest2(csd_resp(csd_group_eye == 4),csd_resp(csd_group_eye == 5)); 
P(4,5)    = pval;
P(5,4)    = pval;
[~,pval]  = ttest2(csd_resp(csd_group_eye == 4),csd_resp(csd_group_eye == 6)); 
P(4,6)    = pval;
P(6,4)    = pval;
[~,pval]  = ttest2(csd_resp(csd_group_eye == 5),csd_resp(csd_group_eye == 6)); 
P(5,6)    = pval;
P(6,5)    = pval;


superbar(X,abs(Y),'E',E,'P',P,'BarFaceColor',[0 0 0;1 0 0; 0.5 0 0.5;1 1 1; 1 1 1;1 1 1], ...
    'BarEdgeColor',[0 0 0;1 0 0; 0.5 0 0.5; 0 0 0;1 0 0; 0.5 0 0.5],'BarLineWidth',3,...
    'PLineBackingColor','none','PStarColor',[0 1 0],'PLineColor',[0 1 0],'PStarShowNS',false,'PStarFontSize',30); hold on;
ylabel('-1*nA/mm^3'); 

% MUA

[~,tbl,stats]  = anovan(csd_resp,{csd_group_eye});
% eye
eye_c = multcompare(stats,'alpha',.05,'ctype','bonferroni','dimension',1,'display','off');

idx1   = sdf_group_eye == 1 & sdf_group_drug == 0; 
idx2   = sdf_group_eye == 2 & sdf_group_drug == 0; 
idx3   = sdf_group_eye == 3 & sdf_group_drug == 0; 

idx4   = sdf_group_eye == 4 & sdf_group_drug == 1; 
idx5   = sdf_group_eye == 5 & sdf_group_drug == 1; 
idx6   = sdf_group_eye == 6 & sdf_group_drug == 1; 
Y      = [nanmean(sdf_resp(idx1)) nanmean(sdf_resp(idx2)) nanmean(sdf_resp(idx3)) ...
    nanmean(sdf_resp(idx4)) nanmean(sdf_resp(idx5)) nanmean(sdf_resp(idx6))];
X      = 1:6;
E      =  [nanstd(sdf_resp(idx1))/sqrt(length(idx1)) ...
    nanstd(sdf_resp(idx2))/sqrt(length(idx2)) ....
    nanstd(sdf_resp(idx3))/sqrt(length(idx3))  ...
    nanstd(sdf_resp(idx4))/sqrt(length(idx4))  ...
    nanstd(sdf_resp(idx5))/sqrt(length(idx5))  ...
    nanstd(sdf_resp(idx6))/sqrt(length(idx6))];


P         = nan(length(X),length(X));
[~,pval]  = ttest2(sdf_resp(sdf_group_eye == 2),sdf_resp(sdf_group_eye == 5)); 
P(2,5)    = pval;
P(5,2)    = pval;
[~,pval]  = ttest2(sdf_resp(sdf_group_eye == 3),sdf_resp(sdf_group_eye == 6)); 
P(3,6)    = pval;
P(6,3)    = pval;
[~,pval]  = ttest2(sdf_resp(sdf_group_eye == 1),sdf_resp(sdf_group_eye == 4)); 
P(1,4)    = pval;
P(4,1)    = pval;


figure,set(gcf,'color','w','position',[1500 500 800 400]); 
superbar(X,abs(Y),'E',E,'P',P,'BarFaceColor',[0 0 0;1 0 0; 0.5 0 0.5;1 1 1; 1 1 1;1 1 1], ...
    'BarEdgeColor',[0 0 0;1 0 0; 0.5 0 0.5; 0 0 0;1 0 0; 0.5 0 0.5],'BarLineWidth',3,...
    'PLineBackingColor','none','PStarShowNS',false,'PStarColor',[1 0 1],'PLineColor',[1 0 1],'PStarFontSize',30); hold on;

set(gca,'tickdir','out','box','off');
set(gca,'linewidth',2);

P         = nan(length(X),length(X));
[~,pval]  = ttest2(sdf_resp(sdf_group_eye == 1),sdf_resp(sdf_group_eye == 2)); 
P(1,2)    = pval;
P(2,1)    = pval;
[~,pval]  = ttest2(sdf_resp(sdf_group_eye == 1),sdf_resp(sdf_group_eye == 3)); 
P(1,3)    = pval;
P(3,1)    = pval;
[~,pval]  = ttest2(sdf_resp(sdf_group_eye == 2),sdf_resp(sdf_group_eye == 3)); 
P(2,3)    = pval;
P(3,2)    = pval;

superbar(X,abs(Y),'E',E,'P',P,'BarFaceColor',[0 0 0;1 0 0; 0.5 0 0.5;1 1 1; 1 1 1;1 1 1], ...
    'BarEdgeColor',[0 0 0;1 0 0; 0.5 0 0.5; 0 0 0;1 0 0; 0.5 0 0.5],'BarLineWidth',3,...
    'PLineBackingColor','none','PStarColor',[0 0 1],'PLineColor',[0 0 1],'PStarShowNS',false,'PStarFontSize',30); hold on;


P         = nan(length(X),length(X));
[~,pval]  = ttest2(sdf_resp(sdf_group_eye == 4),sdf_resp(sdf_group_eye == 5)); 
P(4,5)    = pval;
P(5,4)    = pval;
[~,pval]  = ttest2(sdf_resp(sdf_group_eye == 4),sdf_resp(sdf_group_eye == 6)); 
P(4,6)    = pval;
P(6,4)    = pval;
[~,pval]  = ttest2(sdf_resp(sdf_group_eye == 5),sdf_resp(sdf_group_eye == 6)); 
P(5,6)    = pval;
P(6,5)    = pval;


superbar(X,abs(Y),'E',E,'P',P,'BarFaceColor',[0 0 0;1 0 0; 0.5 0 0.5;1 1 1; 1 1 1;1 1 1], ...
    'BarEdgeColor',[0 0 0;1 0 0; 0.5 0 0.5; 0 0 0;1 0 0; 0.5 0 0.5],'BarLineWidth',3,...
    'PLineBackingColor','none','PStarColor',[0 1 0],'PLineColor',[0 1 0],'PStarShowNS',false,'PStarFontSize',30); hold on;
ylabel('spks/s'); 
%%

for i = 1:2
    switch i
        case 1
            REbc = preCSD(:,ch,preSTIM.eye == 2) - nanmean(preCSD(1:pre,ch,preSTIM.eye == 2),1);
            LEbc = preCSD(:,ch,preSTIM.eye == 3) - nanmean(preCSD(1:pre,ch,preSTIM.eye == 3),1);
            BIbc = preCSD(:,ch,preSTIM.eye == 1) - nanmean(preCSD(1:pre,ch,preSTIM.eye == 1),1);
            ylab = 'naA/mm^3';
        case 2
            REbc = postCSD(:,ch,postSTIM.eye == 2) - nanmean(postCSD(1:pre,ch,postSTIM.eye == 2),1);
            LEbc = postCSD(:,ch,postSTIM.eye == 3) - nanmean(postCSD(1:pre,ch,postSTIM.eye == 3),1);
            BIbc = postCSD(:,ch,postSTIM.eye == 1) - nanmean(postCSD(1:pre,ch,postSTIM.eye == 1),1);
            ylab = 'naA/mm^3';
        case 3
            REbc = preSDF(:,ch,preSTIM.eye == 2) - nanmean(preSDF(1:pre,ch,preSTIM.eye == 2),1);
            LEbc = preSDF(:,ch,preSTIM.eye == 3) - nanmean(preSDF(1:pre,ch,preSTIM.eye == 3),1);
            BIbc = preSDF(:,ch,preSTIM.eye == 1) - nanmean(preSDF(1:pre,ch,preSTIM.eye == 1),1);
            ylab = 'spks/s';
        case 4
            REbc = postSDF(:,ch,postSTIM.eye == 2) - nanmean(postSDF(1:pre,ch,postSTIM.eye == 2),1);
            LEbc = postSDF(:,ch,postSTIM.eye == 3) - nanmean(postSDF(1:pre,ch,postSTIM.eye == 3),1);
            BIbc = postSDF(:,ch,postSTIM.eye == 1) - nanmean(postSDF(1:pre,ch,postSTIM.eye == 1),1);
            ylab = 'spks/s';
    end
    
    if i == 1
        figure,set(gcf,'color','w'); 
        subplot(1,2,1)
        plot(-pre:post,nanmean(BIbc,3),'color','k','linewidth',2); hold on;
        plot(-pre:post,nanmean(REbc,3),'color',[1 0 0],'linewidth',2); hold on; 
        plot(-pre:post,nanmean(LEbc,3),'color',[0.5 0 0.5],'linewidth',2); hold on; 
       
        re_pre = nanmean(REbc,3);  le_pre = nanmean(LEbc,3); bi_pre = nanmean(BIbc,3); 
        plot(-pre:post,-1*sqrt((nanmean(REbc,3).^2) + nanmean(LEbc,3).^2),'color','b','linestyle','--','linewidth',2); hold on; 
        plot(-pre:post,(nanmean(REbc,3) + nanmean(LEbc,3)),'color','g','linestyle','--','linewidth',2); hold on; 
          set(gca,'tickdir','out','box','off','linewidth',2); vline(0); xlim([-pre post]); ylim([-6000 1000]); 
    else
        
        subplot(1,2,2)
        plot(-pre:post,nanmean(BIbc,3),'color','k','linewidth',2); hold on;
        plot(-pre:post,nanmean(REbc,3),'color',[1 0 0],'linewidth',2); hold on; 
        plot(-pre:post,nanmean(LEbc,3),'color',[0.5 0 0.5],'linewidth',2); hold on; 
        plot(-pre:post,bi_pre,'color','k','linewidth',3); hold on;
        
        plot(-pre:post, nanmean(LEbc,3) + nanmean(REbc,3) + (re_pre - nanmean(REbc,3)),'color','b','linestyle','--','linewidth',2); hold on; 

      set(gca,'tickdir','out','box','off','linewidth',2);  xlim([-pre post]); ylim([-6000 1000]); vline(0);
    end


    
end

%% RFORI TUNING?

ch = chans ==22; 
oris = nanunique(preSTIM.tilt);
for i = 1:4
    for o = 1:length(oris)
        switch i
            case 1
                STIM = preSTIM;
                data = preCSD;
                preoriCSD{o} =data(:,ch,STIM.tilt == oris(o)) - nanmean(data(1:pre,ch,STIM.tilt == oris(o)),1);
            case 2
                STIM = postSTIM;
                data = postCSD;
                postoriCSD{o} =data(:,ch,STIM.tilt == oris(o)) - nanmean(data(1:pre,ch,STIM.tilt == oris(o)),1);
            case 3
                STIM = preSTIM;
                data = preSDF;
                preoriSDF{o} =data(:,ch,STIM.tilt == oris(o)) - nanmean(data(1:pre,ch,STIM.tilt == oris(o)),1);
            case 4
                STIM = postSTIM;
                data = postSDF;
                postoriSDF{o} =data(:,ch,STIM.tilt == oris(o)) - nanmean(data(1:pre,ch,STIM.tilt == oris(o)),1);
        end
    end
end


figure,set(gcf,'color','w'); 
clear resp feature x y f
for i = 1:4
    clear data
    switch i
        case 1
             STIM = preSTIM;
            data = preoriCSD;
        case 2
             STIM = postSTIM;
            data = postoriCSD;
        case 3
            STIM = preSTIM;
            data = preoriSDF;
        case 4
            STIM = postSTIM;
            data = postoriSDF;
    end
    h(i) = subplot(2,2,i);
    refwin = ismember(-pre:post,50:250);
    mns    = nan(length(oris),1);   sem   = nan(length(oris),1); 
    resp = []; group = []; 
    for o = 1:length(oris)
        mns(o)  = nanmean(squeeze(nanmean(data{o}(refwin,:,:),3))); 
        sem(o)  = nanstd(squeeze(nanmean(data{o}(refwin,:,:),3)))./sqrt(size(data{o},3)); 
        resp    = [resp; squeeze(nanmean(data{o}(refwin,:,:),1))];
        group   = [group; repmat(oris(o),size(data{o},3),1)]; 
    end
    
     clear tbl stats
    [~,tbl] = anovan(resp,{group},'display','off'); 
 
 
    plot(oris,mns,'-ko','linewidth',2); hold on; 
    errorbar(oris,mns,sem,'linestyle','none'); hold on; 
    set(gca,'box','off','tickdir','out','linewidth',2); 
    title(gca,sprintf('p = %0.2f',tbl{2,7})); hold on; 
    xlim([0 180]); 
end
preylm = get(h(1),'Ylim'); postylm = get(h(2),'Ylim'); 
set(h(1),'YLim',[min([preylm postylm]) max([preylm postylm])]); 
set(h(2),'YLim',[min([preylm postylm]) max([preylm postylm])]);

preylm = get(h(3),'Ylim'); postylm = get(h(4),'Ylim'); 
set(h(3),'YLim',[min([preylm postylm]) max([preylm postylm])]); 
set(h(4),'YLim',[min([preylm postylm]) max([preylm postylm])]);

h(3).XLabel.String = 'ori (deg)'; 
h(3).YLabel.String = 'spks/s'; 
h(1).YLabel.String = 'nA/mm^3'; 
h(4).XLabel.String = 'ori (deg)'; 
 

%%

for i = 1
    clear data
    switch i
        case 1
             STIM = preSTIM;
            data = preoriCSD;
        case 2
             STIM = postSTIM;
            data = postoriCSD;
        case 3
            STIM = preSTIM;
            data = preoriSDF;
        case 4
            STIM = postSTIM;
            data = postoriSDF;
    end

    figure, set(gcf,'color','w','position',[50 50 800 800])
    for o = 1:length(oris)
        subplot(4,4,o)
        plot(-100:200,nanmedian(squeeze(data{o}),2),'color',[0.4 0.4 0.4],'linewidth',2); hold on; 
        set(gca,'box','off','tickdir','out'); 
    end
  
    

end
 


