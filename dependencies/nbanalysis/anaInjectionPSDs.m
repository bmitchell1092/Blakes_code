clear 
addpath(genpath('/users/kaciedougherty/documents/code/nbanalysis/')); 
addpath(genpath('/users/kaciedougherty/documents/code/fNPMK')); 
addpath('/Users/kaciedougherty/Documents/Code/IRASA');
addpath(genpath('/users/kaciedougherty/documents/code/lgn-dichoptic'));

preBRdatafile   = '190606_B_rforidrft001';
postBRdatafile  = '190606_B_rforidrft002';
drname          = ['/volumes/bigdrobo2/drobo2/data/rig022/' preBRdatafile(1:8) '/'];

el              = 'eA'; %'eA';
pre             = 50;
post            = 1000;
chans           = 10:32; 

%% get baseline and post injection data 

clear BRdatafile; 
BRdatafile                                          = [drname preBRdatafile];
[preLFP,preSDF,presdf,prelfp,preCSD,~,preSTIM]      = getInjectionData(BRdatafile,el,chans,pre,post); 

clear BRdatafile; 
BRdatafile                                             = [drname postBRdatafile];
[postLFP,postSDF,postsdf,postlfp,postCSD,TTL,postSTIM] = getInjectionData(BRdatafile,el,chans,pre,post); 

%% remove trials when spiking comes back after injection 

rmtrials   = []; 
if ~isempty(rmtrials)
h_postSDF  = postSDF; 
h_postCSD  = postCSD; 
h_postSTIM = postSTIM; 
h_postLFP  = postLFP ;

postSDF(:,:,rmtrials) = []; 
postCSD(:,:,rmtrials) = []; 
postLFP(:,:,rmtrials) = []; 
fds = fields(postSTIM); 
for f = 1:length(fds)
postSTIM.(fds{f})(:,:,rmtrials) = []; 
end
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

%% FFTs 

for ch = 1:length(chans)
    
L4pre(ch)       = getFractalOscillatoryComponents(squeeze(preLFP(:,ch,:)));

L4post(ch)      = getFractalOscillatoryComponents(squeeze(postLFP(:,ch,:))); 

end
%%
figure, set(gcf,'color','w','position',[1 1 300 800]); 
rows = ceil(length(chans)./2); cols = 2; 
for ch = 1:length(chans)
    
    subplot(rows,cols,ch)
    
    plot(L4pre(ch).freq,nanmean(L4pre(ch).osci,2),'color','r','linewidth',2); hold on;
    plot(L4post(ch).freq,nanmean(L4post(ch).osci,2),'color',[.5 0 .5],'linewidth',2); hold on;
    
    plot(L4pre(ch).freq,nanmean(L4pre(ch).osci,2) + (nanstd(L4pre(ch).osci,0,2)./sqrt(size(L4pre(ch).osci,2))),'color','r','linewidth',1); hold on;
    plot(L4pre(ch).freq,nanmean(L4pre(ch).osci,2) - (nanstd(L4pre(ch).osci,0,2)./sqrt(size(L4pre(ch).osci,2))),'color','r','linewidth',1); hold on;
    
    plot(L4post(ch).freq,nanmean(L4post(ch).osci,2) + (nanstd(L4post(ch).osci,0,2)./sqrt(size(L4post(ch).osci,2))),'color',[.5 0 .5],'linewidth',1); hold on;
    plot(L4post(ch).freq,nanmean(L4post(ch).osci,2) - (nanstd(L4post(ch).osci,0,2)./sqrt(size(L4post(ch).osci,2))),'color',[.5 0 .5],'linewidth',1); hold on;
    
    xlim([2 30]);
    set(gca,'box','off','tickdir','out','fontsize',14); 
    
    ylms(:,ch) = get(gca,'YLim'); 
    
end
L4sink  = 25; 
sinkid  = find(chans == L4sink); 
depths  =  ((sinkid-1)/10):-0.1: 0 : -0.1: (-1*(length(chans) - sinkid)/10);
for ch = 1:length(chans)
      subplot(rows,cols,ch)
      ylim([min(min(ylms)) max(max(ylms))]); 
      title(gca,num2str(depths(ch))); 
end
ylabel('power (a.u.)'); xlabel('freq (Hz)'); 






