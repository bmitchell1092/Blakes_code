function  makealignCSDfig(preDAT,postDAT,chans,fullch,depths,presink,postsink,PRE,POST)


%% calculate CSD
chans = [1:24];
for tr = 1:size(preDAT,3)
    precsd(:,:,tr) = mod_iCSD(preDAT(:,chans,tr)');
end

for tr = 1:size(postDAT,3)
    postcsd(:,:,tr) = mod_iCSD(postDAT(:,chans,tr)');
end

%depths  = depths(chans); 
totchan = length(chans); 

%% baseline correct, pad array and filter CSD
h_trcsd = bsxfun(@minus,precsd,mean(precsd(:,1:abs(PRE),:),2)); bs_h_pretrcsd = h_trcsd; 
h_CSD = padarray(h_trcsd,[1 0],NaN,'replicate');
preCSDf = filterCSD(mean(h_CSD,3)); lm = max(max(abs(preCSDf))).*.7;

clear h_trcsd h_CSD 
h_trcsd = bsxfun(@minus,postcsd,mean(postcsd(:,1:abs(PRE),:),2)); bs_h_posttrcsd = h_trcsd;
h_CSD = padarray(h_trcsd,[1 0],NaN,'replicate');
postCSDf = filterCSD(mean(h_CSD,3)); lm = max(max(abs(preCSDf)));

%% enter filtered CSDs in larger arrays 
emparray = nan(size(fullch,1)*10,size(postCSDf,2),2); 

start = find(fullch(:,1) == chans(1)).*10-5; prestart = start; 
emparray([start:(start + size(preCSDf,1)-1)],:,1) = preCSDf; 
preidx = (find(ismember(fullch(:,1),chans)).*10) -5; 

clear start;
start = find(fullch(:,2) == chans(1)).*10; poststart = start; 
emparray([start:(start + size(postCSDf,1)-1)],:,2) = postCSDf; 
postidx = (find(ismember(fullch(:,2),chans)).*10) -5;  

%% plot figures in reference to channel depth
figure,set(gcf,'Color','w','PaperPositionMode','auto','Position',[1 1 600 1000])

subplot(1,2,1)
h = imagesc([PRE:POST],1:size(emparray,1),emparray(:,:,1)); 
set(gca,'YTick',preidx(2:4:end-1),'YTickLabel',chans(2:4:end-1),'FontSize',20);
set(h,'alphadata',~isnan(emparray(:,:,1)));
caxis([-lm lm]);
vline(0);colormap(flipud(colormap('jet')));colorbar;
set(gca,'TickDir','out','Box','off'); 
xlabel('t (ms)'); ylabel('channel'); 
h = vline(0); set(h,'LineWidth',2,'Color','k');
 

subplot(1,2,2)
h = imagesc([PRE:POST],1:size(emparray,1),emparray(:,:,2)); 
set(gca,'YTick',postidx(2:4:end-1),'YTickLabel',chans(2:4:end-1),'Box','off','FontSize',20);
set(h,'alphadata',~isnan(emparray(:,:,2)));
caxis([-lm lm]);
vline(0);colormap(flipud(colormap('jet')));colorbar;
set(gca,'TickDir','out'); 
xlabel('t (ms)'); ylabel('channel'); 
h = vline(0); set(h,'LineWidth',2,'Color','k');

%% plot figures in reference to cortical depth 
% this time add channels into matrix so that the two files align
clear emparray h_CSD premat postmat

emparray   = nan(size(fullch,1),size(postCSDf,2),2); 
sinkchan   = size(emparray,1)./2; 

prechans = fullch(find(sum(~isnan(fullch),2) == 2),1);
prechans = prechans(prechans>1 & prechans <= size(bs_h_pretrcsd,1)); 
emparray([((sinkchan - presink + 1 ) +1)   : ((sinkchan +  (length(depths) - presink))-1)],:,1)   = nanmean(bs_h_pretrcsd(prechans,:,:),3);

postchans = fullch(find(sum(~isnan(fullch),2) == 2),2);
postchans = postchans(postchans>1 & postchans <= size(bs_h_posttrcsd,1)); 
emparray([((sinkchan - postsink + 1 )+1)   : ((sinkchan +  (length(depths) - postsink))-1)],:,2)  = nanmean(bs_h_posttrcsd(postchans(1:end-1),:,:),3);

premat     = emparray([(sinkchan - presink + 1)+1  : (sinkchan + (length(prechans) - presink))-1],:,1); 
postmat    = emparray([(sinkchan - postsink + 1)+1 : (sinkchan + (length(prechans) - postsink))-1],:,2);

h_CSD      = padarray(premat,[1 0],NaN,'replicate');
A_preCSDf  = filterCSD(mean(h_CSD,3)); 

clear h_CSD
h_CSD      = padarray(postmat,[1 0],NaN,'replicate');
A_postCSDf = filterCSD(mean(h_CSD,3)); 
%%
figure,set(gcf,'Color','w','PaperPositionMode','auto','Position',[1 1 600 1000])
subplot(1,2,1)
h = imagesc([PRE:POST],depths(2:end-1),A_preCSDf); 
caxis([-lm lm]);
vline(0);colormap(flipud(colormap('jet')));colorbar;
set(gca,'TickDir','out','Box','off','YDir','normal'); 
xlabel('t (ms)'); 
h = vline(0); set(h,'LineWidth',2,'Color','k');
ylabel('cortical depth (mm)'); 

subplot(1,2,2)
h = imagesc([PRE:POST],depths(2:end-1),A_postCSDf); 
caxis([-lm lm]);
vline(0);colormap(flipud(colormap('jet')));colorbar;
set(gca,'TickDir','out','Box','off','YDir','normal'); 
xlabel('t (ms)'); 
h = vline(0); set(h,'LineWidth',2,'Color','k');
ylabel('cortical depth (mm)'); 

%%
end