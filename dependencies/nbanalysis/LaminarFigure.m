%% RF MAPPING
clear
addpath(genpath('/Users/kaciedougherty/Documents/Code/v1-dichoptic-collab'))
addpath(genpath('/users/kaciedougherty/documents/code/nbanalysis'));

plotdr = '/users/kaciedougherty/documents/plotsfigures/decoding_dotmap';
brdrname = sprintf('/users/kaciedougherty/documents/neurophysdata/');


%%
BRdatafile   = '161004_E_dotmapping003'; load(strcat('/Users/kaciedougherty/Documents/neurophysdata/',BRdatafile,'.ppnev'),'-MAT');

NEV          = ppNEV; clear ppNEV;
el           = 'eD';
el_array     = 1:1:24;
badobs       = [];
flag_RFM     = true;    % newer method, a bit slower [set false if line below is true]
flag_2Dfit   = false;   % prior method [set false if line above is true]
pre          = -20;
post         = 150;
tvec         = pre:post;

if isempty(el_array), nel = 1; 
else nel = length(el_array);
end

 
for i = 1:nel

    clear SPK Fs I r_sua x y d eye sdf STIM sdftr X Y Z iZ uZ
    
    e              = el_array(i);
    elabel         = sprintf('%s%02u',el,e);
    eidx           = find(cell2mat(cellfun(@(x) ~isempty(strfind(x',elabel)),{NEV.ElectrodesInfo.ElectrodeLabel},'UniformOutput',0)));
    I              =  NEV.Data.Spikes.Electrode == eidx;
    SPK            = double(NEV.Data.Spikes.TimeStamp(I));
    Fs             = double(NEV.MetaTags.SampleRes);
    ct             = 0;
    STIM           = getEventTimeInfo(brdrname,BRdatafile,0);

    [sdf, sua, tm] = spk2sdf(SPK,Fs);
    sdftr          = squeeze(trigData(sdf',floor(STIM.onsets./30),-pre,post));
    sdftr          = sdftr - repmat(nanmean(sdftr(tvec < 0,:),1),[length(tvec) 1]);

    rfcomputation  = 'halfmax';
    res            = 0.05; % dva per pix in matrix
    dd             = max(STIM.diameter);
    [X,Y]          = meshgrid(min(STIM.dot_x)-dd:res:max(STIM.dot_x)+dd, min(STIM.dot_y)-dd:res:max(STIM.dot_y)+dd);
    Z              = NaN(size(X,1),size(X,2),length(STIM.dot_x));
    r_sua          = nanmean(sdftr(tvec > 30,:),1);
    
    for obs = 1:length(STIM.dot_x)
        xx   = STIM.dot_x(obs);
        yy   = STIM.dot_y(obs);
        dd   = STIM.diameter(obs);
        fill = sqrt(abs(X-xx).^2 + abs(Y-yy).^2) < dd/2;
        if ~any(any(fill)), error('check matrix'), end
        trldat       = Z(:,:,obs);
        trldat(fill) = r_sua(obs);
        Z(:,:,obs)   = trldat;
    end
    
    uZ        = nanmean(Z,3);
    
    % determine boundary of RF, using half max OR t-test
    baseline  = nanmedian(r_sua);
    maxspking = max(max(uZ));
    switch rfcomputation
        case 'ttest'
            IslandMap = ttest(Z,baseline,'dim',3,'alpha',0.001,'tail','right');
            IslandMap(isnan(IslandMap)) = 0;
        case 'halfmax'
            IslandMap = uZ > maxspking*0.5;
    end
    
    % use imagetoobox, same as PNAS
    CC = bwconncomp(IslandMap);
    [~, index] = max(cellfun('size',CC.PixelIdxList,1)); %the number of squares that surpass criteria and arein the largest connected patch
    STATS = regionprops(CC,'BoundingBox','Centroid');
    width = (STATS(index).BoundingBox(3:4) .* res) ;
    centroid = round(STATS(index).Centroid);
    centroid = [X(1,centroid(1)) Y(centroid(2),1)];
    rfboundary = [centroid(1)-width(1)/2,centroid(2)-width(2)/2, width(1), width(2)];
    
    uZ  = (uZ - baseline)./(maxspking - baseline);
    
    
    %  smooth map.
    sigma = unique(STIM.diameter)/2/3 /res;
    if sigma > 0
        iNaN = isnan(uZ);
        uZ(iNaN) = 0;
        iZ = imgaussfilt(uZ,sigma);
        iZ(iNaN) = NaN;
    else
        iZ = uZ;
    end
    
    %subplot(length(el_array),1,i)
    figure, set(gcf,'color','w','position',[1 1 200 200]); 
    b       = imagesc(X(1,:),Y(:,1),iZ); hold on
    r_sua   = rectangle('Position',rfboundary,'Curvature',[1,1]); set(b,'AlphaData',~isnan(iZ));
    
    %adjust map to accomidate NaNs
    colormap jet
    map   = colormap;
    cstep = range(get(gca,'CLim')) / length(map);
    set(gca,'CLim', get(gca,'CLim') + [-1*cstep 0]);
    colormap(map);
    title(gca,num2str(e),'fontsize',10);
    set(gca,'Ydir','normal','Box', 'off','fontsize',10,'tickdir','out','ticklength',[.03 .03])
    grid off
    caxis([0 1]);
    colorbar; 
    saveas(gcf,strcat(plotdr,'/','dot_',num2str(el_array(i))),'eps2c');

end

%% LFP and CSD from RFORI FILE 
clearvars -except brdrname plotdr el_array el
BRdatafile   = '170711_I_rfori005';
el = 'eD'; 
el_array     = 7:28;
pre          = -50;
post         = 400;
STIM         = getEventTimeInfo(brdrname,BRdatafile,0); 
lfp          = getLFP([brdrname '/' BRdatafile],'ns2','eD','ascending'); 
lfp          = lfp(:,el_array); 
[DAT, TM]    = trigData(lfp, floor(STIM.onsets./30) , -pre, post); 


clear csd bcsd fcsd
csd = nan(length(pre:post),size(DAT,2)-2,size(DAT,3)); 
for tr = 1:size(DAT,3)
csd(:,:,tr) = mod_iCSD(DAT(:,:,tr)')';
end
mcsd = nanmean(csd,3); 
bcsd = mcsd - repmat(nanmean(mcsd(1:abs(pre),:),1),[size(mcsd,1) 1]); 
fcsd = filterCSD(bcsd')';  

load(['/volumes/Toshiba External/' BRdatafile '.ppnev'],'-MAT');
Fs             = double(ppNEV.MetaTags.SampleRes);
clear sdftr spk_bin
for i = 1:length(el_array)
    clear e elabel eidx I SPK sdf
    e              = el_array(i);
    elabel         = sprintf('%s%02u',el,e);
    eidx           = find(cell2mat(cellfun(@(x) ~isempty(strfind(x',elabel)),{ppNEV.ElectrodesInfo.ElectrodeLabel},'UniformOutput',0)));
    I              = ppNEV.Data.Spikes.Electrode == eidx;
    SPK            = double(ppNEV.Data.Spikes.TimeStamp(I));
    sdf            = spk2sdf(SPK,Fs); 
    sdftr(:,:,i)   = squeeze(trigData(sdf',floor(STIM.onsets./30),-pre,post));
    spk_bin(:,:,i) = squeeze(trigBinaryData(SPK,pre,post,floor(STIM.onsets./30)));
end
sdftr   = permute(sdftr,[1 3 2]); 
spk_bin = permute(sdftr,[1 3 2]); 

mn_sdf       = nanmean(sdftr,3); 
cols       = 4; 
subplts    = 1:length(el_array)*cols; 
first_col  = subplts(1:cols:length(subplts)); 
second_col = subplts(2:cols:length(subplts));
third_col  = subplts(3:cols:length(subplts));
fourth_col = subplts(4:cols:length(subplts));

figure,set(gcf,'color','w','position',[1 1 1600 800]); 
msdf = mn_sdf - repmat(mean(mn_sdf(1:abs(pre),:),1),[length(pre:post) 1]); 
ylm = [min(min(msdf)) max(max(msdf))]; 

for ch = 1:size(msdf,2)
    h = subplot(length(first_col),cols,first_col(ch));
    plot(TM,msdf(:,ch),'color','k','linewidth',2);
    set(gca,'tickdir','out','box','off','xtick',[],'ytick',[]);
    set(gca,'ycolor','w'); 
    sp = get(h,'Position'); 
    
    ylim(ylm); xlim([pre post]); 
    vline(0); 
    set(h,'Position',sp.*[1 1 1.3 1.3]); 
end

mn_lfp = nanmean(DAT,3);
mlfp = mn_lfp - repmat(mean(mn_lfp(1:abs(pre),:),1),[length(pre:post) 1]); 
ylm = [min(min(mlfp)) max(max(mlfp))]; 

for ch = 1:size(mlfp,2)
    h = subplot(length(first_col),cols,second_col(ch));
    plot(TM,mlfp(:,ch),'color','k','linewidth',2);
    set(gca,'ycolor','w'); 
    set(gca,'tickdir','out','box','off','xtick',[],'ytick',[]);
    sp = get(h,'Position'); 
    set(h,'Position',sp.*[1 1 1.3 1.3]); 
    ylim(ylm); xlim([pre post]); 
    vline(0); 
end

ylm = [min(min(bcsd))+10 max(max(bcsd))-10]; 
for ch = 2:size(mlfp,2)-1
    h = subplot(length(first_col),cols,third_col(ch));
    plot(TM,bcsd(:,ch-1),'color','k','linewidth',2);
    set(gca,'ycolor','w'); 
    set(gca,'tickdir','out','box','off','xtick',[],'ytick',[]);
    sp = get(h,'Position'); 
    set(h,'Position',sp.*[1 1 1.3 1.3]); 
      ylim(ylm); xlim([pre post]); 
    vline(0); 
end


subplot(length(first_col),cols,fourth_col(2:end-1));
imagesc(TM,[],fcsd'); 
colormap(flipud(colormap('jet'))); %colorbar; 
caxis([-600 600]); 
set(gca,'box','off','tickdir','out'); 
colorbar















%%

  
          
                       
            % plot CSD and MUA for each electrode:
            pre  = -200;
            post = 300;
            time = pre:post; 
            for tr = 1:length(onsetT)
                emuaA(:,:,tr) = muaA(onsetT(tr)+pre:onsetT(tr)+post,:);
                emuaB(:,:,tr) = muaB(onsetT(tr)+pre:onsetT(tr)+post,:);
                elfpA(:,:,tr) = lfpA(onsetT(tr)+pre:onsetT(tr)+post,:);
                elfpB(:,:,tr) = lfpB(onsetT(tr)+pre:onsetT(tr)+post,:);
                
            end
            for tr = 1:length(onsetT)
                csdA(:,:,tr) = mod_iCSD(elfpA(:,:,tr)')';
                csdB(:,:,tr) = mod_iCSD(elfpB(:,:,tr)')';
            end
            
     
            mcsdA = nanmean(csdA,3);
            mcsdB = nanmean(csdB,3);
            mmuaA = nanmean(emuaA,3);
            mmuaB = nanmean(emuaB,3);
            
            figure
            set(gcf,'Color',[1 1 1],'Position',[1 1 700 800],'PaperPositionMode','auto');
            
            whichElectrode = el{ID(1)}; 
            % MUA
            subplot(1,2,[1])
            
            mmultA = (mmuaA - repmat(mean(mmuaA,1),size(mmuaA,1),1))./(repmat(mean(mmuaA,1),size(mmuaA,1),1));
         
            cst = .2;
            [rsmuaA] = rsEqChans(mmultA,cst);
            for chan = 1:size(rsmuaA,2)
                
                plot(time,rsmuaA(:,chan),'Color','k','LineWidth',2);
                axis tight;
                h = text(pre-80,rsmuaA(1,chan),sprintf('%02d',chan)); set(h,'Color','k','FontSize',14);
                
                hold on
            end
            axis tight
            set(gca,'FontSize',14,'YTickLabel',[],'YTick',[],'Box','off','TickDir','out');
            xlabel('t (ms) from flash');
        
            t = title(gca,sprintf('%s: MUA %s',filesforday{fd}(1:end-4),whichElectrode),'interpreter','none');
            set(t,'FontWeight','bold');
            v = vline(0);
            set(v,'Color',[0.5 0.5 0.5],'LineWidth',1);
            
            %CSD
            
            subplot(1,2,[2])
           
            [fcsd] = filterCSD(mcsdA');
            nn = nan(10,size(fcsd,2));
            pd_fcsd = [nn ;fcsd; nn]; % pad with nans for lost channels
            was = 1:size(mmuaA,2); % n contacts
          
            olddepths = was;
            ticks     = linspace(1,size(pd_fcsd,1),length(was));
            
            pr = 1;  %set lims for colorbar
            mn = ((min(min(fcsd))).*pr) ;
            mx = ((max(max(fcsd))).*pr) ;
            
            h = imagesc(time,1:size(pd_fcsd,1),pd_fcsd,[mn -mn]); cb = colorbar;
            set(h,'AlphaData',~isnan(pd_fcsd)); 
            cm =  size(colormap('jet'),1);
            
            map = flipud(colormap('jet'));
            colormap(map); 
            ylim(cb,[mn -mn]);  %set lims for colorbar and add in white for nans
            set(cb,'YDir','reverse');
            set(gca,'YTick',ticks(2:2:end-1),'YTicklabel',olddepths(2:2:end-1),'FontSize',14); set(gca,'TickDir','out');
            xlabel('t (ms) from flash'); ylabel('channel number'); 
            set(get(cb,'ylabel'),'string','nA/mm^3','FontSize',14); 
            
            t = title(gca,sprintf('%s CSD %s',filesforday{fd}(1:end-4),whichElectrode),'interpreter','none');
            set(t,'FontWeight','bold');
            v = vline(0); set(v,'Color',[0.5 0.5 0.5],'LineWidth',1);
            dur = length(MUA)./1000/60; 
            h = text(post+80,csd(1,1)-15,sprintf('n = %d trs, dur %d min',length(onsetT),dur)); set(h,'Color','k','FontSize',14);

            % 2nd electrode: 
            figure
            set(gcf,'Color',[1 1 1],'Position',[1 1 700 800],'PaperPositionMode','auto');
            
            whichElectrode = el{ID(2)}; 
            % MUA
            subplot(1,2,[1])
            
            mmultB = mmuaB - repmat(mean(mmuaB,1),size(mmuaB,1),1)./(repmat(mean(mmuaB,1),size(mmuaB,1),1));
         
            cst = .2;
            [rsmuaB] = rsEqChans(mmultB,cst);
            for chan = 1:size(rsmuaB,2)
                
                plot(time,rsmuaB(:,chan),'Color','k','LineWidth',2);
                axis tight;
                h = text(pre-80,rsmuaB(1,chan),sprintf('%02d',chan)); set(h,'Color','k','FontSize',14);
                
                hold on
            end
            axis tight
            set(gca,'FontSize',14,'YTickLabel',[],'YTick',[],'Box','off','TickDir','out');
            xlabel('t (ms) from flash');  
            t = title(gca,sprintf('%s: MUA %s',filesforday{fd}(1:end-4),whichElectrode),'interpreter','none');
            set(t,'FontWeight','bold');
            v = vline(0);
            set(v,'Color',[0.5 0.5 0.5],'LineWidth',1);
            
            %CSD
            
            subplot(1,2,[2])
           clear fcsd 
            [fcsd] = filterCSD(mcsdB');
            nn = nan(10,size(fcsd,2));
            pd_fcsd = [nn ;fcsd; nn]; % pad with nans for lost channels
            was = 1:size(mmuaB,2); % n contacts
          
            olddepths = was;
            ticks     = linspace(1,size(pd_fcsd,1),length(was));
            
            pr = 1;  %set lims for colorbar
            mn = ((min(min(fcsd))).*pr) ;
            mx = ((max(max(fcsd))).*pr) ;
            
            h = imagesc(time,1:size(pd_fcsd,1),pd_fcsd,[mn -mn]); cb = colorbar;
            set(h,'AlphaData',~isnan(pd_fcsd));

            map = flipud(colormap('jet'));
            colormap(map);
            ylim(cb,[mn -mn]);  %set lims for colorbar and add in white for nans
            set(cb,'YDir','reverse');
            set(gca,'YTick',ticks(2:2:end-1),'YTicklabel',olddepths(2:2:end-1),'FontSize',14); set(gca,'TickDir','out');
            xlabel('t (ms) from flash'); ylabel('channel number'); 
            set(get(cb,'ylabel'),'string','nA/mm^3','FontSize',14); 
            
            t = title(gca,sprintf('%s CSD %s',filesforday{fd}(1:end-4),whichElectrode),'interpreter','none');
            set(t,'FontWeight','bold');
            v = vline(0); set(v,'Color',[0.5 0.5 0.5],'LineWidth',1);
            h = text(post+80,csd(1,1)-15,sprintf('n = %d trs',length(onsetT))); set(h,'Color','k','FontSize',14);
            
            
        end
        
    end
end











