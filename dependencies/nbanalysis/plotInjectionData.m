
function plotInjectionData(DATA,pre,post,plot_types)


types = {'SDF','CSD','LFP','cTTL','cLFP','cSDF'};


preSDF  = DATA.preSDF; 
postSDF = DATA.postSDF; 
preCSD  = DATA.preCSD; 
postCSD = DATA.postCSD; 
preLFP  = DATA.preLFP; 
postLFP = DATA.postLFP; 
TTL     = DATA.TTL; 
cLFP    = DATA.cLFP; 
cSDF    = DATA.cSDF; 

pval_postSDF = DATA.pval_postSDF;
pval_postLFP = DATA.pval_postLFP; 
pvalSDF      = DATA.pvalSDF;
pvalLFP      = DATA.pvalLFP;
pval_clfp    = DATA.pval_clfp; 
pval_csdf    = DATA.pval_csdf; 
chans        = DATA.chans;



for t = 1:length(types)
    type = types{t};
    if ~any(ismember(plot_types,type)), continue, end
    switch type
        case 'SDF'
            
            [qts] = quantile(reshape(nanmean(preSDF,3),[numel(nanmean(preSDF,3)) 1]),[.1 .9]); 
            mn = qts(1); mx = qts(2); 
            figure
            set(gcf,'color','w','position',[1 1 1000 800]);
            
            subplot(1,2,1)
            imagesc(-pre:post,chans,nanmean(preSDF,3)');
            colorbar; caxis([mn mx]);
            vline(0);
            set(gca,'tickdir','out','ytick',chans);
            
            for ch = 1:length(chans)
                if pvalSDF(ch) < 0.05
                    text(-pre-30,chans(ch),'*');
                end
            end
            title(DATA.prefile,'interpreter','none'); 
            
            subplot(1,2,2)
            imagesc(-pre:post,chans,nanmean(postSDF,3)');
            cb = colorbar; vline(0); ylabel(cb,DATA.tags{ismember(plot_types,'SDF')}); 
            caxis([mn mx]); 
            set(gca,'tickdir','out','ytick',chans); colormap(gray)
            
            for ch = 1:length(chans)
                if pval_postSDF(ch) < 0.05
                    text(-pre-30,chans(ch),'*');
                end
            end
              title(DATA.pstfile,'interpreter','none'); 
            
        case 'CSD'


            [qts] = quantile(reshape(nanmean(preCSD,3),[numel(nanmean(preCSD,3)) 1]),[.1 .9]); 
          
            mn = -(max(abs(qts))); mx = -mn; 
            figure
            set(gcf,'color','w','position',[1 1 1000 800]);
            
            subplot(1,2,1)
            clear fcsd
            fcsd = [nan(length(-pre:post),10) filterCSD(nanmean(preCSD(:,2:end-1,:),3)')' nan(length(-pre:post),10)];
            h    = imagesc(-pre:post,chans,fcsd');
            colorbar; vline(0);
            set(gca,'tickdir','out','ytick',chans);
            set(h,'AlphaData',~isnan(fcsd'));
            for ch = 1:length(chans)
                if pvalLFP(ch) < 0.05
                    text(-pre-30,chans(ch),'*');
                end
            end
            caxis([mn mx]); 
             title(DATA.prefile,'interpreter','none'); 
             
            subplot(1,2,2)
            clear fcsd
            fcsd = [nan(length(-pre:post),10) filterCSD(nanmean(postCSD(:,2:end-1,:),3)')' nan(length(-pre:post),10)];
            h = imagesc(-pre:post,chans,fcsd');
            cb = colorbar; vline(0); ylabel(cb,DATA.tags{ismember(plot_types,'CSD')}); 
            colormap(flipud(colormap('jet')))
            set(gca,'tickdir','out','ytick',chans);
            set(h,'AlphaData',~isnan(fcsd'));
            for ch = 1:length(chans)
                if pval_postLFP(ch) < 0.05
                    text(-pre-30,chans(ch),'*');
                end
            end
            caxis([mn mx]); 
             title(DATA.pstfile,'interpreter','none'); 
        case 'LFP'
            [qts] = quantile(reshape(nanmean(preLFP,3),[numel(nanmean(preLFP,3)) 1]),[.1 .9]); 
            mn = -(max(abs(qts))); mx = -mn; 
            
            figure
            set(gcf,'color','w','position',[1 1 1000 800]);
            
            subplot(1,2,1)
            imagesc(-pre:post,chans,nanmean(preLFP,3)');
            colorbar; vline(0); caxis([mn mx]); 
            set(gca,'tickdir','out','ytick',chans);
            for ch = 1:length(chans)
                if pvalLFP(ch) < 0.05
                    text(-pre-30,chans(ch),'*');
                end
            end
            title(DATA.prefile,'interpreter','none'); 
             
            subplot(1,2,2)
            imagesc(-pre:post,chans,nanmean(postLFP,3)');
            cb = colorbar; ylabel(cb,DATA.tags{ismember(plot_types,'LFP')}); 
            vline(0); caxis([mn mx]); 
            set(gca,'tickdir','out','ytick',chans);
            for ch = 1:length(chans)
                if pval_postLFP(ch) < 0.05
                    text(-pre-30,chans(ch),'*');
                end
            end
            title(DATA.pstfile,'interpreter','none'); 
            
        case 'cTTL'
            
            pulses  = TTL > mean(TTL) + std(TTL);
            start_p = find(pulses,1,'first');
            end_p   = find(pulses,1,'last');
            
            minx    = (1:length(TTL))/1000/60;
            figure,set(gcf,'position',[90 600 1250 300],'color','w');
            plot(minx,TTL,'k'); vline(minx(start_p)); vline(minx(end_p));
            set(gca,'tickdir','out','box','off','linewidth',2); xlim([minx(1) minx(end)]);
            title(sprintf('injection duration %0.2f min\n%s', minx(end_p)-minx(start_p),DATA.pstfile),'interpreter','none'); 

        case 'cLFP'
            window = (1:size(cLFP,1))/1000/60;
            qts    = quantile(reshape(cLFP,[numel(cLFP) 1]),[.05 .95]); 
            mx     = qts(2); 
            mn     = qts(1);
            figure, set(gcf,'color','w','position',[1 1 1250 800]);
            for i = 1:length(chans)
                
                subplot(floor(length(chans)/4),4,i)
                plot(window,cLFP(:,i),'color',[.7 .7 .7])
                xlim([window(1) window(end)]);
                set(gca,'box','off','tickdir','out','linewidth',2);
                ylim([mn mx]); vline(window(start_p)); vline(window(end_p));
                if i == (4*((length(chans)/4)-1))+1
                    xlabel(sprintf('min\n%s',DATA.pstfile),'interpreter','none'); ylabel(DATA.tags{ismember(plot_types,'cLFP')});
                end
                pos = window(end) + (window(end)*.05);  
                if pval_clfp(i) < .05
                text(pos,mx.*.8,['*' num2str(chans(i))],'color','r'); 
                else
                text(pos,mx.*.8,num2str(chans(i))); 
                end
            end
            
        case 'cSDF'
            window = (1:size(cSDF,1))/1000/60;
            qts    = quantile(reshape(cSDF,[numel(cSDF) 1]),[.001 .999]); 
            mx     = qts(2); 
            mn     = qts(1);
            figure, set(gcf,'color','w','position',[1 1 1250 800]);
            for i = 1:length(chans)
                
                subplot(floor(length(chans)/4),4,i)
                plot(window,cSDF(:,i),'color',[.7 .7 .7])
                xlim([window(1) window(end)]);
                set(gca,'box','off','tickdir','out','linewidth',2);
                ylim([mn mx]); vline(window(start_p)); vline(window(end_p));
                if i == (4*((length(chans)/4)-1))+1
                    xlabel(sprintf('min\n%s',DATA.pstfile),'interpreter','none'); ylabel(DATA.tags{ismember(plot_types,'cSDF')});
                end
                   pos = window(end) + (window(end)*.05);  
                if pval_csdf(i) < .05
                text(pos,mx.*.8,['*' num2str(chans(i))],'color','r'); 
                else
                text(pos,mx.*.8,num2str(chans(i))); 
                end
            end
            
            
            
    end
    
    
end