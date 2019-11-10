function plotCRFs(UNITDAT,id,anawin)

colors      = parula(10); 
unqc        = UNITDAT.contrast{id}.unqc;
nunqc       = UNITDAT.contrast{id}.nunqc;
spkdata     = UNITDAT.contrast{id}.data;
Ndata       = UNITDAT.contrast{id}.Ndata;
trspkdata   = UNITDAT.contrast{id}.trdata;
stats       = UNITDAT.contrast{id}.ttestpvals;


for c = 1:length(unqc)
    if ~isnan(trspkdata{c,1})
    DE_example(c)  = nanmean(nanmean(trspkdata{c,1}(anawin,:),1),2);
    BIN_example(c) = nanmean(nanmean(trspkdata{c,end}(anawin,:),1),2);
    end
end

[de_prd,de_xprd]   = runCRFFit(DE_example',unqc); 
[bin_prd,bin_xprd] = runCRFFit(BIN_example',unqc);


figure, 
plot(de_xprd,de_prd,'Color',colors(1,:)); 
hold on; 
plot(bin_xprd,bin_prd,'Color',colors(5,:)); 
hold on; 
plot(unqc.*100',DE_example,'o','Color',colors(1,:)); 
hold on; 
plot(unqc.*100',BIN_example,'o','Color',colors(5,:)); 
set(gca,'xscale','log','box','off','tickdir','out','linewidth',1.5); 
title(gca,num2str(id)); 