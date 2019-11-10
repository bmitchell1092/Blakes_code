% analyze laminar V1 data for kanizsa paper:

clear;
set(0, 'DefaultAxesFontName','TimesRoman','DefaultAxesFontSize',20)
fdates = {'140724','140725','140807','140818','140731','140729'}; % '140820'

int    = 1;
saveon = 1;
kanpdfname = '/users/kaciedougherty/documents/analyzeddata/kanizsaMI_center_summary';

for s = 1:length(fdates)
    close all;
    clearvars -except fdates int saveon s kanpdfname ANADATA
    
    sessdate = fdates{s};
    [fnames] = getKanizsaFileList(sessdate);
    BRdatafile = fnames.evp;
    if ispc
        brdrname = strcat('\\129.59.230.179\CerebrusData\',BRdatafile(1:8));
        mldrname = strcat('\\129.59.230.116\MLData\',BRdatafile(1:8));
        
    elseif exist('/users/kaciedougherty/documents/code')
        %brdrname   = sprintf('/volumes/Toshiba External/%s',BRdatafile(1:8));
        brdrname   = sprintf('/users/kaciedougherty/documents/neurophysdata/%s',BRdatafile(1:8));
        mldrname   = brdrname;
        savedrname = '/users/kaciedougherty/documents/analyzeddata/';
    else
        if str2num(BRdatafile(1:6))>160611
            brdrname = sprintf('/Volumes/Drobo2/DATA/NEUROPHYS/rig021/%s',BRdatafile(1:8));
        else
            brdrname = sprintf('/Volumes/Drobo/DATA/NEUROPHYS/rig021/%s',BRdatafile(1:8));
        end
        mldrname = brdrname;
        spath    = '/volumes/drobo/users/kacie/analysis/grcposter';
    end
    %%
    extension = 'ns2';
    el = fnames.probe;
    
%     plotEVP([fnames.evp],brdrname,extension,el);
%     title(gca,BRdatafile(1:6));
%     subplot(1,4,1), title(gca,'single trial LFP');
%     subplot(1,4,2), title(gca,'avg LFP');
%     subplot(1,4,3), title(gca,'avg CSD');
% %     if ~exist(strcat(kanpdfname,'.pdf'))
% %         export_fig(kanpdfname,'-pdf','-nocrop');
% %     else
% %         export_fig(kanpdfname,'-pdf','-nocrop','-append');
% %     end
%%
% center file one:
[trLFP,trTHETA,trGAMMA,onsets,conds,PRE,POST,Fs] = loadKanizsaFile([fnames.kanizsa{1}],brdrname,extension,el);


[~,~,sink,chans] = getSink(BRdatafile(1:6));
trLFP   = trLFP(:,chans,:);
trTHETA = trTHETA(:,chans,:);
trGAMMA = trGAMMA(:,chans,:);

    %% BOTH FILES TOGETHER
    for c = 1:2
        clear thesetrs hMI
        thesetrs = find(conds == c);
        tic
        [hMI] = runMILFP(trTHETA(:,:,thesetrs),trGAMMA(:,:,thesetrs),Fs);
        toc
        MI(:,:,c) = squeeze(hMI);
    end
    
    if c == 2
        ANADATA(s).center = MI;
    end
    
    if int == 1
        kMI   = fliplr(rot90(interp2(MI(:,:,1),10,'cubic'),3));
        conMI = fliplr(rot90(interp2(MI(:,:,2),10,'cubic'),3));
    else
        kMI   = fliplr(rot90(MI(:,:,1),3));
        conMI = fliplr(rot90(MI(:,:,2),3));
    end
    
    [~,~,sink,chans] = getSink(BRdatafile(1:6));
    totchan    = length(chans);
    depths     = [fliplr([chans(1):sink-1]*100) 0 : -100: ((chans(end)-sink)*-100)];
    
    if c == 2
        ANADATA(s).depths = depths;
    end
    
    figure, set(gcf,'Color','w','Position',[1 1 1500 600]);
    subplot(1,2,1)
    imagesc(depths,depths,kMI),
    colormap('jet'),cb1 = colorbar;
    title(gca,strcat('kanizsa 1',' n=',num2str(length(find(conds == 1))),' sess ',BRdatafile(1:6)));
    set(gca,'YDir','normal','XDir','reverse');
    xlabel('depth for theta (microns)');
    v = vline(0); h = hline(0); set(v,'LineWidth',2,'Color','w'); set(h,'LineWidth',2,'Color','w');
    ylabel('depth for gamma (microns)');axis square
    set(gca,'FontSize',20); 
 
    subplot(1,2,2)
    imagesc(depths,depths,conMI),cb2 = colorbar;
  set(get(cb2,'ylabel'),'string','MI');
    title(gca,strcat('control 2',' n=',num2str(length(find(conds == 2))),' sess ',BRdatafile(1:6)));
    set(gca,'YDir','normal','XDir','reverse');
    first = get(cb1,'Limits'); caxis([first]);
    xlabel('depth for theta (microns)');
    v = vline(0); h = hline(0); set(v,'LineWidth',2,'Color','w'); set(h,'LineWidth',2,'Color','w');
    ylabel('depth for gamma (microns)');axis square
       set(gca,'FontSize',20); 
 
    if ~exist(strcat(kanpdfname,'.pdf'))
        export_fig(kanpdfname,'-pdf','-nocrop');
    else
        export_fig(kanpdfname,'-pdf','-nocrop','-append');
    end
    
    avgkan = nanmean(nanmean(MI(:,:,1),1),2);
    avgcon = nanmean(nanmean(MI(:,:,2),1),2);
    
    figure, set(gcf,'Color','w','Position',[1 1 620 300]);
    h1 = bar([nanmean(avgkan)],'FaceColor',[.8 .2 .2]);
    hold on
    h2 = bar(2,[nanmean(avgcon)],'FaceColor',[.2 .2 .8]);
    set(gca,'Box','off','TickDir','out','FontSize',20,'XTick',[]);
    legend([h1 h2],'kanizsa','control'); ylabel('MI');
    title(gca,strcat('sess ',BRdatafile(1:6)));  set(gca,'FontSize',20); 
    export_fig(kanpdfname,'-pdf','-nocrop','-append');

    
end
dfdf
%%

clear tlab emparray  empdepths datadepths plotdepths mat


cond       = size(ANADATA(1).center,3);
emparray   = nan(totchan*2,totchan*2,cond,length(ANADATA));

empdepths    = [2400:-100:-2400];
nsess = 0;
for s = 1:length(ANADATA)
    
    
            [~,~,sink,chans] = getSink(fdates{s});
            totchan    = length(chans);
            datadepths     = [fliplr([chans(1):sink-1]*100) 0 : -100: ((chans(end)-sink)*-100)];
   
    
            [~,dataid]    = ismember(empdepths,datadepths);
    
    
            emparray(find(dataid),find(dataid),:,s)   = ANADATA(s).center;
    
    avgkan(s) = nanmean(nanmean(ANADATA(s).center(:,:,1),1),2);
    avgcon(s) = nanmean(nanmean(ANADATA(s).center(:,:,2),1),2);
    nsess = nsess + 1;
end
tlabel = 'center';

figure, set(gcf,'Color','w','Position',[1 1 620 300]);
h1 = bar([nanmean(avgkan)],'FaceColor',[.8 .2 .2]);
hold on;
errorbar(1,[nanmean(avgkan)],[std(avgkan,0,2)],'Color',[.8 .2 .2]);
hold on;
h2 = bar(2,[nanmean(avgcon)],'FaceColor',[.2 .2 .8]);
hold on;
errorbar(2,[nanmean(avgcon)],[std(avgcon,0,2)],'Color',[.2 .2 .8]);
set(gca,'Box','off','TickDir','out','FontSize',20,'XTick',[]); 
legend([h1 h2],'kanizsa','control'); ylabel('MI'); 


plotdepths = [1200:-100:-600];
[~,pid] = ismember(empdepths,plotdepths);

figure, set(gcf,'Color','w');
subplot(1,2,1)
mat = emparray(find(pid),find(pid),1,:);

imagesc(plotdepths,plotdepths,fliplr(rot90(nanmean(mat,4),3)));
colormap('jet'),cb1 = colorbar;
set(gca,'YDir','normal','XDir','reverse');
xlabel('depth for theta (microns)');ylabel('depth for gamma (microns)');
v = vline(0); h = hline(0); set(v,'LineWidth',2,'Color','w'); set(h,'LineWidth',2,'Color','w');
axis square
title(gca,strcat(tlabel,' KANIZSA',' n =',num2str(nsess),'sessions'));

% figure, surf(plotdepths,plotdepths,fliplr(rot90(nanmean(mat,4),3)));

subplot(1,2,2)
clear mat
mat = emparray(find(pid),find(pid),2,:);
imagesc(plotdepths,plotdepths,fliplr(rot90(nanmean(mat,4),3)));
colormap('jet'),cb2 = colorbar;
set(gca,'YDir','normal','XDir','reverse');
xlabel('depth for theta (microns)');ylabel('depth for gamma (microns)');
v = vline(0); h = hline(0); set(v,'LineWidth',2,'Color','w'); set(h,'LineWidth',2,'Color','w');
axis square
title(gca,strcat(tlabel,' CONTROL'));
subplot(1,2,1)
caxis([0 max([cb1.Limits(2) cb2.Limits(2)])]);
subplot(1,2,2)
caxis([0 max([cb1.Limits(2) cb2.Limits(2)])]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure, set(gcf,'Color','w');
subplot(1,2,1)
mat = interp2(nanmean(emparray(find(pid),find(pid),1,:),4),10,'cubic');
imagesc(plotdepths,plotdepths,fliplr(rot90(nanmean(mat,4),3)));
colormap('jet'),cb1 = colorbar;
set(gca,'YDir','normal','XDir','reverse');
xlabel('depth for theta (microns)');ylabel('depth for gamma (microns)');
v = vline(0); h = hline(0); set(v,'LineWidth',2,'Color','w'); set(h,'LineWidth',2,'Color','w');
axis square
title(gca,strcat(tlabel,' KANIZSA',' n =',num2str(num2str(nsess)),'sessions'));

subplot(1,2,2)
clear mat
mat = interp2(nanmean(emparray(find(pid),find(pid),2,:),4),10,'cubic');
imagesc(plotdepths,plotdepths,fliplr(rot90(nanmean(mat,4),3)));
colormap('jet'),cb2 = colorbar;
set(gca,'YDir','normal','XDir','reverse');
xlabel('depth for theta (microns)');ylabel('depth for gamma (microns)');
v = vline(0); h = hline(0); set(v,'LineWidth',2,'Color','w'); set(h,'LineWidth',2,'Color','w');
axis square
title(gca,strcat(tlabel,' CONTROL'));
subplot(1,2,1)
caxis([0 max([cb1.Limits(2) cb2.Limits(2)])]);
subplot(1,2,2)
caxis([0 max([cb1.Limits(2) cb2.Limits(2)])]);

saveas(figure(1),'AVG_center','jpg');
saveas(figure(2),'AVG_center_interp','jpg');

close all;



