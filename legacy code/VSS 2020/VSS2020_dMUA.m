%% PANEL 
% choices
tw = 1; 
cont = 4;

% data
clear mon bin
for L = 1:3  %store mon and bin units into struct
mon(L).units = LAY.MON.DE_PS(L).RESP;
bin(L).units = LAY.BIN.PS(L).RESP;
end

for L = 1:3 %retrieve avg and error
    mon(L).avg = nanmean(mon(L).units,3);
    mon(L).err = nanstd(mon(L).units,[],3)./sqrt(size(mon(L).units,3));
    bin(L).avg = nanmean(bin(L).units,3);
    bin(L).err = nanstd(bin(L).units,[],3)./sqrt(size(bin(L).units,3));
end

% Plot
figure('position',[573,347,892,510]);
for L = 1:3
subplot(3,1,L);
distributionPlot([squeeze(mon(L).units(cont,tw,:)) squeeze(bin(L).units(cont,tw,:))],...
    'showMM',4,'color',{[0.2 0.2 0.2],[0.1, 0.2, 0.3]},...
    'xNames',{'MON','BIN'},'globalNorm',2,'addSpread',1)
hold on
set(gca,'box','off','linewidth',2,'FontSize',16,...
    'ylim',[-100 max(bin(1).units(cont,tw,:))*1.2]);

    if L == 3
        title(sprintf('Upper',layerLengths(1)),'FontSize',18);
        ylabel('impulses per sec','Fontsize',18);
    elseif L == 2
            title(sprintf('Middle',layerLengths(2)),'FontSize',18);
            yticklabels([]);
        else 
            title(sprintf('Deep',layerLengths(3)),'FontSize',18);  
            yticklabels([]);
    end
end

% mCt = sum(Trls.mon(:,1,:),'all'); bCt = sum(Trls.bin(:,1,:),'all');

if flag_figsave == 1
    cd(strcat(figDir,'layers\'));
    saveas(gcf, strcat('MONvsBIN_violin', '.svg'));
    disp("Figure saved");
else
    disp("Figure was not saved");
end



%% Panel 1: figure 1
clear LSM AVE QSM SUP mdl Model1 Model2
cont = 3;
tw = 3;
pen = 12;
Model1 = 'QSM';
Model2 = 'SUP';

clear *de* *nde* *bin*
switch cont
    case 2
        cLevel = '22';
    case 3
        cLevel = '45';
    case 4
        cLevel = '90';
end

depth = 12:-1:-4;

% data
de.units = PEN(12).MON.DE_PS.RESP;
nde.units = PEN(pen).MON.NDE_PS.RESP;
bin.units = PEN(pen).BIN.PS.RESP;

% models
[~,~,QSM,SUP] = modelAnalysis(de.units,nde.units);

% plot 1
figure('position',[182,227,1019,563]);
subplot(1,4,1)
plot(squeeze(de.units(cont,tw,:)), depth,'linewidth',1.5,'Color','k'); hold on
plot(squeeze(nde.units(cont,tw,:)), depth,'linewidth',1.5,'Color',[0, 0, 0]+0.8);
grid off

% gca
set(gca,'box','off','linewidth',1.5,'fontsize',14,...
    'xlim',[-5 350],'ylim',[-4 12]);

% labels
yticklabels({'-0.4','-0.2','0','0.2','0.4','0.6','0.8','1.0','1.2'})
xlabel('impulses per sec','FontSize',16);
ylabel('Depth (mm) relative to layer 4/5 boundary','FontSize',16);
legend(sprintf('%s DE',cLevel),sprintf('%s NDE',cLevel),'Location','southeast','orientation','vertical'); legend boxoff
title(sprintf('Monocular\n responses'),'FontSize',16);

hold off

% plot 2
subplot(1,4,2)
plot(squeeze(bin.units(cont,tw,:)), depth,'linewidth',1.5,'Color','b');
hold on 
plot(squeeze(QSM(cont,tw,:)), depth,'linestyle',':','linewidth',2,'Color',[0, 0, 0]+0.3);
grid off
xlim([-5 350]);

% gca
set(gca,'box','off','linewidth',1.5,'fontsize',14,...
    'xlim',[-5 350],'ylim',[-4 12]);

% labels
yticklabels({'-0.4','-0.2','0','0.2','0.4','0.6','0.8','1.0','1.2'})
yticklabels([]); xticklabels([]);
%xlabel('impulses per sec','FontSize',16);
%ylabel('Depth (mm) relative to layer 4/5 boundary','FontSize',16);
legend('BIN','QSM','Location','southeast','orientation','vertical'); legend box off
title(sprintf('Quadratic \nSummation (QS)'),'FontSize',16);
hold off

% plot 3
subplot(1,4,3)
plot(squeeze(bin.units(cont,tw,:)), depth,'linewidth',1.5,'Color','b');
hold on 
plot(squeeze(SUP(cont,tw,:)), depth,'linestyle','-.','linewidth',2,'Color',[0, 0, 0]+0.3);
grid off
xlim([-5 350]);

% gca
set(gca,'box','off','linewidth',1.5,'fontsize',14,...
    'xlim',[-5 350],'ylim',[-4 12]);

% labels
yticklabels({'-0.4','-0.2','0','0.2','0.4','0.6','0.8','1.0','1.2'})
yticklabels([]); xticklabels([]);
%xlabel('impulses per sec','FontSize',16);
%ylabel('Depth (mm) relative to layer 4/5 boundary','FontSize',16);
legend('BIN','SUP','Location','southeast','orientation','vertical'); legend boxoff
title(sprintf('QS with\nsuppressive term'),'FontSize',16);
hold off

% plot 4
subplot(1,4,4)
plot(squeeze(QSM(cont,tw,:)) - squeeze(bin.units(cont,tw,:)), depth,'linestyle',':','linewidth',2,'Color',[0,0,0]+0.3);
hold on 
plot(squeeze(SUP(cont,tw,:)) - squeeze(bin.units(cont,tw,:)), depth,'linestyle','--','linewidth',2,'Color',[0,0,0]+0.3);
grid off
xlim([-5 350]);

% gca
set(gca,'box','off','linewidth',1.5,'fontsize',14,...
    'xlim',[-140 140],'ylim',[-4 12]);
yticklabels([]);
vl = vline(0, 'k');
set(vl,'linewidth',1);

% labels
yticklabels({'-0.4','-0.2','0','0.2','0.4','0.6','0.8','1.0','1.2'})
yticklabels([]); 
%xlabel('residual','FontSize',16);
%ylabel('Depth (mm) relative to layer 4/5 boundary','FontSize',16);
%legend('QSM error','SUP error','Location','southeast','orientation','vertical'); legend boxoff
title(sprintf('Model\n error'),'FontSize',16);
hold off

if flag_figsave == 1
    cd(strcat(figDir,'probe\'));
    saveas(gcf, strcat('BINvsQSM-SUP', '.svg'));
    disp("Figure saved");
else
    disp("Figure was not saved");
end


%% Panel 2: figure 1
clear LSM AVE QSM SUP mdl Model1 Model2
cont = 3;
tw = 3;
pen = 12;
Model1 = 'LSM';
Model2 = 'AVE';

clear *de* *nde* *bin*
switch cont
    case 2
        cLevel = '22';
    case 3
        cLevel = '45';
    case 4
        cLevel = '90';
end

depth = 12:-1:-4;

% data
de.units = PEN(12).MON.DE_PS.RESP;
nde.units = PEN(pen).MON.NDE_PS.RESP;
bin.units = PEN(pen).BIN.PS.RESP;

% models
[LSM,AVE,~,~] = modelAnalysis(de.units,nde.units);

% plot 1
figure('position',[182,227,1019,563]);
subplot(1,4,1)
plot(squeeze(de.units(cont,tw,:)), depth,'linewidth',1.5,'Color','k'); hold on
plot(squeeze(nde.units(cont,tw,:)), depth,'linewidth',1.5,'Color',[0, 0, 0]+0.8);
grid off

% gca
set(gca,'box','off','linewidth',1.5,'fontsize',14,...
    'xlim',[-5 350],'ylim',[-4 12]);

% labels
yticklabels({'-0.4','-0.2','0','0.2','0.4','0.6','0.8','1.0','1.2'})
xlabel('impulses per sec','FontSize',16);
ylabel('Depth (mm) relative to layer 4/5 boundary','FontSize',16);
legend(sprintf('%s DE',cLevel),sprintf('%s NDE',cLevel),'Location','southeast','orientation','vertical'); legend boxoff
title(sprintf('Monocular\n responses'),'FontSize',16);

hold off

% plot 2
subplot(1,4,2)
plot(squeeze(bin.units(cont,tw,:)), depth,'linewidth',1.5,'Color','b');
hold on 
plot(squeeze(LSM(cont,tw,:)), depth,'linestyle',':','linewidth',2,'Color',[0, 0, 0]+0.3);
grid off
xlim([-5 350]);

% gca
set(gca,'box','off','linewidth',1.5,'fontsize',14,...
    'xlim',[-5 350],'ylim',[-4 12]);

% labels
yticklabels({'-0.4','-0.2','0','0.2','0.4','0.6','0.8','1.0','1.2'})
yticklabels([]); xticklabels([]);
%xlabel('impulses per sec','FontSize',16);
%ylabel('Depth (mm) relative to layer 4/5 boundary','FontSize',16);
legend('BIN','LS','Location','southeast','orientation','vertical'); legend box off
title(sprintf('Linear \nSummation (LS)'),'FontSize',16);
hold off

% plot 3
subplot(1,4,3)
plot(squeeze(bin.units(cont,tw,:)), depth,'linewidth',1.5,'Color','b');
hold on 
plot(squeeze(AVE(cont,tw,:)), depth,'linestyle','--','linewidth',2,'Color',[0, 0, 0]+0.3);
grid off
xlim([-5 350]);

% gca
set(gca,'box','off','linewidth',1.5,'fontsize',14,...
    'xlim',[-5 350],'ylim',[-4 12]);

% labels
yticklabels({'-0.4','-0.2','0','0.2','0.4','0.6','0.8','1.0','1.2'})
yticklabels([]); xticklabels([]);
%xlabel('impulses per sec','FontSize',16);
%ylabel('Depth (mm) relative to layer 4/5 boundary','FontSize',16);
legend('BIN','DIV','Location','southeast','orientation','vertical'); legend boxoff
title(sprintf('Binocular\n Division'),'FontSize',16);
hold off

% plot 4
subplot(1,4,4)
plot(squeeze(LSM(cont,tw,:)) - squeeze(bin.units(cont,tw,:)), depth,'linestyle',':','linewidth',2,'Color',[0,0,0]+0.3);
hold on 
plot(squeeze(AVE(cont,tw,:)) - squeeze(bin.units(cont,tw,:)), depth,'linestyle','--','linewidth',2,'Color',[0,0,0]+0.3);
grid off
xlim([-5 350]);

% gca
set(gca,'box','off','linewidth',1.5,'fontsize',14,...
    'xlim',[-140 140],'ylim',[-4 12]);
yticklabels([]);
vl = vline(0, 'k'); 
set(vl,'linewidth',1);

% labels
yticklabels({'-0.4','-0.2','0','0.2','0.4','0.6','0.8','1.0','1.2'})
yticklabels([]); 
%xlabel('residual','FontSize',16);
%ylabel('Depth (mm) relative to layer 4/5 boundary','FontSize',16);
%legend('QSM error','SUP error','Location','southeast','orientation','vertical'); legend boxoff
title(sprintf('Model\n error'),'FontSize',16);
hold off

if flag_figsave == 1
    cd(strcat(figDir,'probe\'));
    saveas(gcf, strcat('BINvsLSM-AVE', '.svg'));
    disp("Figure saved");
else
    disp("Figure was not saved");
end
%% Model vs SDF by contrast
models = {'LSM','AVE','QSM','SUP'};
clear *LSM* *QSM* *AVE* *SUP* mdl
x = sdfWin;

% choices
cont = 4;
model = 'AVE';


% data
clear de nde bin
for L = 1:3
    de(L).units = LAY.MON.DE_PS(L).SDF;
    nde(L).units = LAY.MON.NDE_PS(L).SDF;
    bin(L).units = LAY.BIN.PS(L).SDF;
end

for L = 1:3
    de(L).avg = squeeze(nanmean(de(L).units,3));
    de(L).err = squeeze(nanstd(de(L).units,[],3)./sqrt(size(de(L).units,3)));
    nde(L).avg = squeeze(nanmean(nde(L).units,3));
    nde(L).err = squeeze(nanstd(nde(L).units,[],3)./sqrt(size(nde(L).units,3)));
    bin(L).avg = nanmean(bin(L).units,3);
    bin(L).err = nanstd(bin(L).units,[],3)./sqrt(size(bin(L).units,3));
end

for L = 1:3
    [LSM(L).units,AVE(L).units,QSM(L).units,SUP(L).units] = modelAnalysis(de(L).units,nde(L).units);
end

for L = 1:3
    LSM(L).avg = nanmean(LSM(L).units,3);
    LSM(L).err = nanstd(LSM(L).units,[],3)./sqrt(size(LSM(L).units,3));
    AVE(L).avg = nanmean(AVE(L).units,3);
    AVE(L).err = nanstd(AVE(L).units,[],3)./sqrt(size(AVE(L).units,3));
    QSM(L).avg = nanmean(QSM(L).units,3);
    QSM(L).err = nanstd(QSM(L).units,[],3)./sqrt(size(QSM(L).units,3));
    SUP(L).avg = nanmean(SUP(L).units,3);
    SUP(L).err = nanstd(SUP(L).units,[],3)./sqrt(size(SUP(L).units,3));
end


switch cont
    case 4
        c = 90;
    case 3
        c = 45;
    case 2 
        c = 22;
end

switch model
    case 'LSM'
        mdl = LSM;
    case 'AVE'
        mdl = AVE;
    case 'QSM'
        mdl = QSM;
    case 'SUP'
        mdl = SUP;
end
% Plot  
figure('position',[784.4285714285713,76.42857142857143,368.5714285714287,847.4285714285713]);
for L = 1:3
    subplot(3,1,L)
plot(x,bin(L).avg(cont,:),'b','linewidth',2); hold on
ci1 = ciplot(bin(L).avg(cont,:)+bin(L).err(cont,:),...
    bin(L).avg(cont,:)-bin(L).err(cont,:),x,'b',0.1); set(ci1,'linestyle','none','handleVisibility','off');
%plot(x,mdl(L).avg(cont,:),'k','linewidth',1); 
ci2 = ciplot(mdl(L).avg(cont,:)+mdl(L).err(cont,:),...
    mdl(L).avg(cont,:)-mdl(L).err(cont,:),x,'k',0.1); set(ci2,'linestyle','none','handleVisibility','off');

[h,p] = ttest_time(bin(L).units,mdl(L).units); tmh = find(h(cont,:));
%scatter(sdfWin(tmh), ones(1,numel(tmh)) * max(bin(L).avg(cont,:)) * 1.2,4,'*k') % significance bar

ylimit = max(bin(L).avg(cont,:));
set(gca,'Box','off','TickDir','out','linewidth',2,'ylim',[-10 250],'xlim',[0 .300],'FontSize',16)
% [-10 ylimit*1.5]
    if L == 1
        %title(sprintf('Upper (n = %d)',layerLengths(1)),'FontSize',16);
        bCt = sum(Trls.bin(cont,1,:),'all'); 
        %legend(sprintf('BIN %d|%d (%d trials)',c,c,bCt),sprintf('Model %s',model),'location','northwest'); legend boxoff
        legend(sprintf('BIN %d|%d',c,c),sprintf('Model %s',model),'location','southeast'); legend boxoff
        ylabel([]);
        ylabel([]);
        xticklabels([]); xlabel([]);
    elseif L == 2
            %title(sprintf('Middle (n = %d)',layerLengths(2)),'FontSize',16);
            ylabel([]);
            xticklabels([]); xlabel([]);
        else 
            %title(sprintf('Deep (n = %d)',layerLengths(3)),'FontSize',16); 
            if cont == 2
            xlabel('time (s) from stim onset','Fontsize',18); 
            ylabel('impulses per sec','Fontsize',18);
            end
    end
    
    if c == 90
        xticklabels([]); xlabel([]);
        yticklabels([]); ylabel([]);
    end

end

%sgtitle(sprintf('dMUA | N = %d | Penetrations: %d \nContrast: %d|%d | Model: %s',uct,N,c,c,model),'Fontsize',16); 
if flag_figsave == 1
    cd(strcat(figDir,'layers\'));
    saveas(gcf, strcat(model,'_',num2str(c),'_sdf', '.svg'));
    disp("Figure saved");
else
    disp("Figure was not saved");
end
