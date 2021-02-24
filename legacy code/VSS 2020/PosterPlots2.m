%% Panel 1: figure 1
figure('position',[151,136,979,496]);
clear i
for i = 1:3
subplot(1,3,i)
plot(squeeze(sAVG.DE.aMUA.pc.coll.aligned(i+1,2,:)),corticaldepth,'linewidth',1.5,'Color',[0, 0.4470, 0.7410]);
hold on 
plot(squeeze(sAVG.BIN.aMUA.pc.coll.aligned(i+1,2,:)),corticaldepth,'linewidth',1.5','color','r');
grid on
xlim([-5 80]);
%hline(0,':','BOL4')
ylim([-7 17])
set(gca,'Box','off','linewidth',1)
hold off
    if i == 1
        title('Low Contrast','FontSize',12);
        legend('One eye','Both eyes','Location','northeast','orientation','vertical');
            xlabel('Multi-unit activity (MUA)','FontSize',12);
            yticklabels({'-0.5','0','0.5','1','1.5','2','2.5'})
            ylabel('Depth (mm) relative to layer 4/5 boundary','FontSize',12);
    elseif i == 2
            title('Medium Contrast','FontSize',12);
            yticklabels([]);
            xticklabels([]);
        else 
            title('High Contrast','FontSize',12);
            xticklabels([]);
            yticklabels([]);
    end
end

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\VSS\')
saveas(gcf, strcat('panel_1_fig1', '.svg'));

%% Panel 1: figure 2
figure('position',[360,41.666666666666664,308.3333333333333,599.3333333333333]);
refwin = -50:600;
clear c L
for L = 1:3
    subplot(3,1,L)
    plot(refwin,sAVG.DE.aMUA.pc.layers.data(2:4,:,L),'color','b','linewidth',1)
    hold on
    plot(refwin,sAVG.BIN.aMUA.pc.layers.data(2:4,:,L),'color','r','linewidth',1)
    xlim([-50 600])
    %ylimit = max(abs(get(gcf,'ylim')));
    ylimit = 90;
    set(gca,'ylim',[-10 ylimit],'Box','off','TickDir','out','linewidth',1)
    
    if L == 1
        title('Supragranular','FontSize',12);
        %legend({'BIN',NaN,'DE'},'location','northwest');
        yticklabels([]);
        xticklabels([]);
    elseif L == 2
            title('Granular','FontSize',12);
            yticklabels([]);
            xticklabels([]);
        else 
            title('Infragranular','FontSize',12);
            xlabel('time (ms)','FontSize',12);
            %xticklabels(labels);
            %set(gca,'xticklabels', 'FontSize',8);
%             lgd = legend('MON','BIN','location','northwest');
%             lgd.FontSize = 8;
            ylabel(sprintf('MUA response'),'FontSize',12);

    end
end

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\VSS\')
saveas(gcf, strcat('panel_1_fig2_v2', '.svg'));

%% Panel 1: figure 3
%rcontrast = round(tmp.STIM.levels,2,'significant');

figure('position',[360,41.666666666666664,259,599.3333333333333]);

labels = {0, 22, 45, 90}; format bank;

clear L c
for L = 1:3
subplot(3,1,L)
bar(sAVG.BIN.aMUA.pc.coll_layers.data(:,2,L),.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar(sAVG.BIN.aMUA.pc.coll_layers.data(:,2,L),sAVG.BIN.aMUA.pc.coll_layers.error(:,2,L),'o','marker','none','color','k');
bar(sAVG.DE.aMUA.pc.coll_layers.data(:,2,L),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
errorbar(sAVG.DE.aMUA.pc.coll_layers.data(:,2,L),sAVG.DE.aMUA.pc.coll_layers.error(:,2,L),'o','marker','none','color','k');
set(gca,'box','off','linewidth',1);
xticklabels([]);
ylim([0 60]);
%xlabel('% contrast')


    if L == 1
        title('Supragranular','FontSize',12);
        %legend({'BIN',NaN,'DE'},'location','northwest');
        yticklabels([]);
        xticklabels([]);
    elseif L == 2
            title('Granular','FontSize',12);
            yticklabels([]);
            xticklabels([]);
        else 
            title('Infragranular','FontSize',12);
            xlabel('stimulus contrast','FontSize',12);
            xticklabels(labels)
%              lgd = legend('BIN','MON','location','northeast');
%              lgd.FontSize = 8;
            ylabel(sprintf('MUA response'),'FontSize',12);
    end
hold off
end

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\VSS\')
saveas(gcf, strcat('panel_1_fig3_v2', '.svg'));

%% Panel 2: figure 1
cIndex = 4;

switch cIndex
    case 2
        cLevel = '.22';
    case 3
        cLevel = '.45';
    case 4
        cLevel = '.90';
end

figure('position',[151,58.333,834.66,574.66]);
clear i
subplot(1,3,1)
plot(squeeze(sAVG.DE.aMUA.pc.coll.aligned(cIndex,2,:)),corticaldepth,'linewidth',1.5,'Color',[0, 0.4470, 0.7410]);
hold on 
plot(squeeze(sAVG.NDE.aMUA.pc.coll.aligned(cIndex,2,:)),corticaldepth,'linewidth',.8','Color', [0, 0.4470, 0.7410]);
grid on
set(gca,'box','off','linewidth',1);
xlim([-5 80]);
%hline(0,':','BOL4')
ylim([-7 17])
yticklabels({'-0.5','0','0.5','1','1.5','2','2.5'})
xlabel('MUA (percent change)','FontSize',12);
legend(sprintf('%s DE',cLevel),sprintf('%s NDE',cLevel),'Location','northeast','orientation','vertical');
title('Monocular responses','FontSize',12);
ylabel('Depth (mm) relative to layer 4/5 boundary','FontSize',12);
hold off

subplot(1,3,2)
plot(squeeze(sAVG.BIN.aMUA.pc.coll.aligned(cIndex,2,:)),corticaldepth,'Color','r','linewidth',1.5);
hold on 
%plot(squeeze(sAVG.BIN.aMUA.pc_LSM.coll.aligned(cIndex,2,:)),corticaldepth,'linewidth',1.5,'linestyle','--','color',[.35 .4 .3]);
plot(squeeze(sAVG.BIN.aMUA.pc_QSM.coll.aligned(cIndex,2,:)),corticaldepth,'linewidth',1.5,'linestyle',':','color','k');
grid on
xlim([-5 80]);
%xticklabels([]);
xlabel('MUA (percent change)','FontSize',12);
yticklabels([]);
ylim([-7 17])
set(gca,'box','off','linewidth',1);
title(sprintf('Binocular response \nand Model prediction'),'FontSize',12);
legend('BIN','QSM','Location','northeast','orientation','vertical');
hold off

subplot(1,3,3)
%plot(squeeze(sAVG.BIN.aMUA.pc_LSM.coll.aligned(cIndex,2,:)-(sAVG.BIN.aMUA.pc.coll.aligned(cIndex,2,:))),corticaldepth,'linewidth',1.5,'linestyle','--','color',[.35 .4 .3]);
plot(squeeze(sAVG.BIN.aMUA.pc_QSM.coll.aligned(cIndex,2,:)-(sAVG.BIN.aMUA.pc.coll.aligned(cIndex,2,:))),corticaldepth,'linewidth',1.5,'linestyle',':','color','k');
hold on
grid on
xlim([-25 25]);
ylim([-7 17])
%xticklabels([]);
xlabel('MUA difference','FontSize',12)
yticklabels([]);
vl = vline(0, 'k');
set(vl,'linewidth',1);
set(gca,'box','off','linewidth',1);
title(sprintf('Model difference \nfrom Binocular response'),'FontSize',12);
legend('QSM','Location','northeast','orientation','vertical');
hold off

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\VSS\')
saveas(gcf, strcat('panel_2_fig1_v2', '.svg'));

%% Panel 2: figure 2

figure('position',[259.6666666666666,41.666666666666664,449.3333333333334,599.3333333333333]);
cIndex = 4;

switch cIndex
    case 2
        cLevel = '.22';
    case 3
        cLevel = '.45';
    case 4
        cLevel = '.90';
end

clear i L
for L = 1:3
subplot(3,1,L)
%plot(smooth(sAVG.BIN.aMUA.pc_LSM.layers.data(cIndex,:,L)-sAVG.BIN.aMUA.pc.layers.data(cIndex,:,L),0.1,'rloess'),'--','linewidth',1.5,'color','k');
plot(smooth(sAVG.BIN.aMUA.pc_QSM.layers.data(cIndex,:,L)-sAVG.BIN.aMUA.pc.layers.data(cIndex,:,L),0.1,'rloess'),':','linewidth',1.5,'color','k');
ylim([-25 25])
xlim([0 600])
hline(0,'k')
set(gca,'linewidth',1,'box','off')
%vline(tmp.STIM.off,'k','offset')
   if L == 1
        title('Supragranular','FontSize',12);
        xticklabels([]);
        yticklabels([]);
    elseif L == 2
            title('Granular','FontSize',12);
            xticklabels([]);
            yticklabels([]);
        else 
            title('Infragranular','FontSize',12);
            lgd = legend('QSM residual','location','southeast');
            lgd.FontSize = 8;
            xlabel('time (ms)','FontSize',12);
            ylabel(sprintf('model error'),'FontSize',12);
    end
end

%sgtitle(sprintf('Model error across stimulus duration\nDioptic: %s DE | %s NDE',cLevel,cLevel));


cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\VSS\')
saveas(gcf, strcat('panel_2_fig2_v2', '.svg'));


%% Panel 2: figure 3
labels = {0, 22, 45, 90,[],0,22,45,90}; format bank;
figure('Position', [140,189,784.3333333333333,436]);

test = sAVG.BIN.aMUA.pc.coll_layers.data(:,:,:)-sAVG.BIN.aMUA.pc.coll_layers.data(1,:,:);
test2 = sAVG.BIN.aMUA.pc_QSM.coll_layers.data(:,:,:)-sAVG.BIN.aMUA.pc_QSM.coll_layers.data(1,:,:);
clear L c
for L = 1:3
subplot(2,3,L)
%bar([sAVG.BIN.aMUA.pc_LSM.coll_layers.data(:,2,L);NaN;sAVG.BIN.aMUA.pc_LSM.coll_layers.data(:,3,L)],0.8,'FaceColor',[1, 1, 1],'linestyle','--','EdgeColor','k','LineWidth',0.8);
%hold on
%bar([sAVG.BIN.aMUA.pc.coll_layers.data(:,2,L);NaN;sAVG.BIN.aMUA.pc.coll_layers.data(:,3,L)],0.6,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',1);
bar(test(:,2,L),.7,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',1);
hold on
errorbar(sAVG.BIN.aMUA.pc.coll_layers.data(:,2,L),sAVG.BIN.aMUA.pc.coll_layers.error(:,2,L),'o','marker','none','color','k');
%bar([sAVG.BIN.aMUA.pc_QSM.coll_layers.data(:,2,L);NaN;sAVG.BIN.aMUA.pc_QSM.coll_layers.data(:,3,L)],0.5,'FaceAlpha',.05,'linestyle',':','EdgeColor','k','LineWidth',1.5);
bar(test2(:,2,L),0.5,'FaceAlpha',.05,'linestyle',':','EdgeColor','k','LineWidth',1.5);
set(gca,'box','off','linewidth',1);
ylim([0 65]);
hold off
    if L == 1
        title('Supragranular','FontSize',12);
        xticklabels(labels)
        xlabel('stimulus contrast','FontSize',12)
        ylabel(sprintf('MUA response\nrelative to baseline'),'FontSize',12);
    elseif L == 2
            title('Granular','FontSize',12);
            xticklabels([])
            yticklabels([])
        else 
            title('Infragranular','FontSize',12);
            xticklabels([])
            yticklabels([])
            lgd = legend('Dioptic','QSM','location','northwest');
            lgd.FontSize = 8;
    end
    
end
    
labels = {'22|45', '22|90', '45|22', '45|90','90|22','90|45',[],'22|45','22|90','45|22','45|90','90|22','90|45'}; format bank;

clear L c
for L = 1:3
subplot(2,3,L+3)
% bar([sAVG.DI.aMUA.pc_LSM.coll_layers.data(:,2,L);NaN;sAVG.DI.aMUA.pc_LSM.coll_layers.data(:,3,L)],0.8,'FaceColor',[1, 1, 1],'linestyle','--','EdgeColor','k','LineWidth',0.8);
% hold on
%bar([sAVG.DI.aMUA.pc.coll_layers.data(:,2,L);NaN;sAVG.DI.aMUA.pc.coll_layers.data(:,3,L)],0.6,'FaceColor',[.4 .7 .4],'EdgeColor','k','LineWidth',1);
bar(sAVG.DI.aMUA.pc.coll_layers.data(:,2,L),0.7,'FaceColor',[.4 .7 .4],'EdgeColor','k','LineWidth',1);
hold on
%bar([sAVG.DI.aMUA.pc_QSM.coll_layers.data(:,2,L);NaN;sAVG.DI.aMUA.pc_QSM.coll_layers.data(:,3,L)],0.4,'FaceAlpha',.05,'linestyle',':','EdgeColor','k','LineWidth',1.5);
bar(sAVG.DI.aMUA.pc_QSM.coll_layers.data(:,2,L),0.5,'FaceAlpha',.05,'linestyle',':','EdgeColor','k','LineWidth',1.5);
set(gca,'box','off','linewidth',1);
ylim([0 65]);
hold off
    if L == 1
        ylabel(sprintf('MUA response\nrelative to baseline'),'FontSize',12);
        xticklabels(labels)
        xtickangle(45)
        xlabel('stimulus contrast','FontSize',12)
    elseif L == 2
            xticklabels([])
            yticklabels([])
        else 
            xticklabels([])
            yticklabels([])
            lgd = legend('Dichoptic','QSM','location','northwest');
            lgd.FontSize = 8;
    end
hold off
end


cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\VSS\')
saveas(gcf, strcat('Panel_2_fig3_v2', '.svg'));

%% Panel 3: figure 1
cDE = 4;
cNDE = 3;
% 1 = DE22NDE45, 2 = DE22NDE90, 3 = DE45NDE22, 4 = DE45NDE90
% 5 = DE90NDE22, 6 = 4;3
di = 6;

switch cDE
    case 2
        DELevel = '.22';
    case 3
        DELevel = '.45';
    case 4
        DELevel = '.90';
end

switch cNDE
    case 2
        NDELevel = '.22';
    case 3
        NDELevel = '.45';
    case 4
        NDELevel = '.90';
end

figure('position',[151,58.333,834.66,574.66]);
subplot(1,3,1)
plot(squeeze(sAVG.DE.aMUA.pc.coll.aligned(cDE,2,:)),corticaldepth,'Color',[0, 0.4470, 0.7410],'linewidth',1.5);
hold on 
plot(squeeze(sAVG.NDE.aMUA.pc.coll.aligned(cNDE,2,:)),corticaldepth,'Color', [0, 0.4470, 0.7410],'linewidth',.8);
grid on
xlim([-5 80]);
%hline(0,'-.','BOL4')
ylim([-7 17])
set(gca,'box','off','linewidth',1);
yticklabels({'-0.5','0','0.5','1','1.5','2','2.5'})
xlabel('MUA (percent change)','FontSize',12);
title('Monocular responses','FontSize',12);
legend(sprintf('%s DE',DELevel),sprintf('%s NDE',NDELevel),'Location','northeast','orientation','vertical');
ylabel('Depth (mm) relative to layer 4/5 boundary','FontSize',12);
hold off

subplot(1,3,2)
plot(squeeze(sAVG.DI.aMUA.pc.coll.aligned(di,2,:)),corticaldepth,'Color',[.4 .8 .4],'linewidth',1.5);
hold on 
plot(squeeze(sAVG.DI.aMUA.pc_NRM.coll.aligned(di,2,:)),corticaldepth,'linewidth',1.5,'linestyle',':','color','k');
plot(squeeze(sAVG.DI.aMUA.pc_LSM.coll.aligned(di,2,:)),corticaldepth,'linewidth',1.5,'linestyle','--','color',[.35 .4 .3]);
grid on
xlim([-5 80]);
%hline(0,'-.','BOL4')
ylim([-7 17])
yticklabels([]);
set(gca,'box','off','linewidth',1);
xlabel('MUA (percent change)','FontSize',12);
title(sprintf('Binocular response \nand Model predictions'),'FontSize',12);
legend('BIN','NRM','LSM','Location','northeast','orientation','vertical');
hold off

subplot(1,3,3)
plot(squeeze(sAVG.DI.aMUA.pc_NRM.coll.aligned(di,2,:)-(sAVG.DI.aMUA.pc.coll.aligned(di,2,:))),corticaldepth,'linewidth',1.5,'linestyle',':','color','k');
hold on 
plot(squeeze(sAVG.DI.aMUA.pc_LSM.coll.aligned(di,2,:)-(sAVG.DI.aMUA.pc.coll.aligned(di,2,:))),corticaldepth,'linewidth',1.5,'linestyle','--','color',[.35 .4 .3]);
grid on
xlim([-25 25]);
%hline(0,'-.k')
ylim([-7 17])
vl = vline(0, 'k');
set(vl,'linewidth',1);
xlabel('MUA (difference)','FontSize',12)
yticklabels([]);
set(gca,'box','off','linewidth',1);
title(sprintf('Model difference \nfrom Binocular response'),'FontSize',12);
legend('NRM','LSM','Location','northeast','orientation','vertical');
hold off

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\VSS\')
saveas(gcf, strcat('panel_3_fig1_v2', '.svg'));

%% Panel 3: figure 2

figure('position',[259.6666666666666,41.666666666666664,449.3333333333334,599.3333333333333]);
% 2 = [.22], 3 = [.45], 4 = [.90]
cDE = 4;
cNDE = 3;
% 1 = DE22NDE45, 2 = DE22NDE90, 3 = DE45NDE22, 4 = DE45NDE90
% 5 = DE90NDE22, 6 = 4;3
di = 6;

switch cDE
    case 2
        DELevel = '.22';
    case 3
        DELevel = '.45';
    case 4
        DELevel = '.90';
end

switch cNDE
    case 2
        NDELevel = '.22';
    case 3
        NDELevel = '.45';
    case 4
        NDELevel = '.90';
end

clear i L
for L = 1:3
subplot(3,1,L)
plot(smooth(sAVG.DI.aMUA.pc_NRM.layers.data(di,:,L)-sAVG.DI.aMUA.pc.layers.data(di,:,L),0.1,'rloess'),':','linewidth',1.5,'color','k');
hold on
plot(smooth(sAVG.DI.aMUA.pc_LSM.layers.data(di,:,L)-sAVG.DI.aMUA.pc.layers.data(di,:,L),0.1,'rloess'),'--','linewidth',1.5,'color','k');
ylim([-25 25])
xlim([0 600])
set(gca,'linewidth',1,'box','off')
hline(0,'k')
   if L == 1
        title('Supragranular','FontSize',12);
        xticklabels([]);
        yticklabels([]);
    elseif L == 2
            title('Granular','FontSize',12);
            xticklabels([]);
            yticklabels([]);
        else 
            title('Infragranular','FontSize',12);
            lgd = legend('NRM residual','LSM residual','location','southeast');
            lgd.FontSize = 8;
            xlabel('time (ms)','FontSize',12);
            ylabel(sprintf('model error'),'FontSize',12);
    end
end

%sgtitle(sprintf('Model error across stimulus duration \nDichoptic: %s DE | %s NDE',DELevel,NDELevel));

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\VSS\')
saveas(gcf, strcat('panel_3_fig2_v2', '.svg'));


%% Panel 3: figure 3
labels = {0, 22, 45, 90,[],0,22,45,90}; format bank;
figure('Position', [140,189,784.3333333333333,436]);

test = sAVG.BIN.aMUA.pc.coll_layers.data(:,:,:)-sAVG.BIN.aMUA.pc.coll_layers.data(1,:,:);
test2 = sAVG.BIN.aMUA.pc_NRM.coll_layers.data(:,:,:)-sAVG.BIN.aMUA.pc_NRM.coll_layers.data(1,:,:);
test3 = sAVG.BIN.aMUA.pc_LSM.coll_layers.data(:,:,:)-sAVG.BIN.aMUA.pc_LSM.coll_layers.data(1,:,:);

clear L c
for L = 1:3
subplot(2,3,L)
%bar([sAVG.BIN.aMUA.pc_LSM.coll_layers.data(:,2,L);NaN;sAVG.BIN.aMUA.pc_LSM.coll_layers.data(:,3,L)],0.8,'FaceColor',[1, 1, 1],'linestyle','--','EdgeColor','k','LineWidth',0.8);
%hold on
%bar([sAVG.BIN.aMUA.pc.coll_layers.data(:,2,L);NaN;sAVG.BIN.aMUA.pc.coll_layers.data(:,3,L)],0.6,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',1);
bar(test(:,2,L),.6,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',1);
hold on
%bar([sAVG.BIN.aMUA.pc_QSM.coll_layers.data(:,2,L);NaN;sAVG.BIN.aMUA.pc_QSM.coll_layers.data(:,3,L)],0.5,'FaceAlpha',.05,'linestyle',':','EdgeColor','k','LineWidth',1.5);
bar(test2(:,2,L),0.4,'FaceAlpha',.05,'linestyle',':','EdgeColor','k','LineWidth',1.3);
bar(test3(:,2,L),0.8,'FaceAlpha',.05,'linestyle','--','EdgeColor','k','LineWidth',1.3);
set(gca,'box','off','linewidth',1);
ylim([0 65]);
hold off
    if L == 1
        title('Supragranular','FontSize',12);
        xticklabels(labels)
        xlabel('stimulus contrast','FontSize',12)
        ylabel(sprintf('MUA response\nrelative to baseline'),'FontSize',12);
    elseif L == 2
            title('Granular','FontSize',12);
            xticklabels([])
            yticklabels([])
        else 
            title('Infragranular','FontSize',12);
            xticklabels([])
            yticklabels([])
            lgd = legend('Dioptic','NRM','LSM','location','northwest');
            lgd.FontSize = 7;
    end
    
end
    
labels = {'22|45', '22|90', '45|22', '45|90','90|22','90|45',[],'22|45','22|90','45|22','45|90','90|22','90|45'}; format bank;

clear L c
for L = 1:3
subplot(2,3,L+3)
% bar([sAVG.DI.aMUA.pc_LSM.coll_layers.data(:,2,L);NaN;sAVG.DI.aMUA.pc_LSM.coll_layers.data(:,3,L)],0.8,'FaceColor',[1, 1, 1],'linestyle','--','EdgeColor','k','LineWidth',0.8);
% hold on
%bar([sAVG.DI.aMUA.pc.coll_layers.data(:,2,L);NaN;sAVG.DI.aMUA.pc.coll_layers.data(:,3,L)],0.6,'FaceColor',[.4 .7 .4],'EdgeColor','k','LineWidth',1);
bar(sAVG.DI.aMUA.pc.coll_layers.data(:,2,L),0.6,'FaceColor',[.4 .7 .4],'EdgeColor','k','LineWidth',1);
hold on
%bar([sAVG.DI.aMUA.pc_QSM.coll_layers.data(:,2,L);NaN;sAVG.DI.aMUA.pc_QSM.coll_layers.data(:,3,L)],0.4,'FaceAlpha',.05,'linestyle',':','EdgeColor','k','LineWidth',1.5);
bar(sAVG.DI.aMUA.pc_NRM.coll_layers.data(:,2,L),0.4,'FaceAlpha',.05,'linestyle',':','EdgeColor','k','LineWidth',1.3);
bar(sAVG.DI.aMUA.pc_LSM.coll_layers.data(:,2,L),0.8,'FaceAlpha',.05,'linestyle','--','EdgeColor','k','LineWidth',1.3);
set(gca,'box','off','linewidth',1);
ylim([0 65]);
hold off
    if L == 1
        ylabel(sprintf('MUA response\nrelative to baseline'),'FontSize',12);
        xticklabels(labels)
        xtickangle(45)
        xlabel('stimulus contrast','FontSize',12)
    elseif L == 2
            xticklabels([])
            yticklabels([])
        else 
            xticklabels([])
            yticklabels([])
            lgd = legend('Dichoptic','NRM','LSM','location','northwest');
            lgd.FontSize = 7;
    end
hold off
end

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\VSS\')
saveas(gcf, strcat('panel_3_fig3_v2', '.svg'));