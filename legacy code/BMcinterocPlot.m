%% BMplot
% Single Session plots

if ~isexist('STIM')
    error('No session loaded')
end

%% Dominant Eye Contrast responses (Contrasts)
% Dominant eye response to stimulus with different contrast levels

figure('position',[15,135,1200,500]);
subplot(1,length(STIM.levels),numel(STIM.levels))
hold on
contrastValue = max(STIM.levels);
global scalingfactor
f_ShadedLinePlotbyDepth(mean(STIM.aMUA.raw(:,:,STIM.conditions.DE(numel(STIM.levels),:)),3),0:(1/(numel(STIM.channels))):1,STIM.refwin,STIM.channels,1,1);
title('.90 contrast in DE');
xlabel('time (ms)');
hold off

clear i 
for i = 1:length(STIM.levels)-1
    subplot(1,length(STIM.levels),i);
    f_ShadedLinePlotbyDepth_BAM(mean(STIM.aMUA.raw(:,:,STIM.conditions.DE(i,:)),3),0:(1/(numel(STIM.channels))):1,STIM.refwin,STIM.channels,1,1,false,scalingfactor);
   
plot([0 0], ylim,'k')
plot([STIM.off STIM.off], ylim,'k')
if i == 1
    title({STIM.levels(i),' contrast in both eyes'});
else 
    title({STIM.levels(i),' contrast in DE'});
end
xlabel('time (ms)')
ylabel('contacts indexed down from surface')
hold off
end

sgtitle({'aMUA | Varying contrast to dominant eye',BRdatafile},'Interpreter','none');

% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_contrasts-DE',BRdatafile), '-jpg', '-transparent');

%% Bar plots of Monocular vs Binocular Stimulation (Bar-contrasts)
% Several bar plots for select contacts showing monocular (blue bar) and
% binocular stimulation (red bar).

figure('Position', [60 211 1100 300]);
selectchannels = [STIM.laminae.supra(1):3:STIM.laminae.infra(end)];
%selectchannels = [19 20 21 22 23 24 25 26];
clear c
for c = 1:length(selectchannels)
    subplot(1,length(selectchannels),c)
    bar(STIM.BIN.aMUA.pc.coll(:,2,selectchannels(c)),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
    hold on
    bar(STIM.DE.aMUA.pc.coll(:,2,selectchannels(c)),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
    set(gca,'box','off');
    ylim([-10 100]);
    xticklabels('')
    xlabel('contrast level')
    title({'Contact',selectchannels(c)});
hold off

    if c == 1
    ylabel('percent change from baseline')
    end
end

sgtitle({'Monocular vs Binocular responses'...
    'Collapsed across transient window (40-100ms)',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_bar-contrasts-transient',BRdatafile), '-jpg', '-transparent');

%% All contacts, collapsed across time, sectioned by contrast level (Tightplot)
figure('position',[185,41.66,645.3333333333333,599.3333333333333]);

[ha, pos] = tight_subplot(2,3,[0.005 .03],[.10 .2],[.05 .05]); %channels, columns, [spacing], [bottom and top margin], [left and right margin]
clear c
for c = 1:3
    
    axes(ha(c)); % ha is a variable that gets the axis of each subplot
    plot(flipud(squeeze(STIM.DE.aMUA.pc.coll(c+1,2,:))),STIM.channels,'b','linewidth',.5);
    hold on
    plot(flipud(squeeze(STIM.BIN.aMUA.pc.coll(c+1,2,:))),STIM.channels,'.-r','linewidth',0.5);
    hline(STIM.channels(end)-STIM.laminae.gran(end),'-.')
    xlim([-10 100])
    yticks('')
    yticklabels(fliplr(1:length(STIM.channels)))
    ylim([1 length(STIM.channels)])
    %yticklabels({flipud(1:length(STIM.channels))});
    grid on
    hold off
    xticklabels('');
    %xlabel('Percent change');
    
    if c == 1
        ylabel('Transient');
        title({STIM.levels(c+1),'contrast'});
    else 
        title({STIM.levels(c+1),'contrast'});
    end   
   
end

clear c
for c = 4:6
    
    axes(ha(c)); % ha is a variable that gets the axis of each subplot
    plot(flipud(squeeze(STIM.DE.aMUA.pc.coll(c-2,3,:))),STIM.channels,'b','linewidth',.5);
    hold on
    plot(flipud(squeeze(STIM.BIN.aMUA.pc.coll(c-2,3,:))),STIM.channels,'.-r','linewidth',0.5);
    hline(STIM.channels(end)-STIM.laminae.gran(end),'-.')
    xlim([-10 100])
    yticks('')
    yticklabels(fliplr(1:length(STIM.channels)))
    ylim([1 length(STIM.channels)])
    %yticklabels({flipud(1:length(STIM.channels))});
    grid on
    hold off
    
    xlabel('Percent change'); 
    if c == 4
        ylabel('Sustained');
    else 
   
    end
end

sgtitle({'Contrast response profiles | Ocularity x contrast x time',BRdatafile},'interpreter','none')
%cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
%export_fig(sprintf('%s_tightplot',BRdatafile), '-jpg', '-transparent');

%% Collapsed contrast responses for monocular and binocular stim (contrast lines)
% Designed to show Monocular vs binocular contrast responses alongside CSD
% for layer reference. 
figure('position',[185 150 887 450]);

h = subplot(1,4,1);
bAVG_iCSD = filterCSD(STIM.CSD.bsl')';
imagesc(STIM.refwin,STIM.channels,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(0); set(v,'color','k','linestyle','-','linewidth',1);
set(gca,'tickdir','out');  
climit = max(abs(get(gca,'CLim'))*1);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
plot([STIM.off STIM.off], ylim,'k','linestyle','-.','linewidth',0.5)
title('iCSD: All Trials')
xlabel('time (ms)')
clrbar = colorbar; %clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts indexed down from surface');
set(h,'position',[0.065388951521984,0.097526988745119,0.145749605022835,0.722586325702473]);
hold off

subplot(1,4,3)
plot(fliplr(squeeze(STIM.BIN.aMUA.pc.coll(:,2,:))),STIM.channels);
hold on
%ylabel('contacts indexed down from surface');
set(gca,'box','off');
%yticklabels({'','20','15','10','5'});
xlabel('Percent change');
grid on
climit = max(abs(get(gca,'xlim'))*1);
xlim([-5 climit*1.2]);
yticks(1:length(STIM.channels))
yticklabels(fliplr(1:length(STIM.channels)))
ylim([1 length(STIM.channels)])
title('Binocular');
hold off

subplot(1,4,2)
hold on
plot(fliplr(squeeze(STIM.DE.aMUA.pc.coll(:,2,:))),STIM.channels);
set(gca,'box','off');
grid on
xlim([-5 climit*1.2]);
yticks(1:length(STIM.channels))
yticklabels(fliplr(1:length(STIM.channels)))
ylim([1 length(STIM.channels)])
xlabel('Percent change');
grid on
title('Monocular');
hold off

subplot(1,4,4)
plot(fliplr(squeeze(STIM.calc.aMUA.pc.subtract.BIN_DE.coll(:,2,:))),STIM.channels);
hold on 
grid on
xlim([-5 25]);
set(gca,'box','off');
yticks(1:length(STIM.channels))
yticklabels(fliplr(1:length(STIM.channels)))
ylim([1 length(STIM.channels)])
xlabel('Percent change');
title('Subtraction (bin - mon)');
%legend(num2str(STIM.levels),'Location','southoutside','orientation','horizontal');
hold off

sgtitle({'Monocular vs binocular aMUA averaged over stimulus duration',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_lineplots_transient',BRdatafile), '-jpg', '-transparent');


%% Contrast lines across time (Timeplots)
% and corresponding binocular minus monocular subtractions 

figure('position',[213.6666666666667,149.6666666666667,724.6666666666666,425.3333333333334]);

clear c L
for L = 1:3
    subplot(2,3,L)
    plot(STIM.refwin,STIM.DE.aMUA.pc.layers(:,:,L),'color','b')
    hold on
    plot(STIM.refwin,STIM.BIN.aMUA.pc.layers(:,:,L),'color','r')
    %ylimit = max(abs(get(gcf,'ylim')));
    ylimit = 100;
    set(gca,'ylim',[-10 ylimit],'Box','off','TickDir','out')
    ylabel({'Percent change'...
    'from baseline'});
    xlabel('time (ms)');
    if L == 1
        title('Supragranular');
    elseif L == 2
            title('Granular');
        else 
            title('Infragranular');
    end
end

clear c L
for L = 1:3
    for c = 1:4
    subplot(2,3,L+3)
    plot(STIM.refwin,((smooth(STIM.BIN.aMUA.pc.layers(c,:,L)))-(smooth(STIM.DE.aMUA.pc.layers(c,:,L)))),'linewidth',.1);
    hold on
    %ylimit = max(abs(get(gcf,'ylim')));
    ylimit = 50;
    set(gca,'ylim',[-10 ylimit],'Box','off','TickDir','out')
    ylabel({'Percent change'...
    'from baseline'});
    xlabel('time (ms)');
    end
end

sgtitle({'V1 laminae contrast responses: monoptic (blue) vs dioptic (red) stimulation'...
   ,BRdatafile},'Interpreter','none');

% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_timeplots',BRdatafile), '-jpg', '-transparent');

%% Bar Graph: Coll_layer; DE vs NDE vs BIN
% Using bar graphs to represent monocular vs binocular response layer
% differences and fold change from one to the other

rcontrast = round(STIM.levels,2,'significant');

figure('Position', [155,98,965,487]);

labels = {0, 22, 45, 90,[],0,22,45,90}; format bank;

clear L c
for L = 1:3
subplot(2,3,L)
x1 = [STIM.BIN.aMUA.pc.coll_layers(:,2,L);NaN;STIM.BIN.aMUA.pc.coll_layers(:,3,L)];
x2 = [STIM.DE.aMUA.pc.coll_layers(:,2,L);NaN;STIM.DE.aMUA.pc.coll_layers(:,3,L)];
bar(x1,.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(x2,0.6,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
%bar([STIM.NDE.aMUA.pc.coll_layers(:,2,L);NaN; STIM.NDE.aMUA.pc.coll_layers(:,3,L)],0.4,'FaceColor',[.2, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-5 55]);
xticklabels(labels)
xlabel('% Contrast');
ylabel('Relative aMUA response');
    if L == 1
        title('Supragranular');
        legend('BIN','DE','NDE','location','northeast');
    elseif L == 2
            title('Granular');
        else 
            title('Infragranular');
    end
hold off
end

clear L c
for L = 1:3
subplot(2,3,L+3)
bar([STIM.calc.aMUA.pc.subtract.BIN_DE.coll_layers(:,2,L);NaN; STIM.calc.aMUA.pc.subtract.BIN_DE.coll_layers(:,3,L)],.8,'FaceColor',[0.2, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
hold on
%bar([STIM.calc.aMUA.pc.subtract.BIN_NDE.coll_layers(:,2,L);NaN; STIM.calc.aMUA.pc.subtract.BIN_NDE.coll_layers(:,3,L)],0.6,'FaceColor',[0.2500, 0.250, 0.2980],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-5 30]);
xticklabels(labels)
xlabel('contrast')
ylabel('difference');
hold off
end

sgtitle({'Binned contacts by layer | aMUA responses',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_BINvsMON',BRdatafile), '-jpg', '-transparent');

%% Bar Graphs: BIN vs Models

figure('Position', [155,41.666666666666664,994.6666666666665,599.3333333333333]);

clear L c
for L = 1:3
subplot(3,3,L)
bar([STIM.DE.aMUA.pc.coll_layers(:,2,L);NaN;STIM.DE.aMUA.pc.coll_layers(:,3,L)],0.8,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
hold on
bar([STIM.NDE.aMUA.pc.coll_layers(:,2,L);NaN;STIM.NDE.aMUA.pc.coll_layers(:,3,L)],0.6,'FaceColor',[.2, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-5 80]);
xticklabels(labels)
xlabel('contrast')
ylabel('aMUA response');
    if L == 1
        title('Supragranular');
        lgd = legend('DE','NDE','location','northwest');
        lgd.FontSize = 4;
    elseif L == 2
            title('Granular');
        else 
            title('Infragranular');
    end
hold off
end

clear L c
for L = 1:3
subplot(3,3,L+3)
bar([STIM.BIN.aMUA.pc_LSM.coll_layers(:,2,L);NaN;STIM.BIN.aMUA.pc_LSM.coll_layers(:,3,L)],0.8,'FaceColor',[1, 1, 1],'linestyle','--','EdgeColor','k','LineWidth',0.8);
hold on
bar([STIM.BIN.aMUA.pc.coll_layers(:,2,L);NaN;STIM.BIN.aMUA.pc.coll_layers(:,3,L)],0.6,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
bar([STIM.BIN.aMUA.pc_QSM.coll_layers(:,2,L);NaN;STIM.BIN.aMUA.pc_QSM.coll_layers(:,3,L)],0.4,'FaceColor',[.85, .325, .098],'linestyle',':','EdgeColor','k','LineWidth',1);
set(gca,'box','off');
ylim([-5 80]);
xticklabels(labels)
xlabel('contrast')
ylabel('aMUA response');
hold off
    if L == 1
        lgd = legend('LSM','BIN','QSM','location','northwest');
        lgd.FontSize = 4;
    end
end

clear L
for L = 1:3
subplot(3,3,L+6)
bar([STIM.BIN.aMUA.pc_LSM.coll_layers(:,2,L)-STIM.BIN.aMUA.pc.coll_layers(:,2,L);NaN;STIM.BIN.aMUA.pc_LSM.coll_layers(:,3,L)-STIM.BIN.aMUA.pc.coll_layers(:,3,L)],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',0.8,'linestyle','--');
hold on
bar([STIM.BIN.aMUA.pc_QSM.coll_layers(:,2,L)-STIM.BIN.aMUA.pc.coll_layers(:,2,L);NaN;STIM.BIN.aMUA.pc_QSM.coll_layers(:,3,L)-STIM.BIN.aMUA.pc.coll_layers(:,3,L)],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',1,'linestyle',':');
set(gca,'box','off');
ylim([-15 30]);
xticklabels(labels);
xlabel('contrast')
ylabel('difference from model');
hold off
end

%sgtitle({'Binned contacts by layer | aMUA responses',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_BINvsModels',BRdatafile), '-jpg', '-transparent');

%% Bar graphs: DIvsModels
figure('Position', [155,98,965,487]);
labels = {'22|45', '22|90', '45|22', '45|90','90|22','90|45',[],'22|45','22|90','45|22','45|90','90|22','90|45'}; format bank;
clear L c
for L = 1:3
subplot(2,3,L)
bar([STIM.DI.aMUA.pc_LSM.coll_layers(:,2,L);NaN;STIM.DI.aMUA.pc_LSM.coll_layers(:,3,L)],0.8,'FaceColor',[1, 1, 1],'linestyle','--','EdgeColor','k','LineWidth',0.8);
hold on
bar([STIM.DI.aMUA.pc.coll_layers(:,2,L);NaN;STIM.DI.aMUA.pc.coll_layers(:,3,L)],0.6,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
bar([STIM.DI.aMUA.pc_QSM.coll_layers(:,2,L);NaN;STIM.DI.aMUA.pc_QSM.coll_layers(:,3,L)],0.4,'FaceAlpha',.05,'linestyle',':','EdgeColor','k','LineWidth',1);
set(gca,'box','off');
ylim([-5 65]);
xticklabels(labels)
xlabel('contrast')
ylabel('aMUA response');
hold off
xtickangle(45)
    if L == 1
        title('Supragranular');
        lgd = legend('LSM','BIN','QSM','location','northwest');
        lgd.FontSize = 5;
    elseif L == 2
            title('Granular');
        else 
            title('Infragranular');
    end
end

clear L
for L = 1:3
subplot(2,3,L+3)
bar([STIM.DI.aMUA.pc_LSM.coll_layers(:,2,L)-STIM.DI.aMUA.pc.coll_layers(:,2,L);NaN;STIM.DI.aMUA.pc_LSM.coll_layers(:,3,L)-STIM.DI.aMUA.pc.coll_layers(:,3,L)],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',0.8,'linestyle','--');
hold on
bar([STIM.DI.aMUA.pc_QSM.coll_layers(:,2,L)-STIM.DI.aMUA.pc.coll_layers(:,2,L);NaN;STIM.DI.aMUA.pc_QSM.coll_layers(:,3,L)-STIM.DI.aMUA.pc.coll_layers(:,3,L)],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',1,'linestyle',':');
set(gca,'box','off');
ylim([-15 25]);
xticklabels(labels);
xtickangle(45)
xlabel('contrast')
ylabel('difference from model');
    if L == 1
        lgd = legend('LSM-BIN','QSM-BIN','location','northwest');
        lgd.FontSize = 5;
    end
hold off
end

sgtitle({'Binned contacts by layer | aMUA responses',BRdatafile},'Interpreter','none');

% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_DIvsModels',BRdatafile), '-jpg', '-transparent');

%% Interesting plots

%% Models and Data over time

figure('position',[145,47,984.6666666666665,582.6666666666665]);
clear i L
for L = 1:3
subplot(3,3,L)
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc.layers(2,:,L),0.1,'rloess'),'-k','linewidth',2);
hold on
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_LSM.layers(2,:,L),0.1,'rloess'),'-g','linewidth',1);
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_QSM.layers(2,:,L),0.1,'rloess'),'-b','linewidth',1);
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_NRM.layers(2,:,L),0.1,'rloess'),'-r','linewidth',1);
ylim([-5 85])
xlim([-50 550])
set(gca,'box','off','linewidth',1,'FontSize',12);
    if L == 1
        title('Supragranular');
%                 lgd = legend('Data','LSM','QSM','AVE','location','northeast');
%                 lgd.FontSize = 6;
                ylabel(sprintf('response\n (percent change)'));
                xlabel('time (ms)');
                %xticklabels([]);
                %yticklabels([]);
    elseif L == 2
            title('Granular');
                xticklabels([]);
                yticklabels([]);
        else 
            title('Infragranular');
            xticklabels([]);
            yticklabels([]);
    end
end

clear i L
for L = 1:3
subplot(3,3,L+3)
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc.layers(3,:,L),0.1,'rloess'),'-k','linewidth',2);
hold on
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_LSM.layers(3,:,L),0.1,'rloess'),'-g','linewidth',1);
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_QSM.layers(3,:,L),0.1,'rloess'),'-b','linewidth',1);
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_NRM.layers(3,:,L),0.1,'rloess'),'-r','linewidth',1);
ylim([-5 85])
xlim([-50 550])
set(gca,'box','off','linewidth',1,'FontSize',12);
xticklabels([]);
yticklabels([]);
end

clear i L
for L = 1:3
subplot(3,3,L+6)
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc.layers(4,:,L),0.1,'rloess'),'-k','linewidth',2);
hold on
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_LSM.layers(4,:,L),0.1,'rloess'),'-g','linewidth',1);
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_QSM.layers(4,:,L),0.1,'rloess'),'-b','linewidth',1);
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_NRM.layers(4,:,L),0.1,'rloess'),'-r','linewidth',1);
ylim([-5 85])
xlim([-50 550])
set(gca,'box','off','linewidth',1,'FontSize',12);
xticklabels([]);
yticklabels([]);
end

%sgtitle("Model difference from Binocular response (midline = no difference)");

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
saveas(gcf,strcat(sprintf('%s_models_over-time-Allori',BRdatafile),'.svg'));

%% Model difference over time

figure('position',[145,47,984.6666666666665,582.6666666666665]);
clear i L
for L = 1:3
subplot(3,3,L)
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_LSM.layers(2,:,L)-STIM.BIN.aMUA.pc.layers(2,:,L),0.1,'rloess'),'-g','linewidth',1);
hold on
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_QSM.layers(2,:,L)-STIM.BIN.aMUA.pc.layers(2,:,L),0.1,'rloess'),'-b','linewidth',1);
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_NRM.layers(2,:,L)-STIM.BIN.aMUA.pc.layers(2,:,L),0.1,'rloess'),'-r','linewidth',1);
ylim([-25 50])
xlim([-50 550])
set(gca,'box','off','linewidth',1,'FontSize',12);
hl = hline(0,'k');
set(hl,'linewidth',1);
    if L == 1
        title('Supragranular');
%                 lgd = legend('LSM','QSM','AVE','location','northeast');
%                 lgd.FontSize = 6;
                ylabel(sprintf('model residual'));
                xlabel('time (ms)');
                %xticklabels([]);
                %yticklabels([]);
    elseif L == 2
            title('Granular');
                xticklabels([]);
                yticklabels([]);
        else 
            title('Infragranular');
            xticklabels([]);
            yticklabels([]);
    end
end

clear i L
for L = 1:3
subplot(3,3,L+3)
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_LSM.layers(3,:,L)-STIM.BIN.aMUA.pc.layers(3,:,L),0.1,'rloess'),'-g','linewidth',1);
hold on
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_QSM.layers(3,:,L)-STIM.BIN.aMUA.pc.layers(3,:,L),0.1,'rloess'),'-b','linewidth',1);
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_NRM.layers(3,:,L)-STIM.BIN.aMUA.pc.layers(3,:,L),0.1,'rloess'),'-r','linewidth',1);
ylim([-25 50])
xlim([-50 550])
set(gca,'box','off','linewidth',1,'FontSize',12);
hl = hline(0,'k');
set(hl,'linewidth',1.5);
xticklabels([]);
yticklabels([]);
%     if L == 3
%         ylim([-10 30])
%     end
end

clear i L
for L = 1:3
subplot(3,3,L+6)
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_LSM.layers(4,:,L)-STIM.BIN.aMUA.pc.layers(4,:,L),0.1,'rloess'),'-g','linewidth',1);
hold on
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_QSM.layers(4,:,L)-STIM.BIN.aMUA.pc.layers(4,:,L),0.1,'rloess'),'-b','linewidth',1);
plot(STIM.refwin,smooth(STIM.BIN.aMUA.pc_NRM.layers(4,:,L)-STIM.BIN.aMUA.pc.layers(4,:,L),0.1,'rloess'),'-r','linewidth',1);
ylim([-25 50])
xlim([-50 550])
set(gca,'box','off','linewidth',1,'FontSize',12);
hl = hline(0,'k');
set(hl,'linewidth',1.5);
xticklabels([]);
yticklabels([]);
%     if L == 3
%         ylim([-10 30])
%     end
end

%sgtitle("Model difference from Binocular response (midline = no difference)");

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
saveas(gcf,strcat(sprintf('%s_model_diff_overtime-Allori',BRdatafile),'.svg'));

%% Model difference over space

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
plot(flipud(squeeze(STIM.DE.aMUA.pc.coll(cIndex,2,:))),STIM.channels,'linewidth',1.5,'Color',[0, 0.4470, 0.7410]);
hold on 
plot(flipud(squeeze(STIM.NDE.aMUA.pc.coll(cIndex,2,:))),STIM.channels,'linewidth',.8','Color', [0, 0.4470, 0.7410]);
grid on
set(gca,'box','off','linewidth',1);
xlim([-5 80]);
%hline(0,':','BOL4')
%ylim([-7 17])
%yticklabels({'-0.5','0','0.5','1','1.5','2','2.5'})
xlabel('MUA (percent change)','FontSize',12);
legend(sprintf('%s DE',cLevel),sprintf('%s NDE',cLevel),'Location','northeast','orientation','vertical');
title('Monocular responses','FontSize',12);
ylabel('Depth (mm) relative to layer 4/5 boundary','FontSize',12);
hold off

subplot(1,3,2)
plot(squeeze(STIM.BIN.aMUA.pc.coll(cIndex,2,:)),STIM.channels,'Color','r','linewidth',1.5);
hold on 
plot(squeeze(STIM.BIN.aMUA.pc_LSM.coll(cIndex,2,:)),STIM.channels,'linewidth',1.5,'linestyle','--','color',[.35 .4 .3]);
plot(squeeze(STIM.BIN.aMUA.pc_QSM.coll(cIndex,2,:)),STIM.channels,'linewidth',1.5,'linestyle',':','color','k');
plot(squeeze(STIM.BIN.aMUA.pc_NRM.coll(cIndex,2,:)),STIM.channels,'linewidth',1.5,'linestyle','-.','color','k');
grid on
xlim([-5 80]);
%xticklabels([]);
xlabel('MUA (percent change)','FontSize',12);
yticklabels([]);
%ylim([-7 17])
set(gca,'box','off','linewidth',1);
title(sprintf('Binocular response \nand Model prediction'),'FontSize',12);
legend('BIN','LSM','QSM','NRM','Location','northeast','orientation','vertical');
hold off

subplot(1,3,3)
plot(squeeze(STIM.BIN.aMUA.pc_LSM.coll(cIndex,2,:)-(STIM.BIN.aMUA.pc.coll(cIndex,2,:))),STIM.channels,'linewidth',1.5,'linestyle','--','color',[.35 .4 .3]);
hold on
plot(squeeze(STIM.BIN.aMUA.pc_QSM.coll(cIndex,2,:)-(STIM.BIN.aMUA.pc.coll(cIndex,2,:))),STIM.channels,'linewidth',1.5,'linestyle',':','color','k');
plot(squeeze(STIM.BIN.aMUA.pc_NRM.coll(cIndex,2,:)-(STIM.BIN.aMUA.pc.coll(cIndex,2,:))),STIM.channels,'linewidth',1.5,'linestyle','-.','color','k');
grid on
xlim([-10 30]);
%ylim([-7 17])
%xticklabels([]);
xlabel('MUA (difference)','FontSize',12)
yticklabels([]);
vl = vline(0, 'k');
set(vl,'linewidth',1);
set(gca,'box','off','linewidth',1);
title(sprintf('Model difference \nfrom Binocular response'),'FontSize',12);
legend('LSM-BIN','QSM-BIN','NRM-BIN','Location','northeast','orientation','vertical');
hold off

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_model-space',BRdatafile), '-jpg', '-transparent');

%% CSD: DE vs BIN
figure('position',[166.3333333333333,85,990.6666666666667,537.3333333333333]);
clear i
for i = 1:4
subplot(2,4,i);
bAVG_iCSD = filterCSD(squeeze(STIM.DE.CSD.bsl(i,:,:))')';
imagesc(STIM.refwin,STIM.channels,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(0); set(v,'color','k','linestyle','-','linewidth',1);
set(gca,'tickdir','out');  
plot([STIM.off STIM.off], ylim,'k','linestyle','-.','linewidth',0.5)
xlabel('time (ms)')
clrbar = colorbar; %clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts indexed down from surface');
set(gca,'CLim',[-500 500],'Box','off','TickDir','out')
%set(h,'position',[0.065388951521984,0.097526988745119,0.145749605022835,0.722586325702473]);
hold off
    if i == 1
        title('0 contrast');
    elseif i == 2
        title('.22 contrast');
    elseif i == 3
        title('.45 contrast');
    else
        title('.90 contrast');
    end
end

for i = 1:4
subplot(2,4,i+4);
bAVG_iCSD = filterCSD(squeeze(STIM.BIN.CSD.bsl(i,:,:))')';
imagesc(STIM.refwin,STIM.channels,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(0); set(v,'color','k','linestyle','-','linewidth',1);
set(gca,'tickdir','out');  
plot([STIM.off STIM.off], ylim,'k','linestyle','-.','linewidth',0.5)
xlabel('time (ms)')
clrbar = colorbar; %clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts indexed down from surface');
set(gca,'CLim',[-500 500],'Box','off','TickDir','out')
%set(h,'position',[0.065388951521984,0.097526988745119,0.145749605022835,0.722586325702473]);
hold off
end
% climit = max(abs(get(gca,'CLim'))*1);
% set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_BINvDE_CSD',BRdatafile), '-jpg', '-transparent');

%% CSD vs models: Dioptic

figure('position',[166.3333333333333,85,990.6666666666667,537.3333333333333]);
clear i
for i = 1:4
subplot(2,4,i);
bAVG_iCSD = filterCSD(squeeze(STIM.BIN.CSD.bsl(i,:,10:30))')';
imagesc(STIM.refwin,10:30,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(0); set(v,'color','k','linestyle','-','linewidth',1);
set(gca,'tickdir','out');  
plot([STIM.off STIM.off], ylim,'k','linestyle','-.','linewidth',0.5)
xlabel('time (ms)')
clrbar = colorbar; %clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts indexed down from surface');
set(gca,'CLim',[-3000 3000],'Box','off','TickDir','out')
%set(h,'position',[0.065388951521984,0.097526988745119,0.145749605022835,0.722586325702473]);
hold off
    if i == 1
        title('0 contrast');
    elseif i == 2
        title('.22 contrast');
    elseif i == 3
        title('.45 contrast');
    else
        title('.90 contrast');
    end
end

for i = 1:4
subplot(2,4,i+4);
bAVG_iCSD = filterCSD(squeeze(STIM.BIN.CSD.LSM(i,:,10:30))')';
imagesc(STIM.refwin,10:30,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(0); set(v,'color','k','linestyle','-','linewidth',1);
set(gca,'tickdir','out');  
plot([STIM.off STIM.off], ylim,'k','linestyle','-.','linewidth',0.5)
xlabel('time (ms)')
clrbar = colorbar; %clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts indexed down from surface');
set(gca,'CLim',[-3000 3000],'Box','off','TickDir','out')
%set(h,'position',[0.065388951521984,0.097526988745119,0.145749605022835,0.722586325702473]);
hold off
end

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_BINvModels_CSD',BRdatafile), '-jpg', '-transparent');

%% CSD vs models: Dichoptic

figure('position',[166.3333333333333,85,990.6666666666667,537.3333333333333]);
clear i
for i = 1:6
subplot(2,6,i);
bAVG_iCSD = filterCSD(squeeze(STIM.DI.CSD.bsl(i,:,10:30))')';
imagesc(STIM.refwin,10:30,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
set(gca,'tickdir','out');  
plot([STIM.off STIM.off], ylim,'k','linestyle','-.','linewidth',0.5)
xlabel('time (ms)')
ylabel('contacts indexed down from surface');
set(gca,'CLim',[-3000 3000],'Box','off','TickDir','out')
%set(h,'position',[0.065388951521984,0.097526988745119,0.145749605022835,0.722586325702473]);
hold off
    if i == 1
        title('22 DE | 45 NDE');
    elseif i == 2
        title('45 DE | 22 NDE');
    elseif i == 3
        title('22 DE | 90 NDE');
    elseif i == 4
        title('90 DE | 22 NDE');
        elseif i == 5
        title('45 DE | 90 NDE');
        elseif i == 6
        title('90 DE | 45 NDE');
        %colorbar; v = vline(0); set(v,'color','k','linestyle','-','linewidth',1);
        %clrbar = colorbar; %clrbar.Label.String = 'nA/mm^3';
        %set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
    end
end

clear i
for i = 1:6
subplot(2,6,i+6);
bAVG_iCSD = filterCSD(squeeze(STIM.DI.CSD.LSM(i,:,10:30))')';
imagesc(STIM.refwin,10:30,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
%colorbar; v = vline(0); set(v,'color','k','linestyle','-','linewidth',1);
set(gca,'tickdir','out');  
plot([STIM.off STIM.off], ylim,'k','linestyle','-.','linewidth',0.5)
xlabel('time (ms)')
%clrbar = colorbar; %clrbar.Label.String = 'nA/mm^3'; 
%set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts indexed down from surface');
set(gca,'CLim',[-3000 3000],'Box','off','TickDir','out')
%set(h,'position',[0.065388951521984,0.097526988745119,0.145749605022835,0.722586325702473]);
hold off
end

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_BINvModels_CSD-dichoptic',BRdatafile), '-jpg', '-transparent');