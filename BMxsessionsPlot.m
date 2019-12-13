%% BMplot_all
%

figure('position',[213.6666666666667,149.6666666666667,724.6666666666666,425.3333333333334]);

clear c L
for L = 1:3
    subplot(2,3,L)
    plot(tmp.STIM.refwin,sAVG.DE.aMUA.pc.layers.data(:,:,L),'color','b')
    hold on
    plot(tmp.STIM.refwin,sAVG.BIN.aMUA.pc.layers.data(:,:,L),'color','r')
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
    plot(tmp.STIM.refwin,((smooth(sAVG.BIN.aMUA.pc.layers.data(c,:,L)))-(smooth(sAVG.DE.aMUA.pc.layers.data(c,:,L)))),'linewidth',.1);
    hold on
    %ylimit = max(abs(get(gcf,'ylim')));
    ylimit = 50;
    set(gca,'ylim',[-10 ylimit],'Box','off','TickDir','out')
    ylabel({'Percent change'...
    'from baseline'});
    xlabel('time (ms)');
    end
end

sgtitle({'V1 laminae contrast responses: monoptic (blue) vs dioptic (red) stimulation'},'Interpreter','none');