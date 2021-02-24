%% CRF ocularity
% Figure save directory
if strcmp(getenv('username'),'bmitc')
    figDir = 'C:\Users\bmitc\OneDrive - Vanderbilt\Maier Lab\Figures\';
elseif strcmp(getenv('username'),'mitchba2')
    figDir = 'C:\Users\mitchba2\OneDrive - Vanderbilt\Maier Lab\Figures\';
elseif strcmp(getenv('username'),'bmitc_000')
    figDir = 'C:\Users\bmitc_000\OneDrive - Vanderbilt\Maier Lab\Figures\';
end

% toggle to save figures
flag_figsave = 1;

%% Curve Fitting: Dioptic data
clear -global Data
clear bin mon bin clear prd prd1 prd2 cv curves a K n b 

% choices
tw = 1; % time window
error = 1;

% strings
curves = {'de','nde','bin'};

try
x = IDX(1).monLevels*100; x(1) = x(1)+1;
catch
    disp('Using average of binned contrast levels');
    x = mean(IDX(1).cbins,2)*100; x(1) = x(1)+1;
end

% data
for o = 1:length(occGroups)
de_all(o).units = OCC.MON.DE_PS(o).RESP;
nde_all(o).units = OCC.MON.NDE_PS(o).RESP;
bin_all(o).units = OCC.BIN.PS(o).RESP;
end

for o = 1:length(occGroups)
    de_all(o).avg = nanmean(de_all(o).units,3);
    de_all(o).err = nanstd(de_all(o).units,[],3)./sqrt(size(de_all(o).units,3));
    nde_all(o).avg = nanmean(nde_all(o).units,3);
    nde_all(o).err = nanstd(nde_all(o).units,[],3)./sqrt(size(nde_all(o).units,3));
    bin_all(o).avg = nanmean(bin_all(o).units,3);
    bin_all(o).err = nanstd(bin_all(o).units,[],3)./sqrt(size(bin_all(o).units,3));
end

global Data
Data(1,:) = x;
figure('position',[275,342,281,507]);
for o = 1:length(occGroups) % for each group of ocularity
    
    de = de_all(o).avg(:,tw);
    de_ref = de_all(o).avg(:,1);
    nde = nde_all(o).avg(:,tw);
    bin = bin_all(o).avg(:,tw);
    
    mn       = min(de_ref);
    mx       = max(de_ref);
    nDE      = (de - mn)./(mx - mn);
    nNDE     = (nde - mn)./(mx - mn);
    nBIN     = (bin - mn)./(mx - mn);
    
    for curve = 1:size(curves,2) % for each curve (bin and mon)
        switch curves{curve}
            case 'de'
                Data(2,:) = nDE;
            case 'nde'
                Data(2,:) = nNDE;
            case 'bin'
                Data(2,:) = nBIN;
        end
        
        [a,K,n,~] = BMfitCRdata;
        predictions = 1:100;
        for c = 1:length(predictions) % generate prediction for mon
            prd(curve,c) = a*[(c^n)/((c^n) + (K^n))]; % mon prediction
        end
    end
    % Plot
    subplot(length(occGroups),1,o)
    s1 = semilogx(predictions,prd(1,:),'color','k','linestyle','-','linewidth',2,'HandleVisibility','off');  hold on
    s2 = semilogx(predictions,prd(2,:),'color',[0 0 0]+0.9,'linestyle','-','linewidth',1.75,'HandleVisibility','off');  hold on
    s3 = semilogx(predictions,prd(3,:),'color','b','linestyle','-','linewidth',2,'HandleVisibility','off');  hold on
    
    s4 = semilogx(Data(1,:),nDE,'o','color','k','linewidth',2,'markersize',5);
    %s5 = semilogx(Data(1,:),nNDE,'o','color',[.4 .4 .4],'linewidth',2,'markersize',5);
    s6 = semilogx(Data(1,:),nBIN,'*','color','b','linewidth',2,'markersize',5);
    
    if tw == 1
        set(s1,'linestyle','-'); ...set(s2,'linestyle','-');
            set(s3,'linestyle','-');
    else
        set(s1,'linestyle','-'); ...set(s2,'linestyle','--');
            set(s3,'linestyle','-');
    end
    
    set(gca,'FontSize',12,'linewidth',1.5,'box','off',...
        'xlim',[1 100],...
        'XTick',x);
    if tw == 1
        set(gca,'ylim',[0 1.4]);
    else
        set(gca,'ylim',[0 0.7]);
    end

    if length(occGroups) == 3
        if o == 1
            %occRange = abs([occValues(occGroups{o,1}(1)),occValues(occGroups{o,1}(end))]);
            %title(sprintf('Lowest (n = %d)\n OccIndex [%.3g : %.3g]',occLengths(o),occRange(1),occRange(2)),'FontSize',16);
            %xlabel('Stimulus contrast'); ylabel(sprintf('Relative response \n (normalized to monocular)'));
            xticklabels([]); ylabel([]); xticklabels([]);
        elseif o == 2
            %occRange = abs([occValues(occGroups{o,1}(1)),occValues(occGroups{o,1}(end))]);
            %title(sprintf('Average (n = %d)\n OccIndex [%.3g : %.3g]',occLengths(o),occRange(1),occRange(2)),'FontSize',16);
            %mCt = sum(Trls.mon(:,1,:),'all'); bCt = sum(Trls.bin(2:end,1,:),'all');
            %ndeCt = sum(Trls.mon(:,2,:),'all');
            %legend(sprintf('DE (%d trials)',mCt),...sprintf('NDE (%d trials)',ndeCt)
                %sprintf('BIN (%d trials)',bCt),'location','northwest'); legend boxoff
            xlabel([]); ylabel([]); xticklabels([]);
        else
            %occRange = abs([occValues(occGroups{o,1}(1)),occValues(occGroups{o,1}(end))]);
            %title(sprintf('Highest (n = %d)\n OccIndex [%.3g : %.3g]',occLengths(o),occRange(1),occRange(2)),'FontSize',16);
            %xlabel([]); ylabel([]); xticklabels([]);
                    xticks([1,5,10,20,50,100]);
                    ylabel('Normalized response','fontsize',12);
                    xlabel('contrast','fontsize',12);
        end
    end

end

if flag_figsave == 1
    cd(strcat(figDir,'ocularity\'));
    saveas(gcf, strcat('CRF-dioptic_', num2str(tw),'.svg'));
    disp("Figure saved");
else
    disp("Figure was not saved");
end

%% Dichoptic
clear -global Data
clear DE NDE BIN DI de nde bin di prd1 prd2 cv curves a K n b
clear nMon nmedNDE nhighNDE nlowNDE

% choices
tw = 1; %1 = transient, 2 = sustained

% strings
curves = {'zeroNDE','lowNDE','medNDE','highNDE'};

% data
try
    x = IDX(1).monLevels*100; x(1) = x(1)+1;
catch
    disp('Not currently set up for this analysis');
end

for o = 1:3
DI(o).units = LAY.DI.PS(o).RESP;
end

for o = 1:3
    DI(o).avg = nanmean(DI(o).units,3);
    DI(o).err = nanstd(DI(o).units,[],3)./sqrt(size(DI(o).units,3));
end

global Data
Data(1,:) = x;
figure('position',[135,404.3333333333333,1420.666666666667,432.6666666666667]);
for o = 1:3
    di  = DI(o).avg(:,tw);
    de_ref = DI(o).avg(1:4,1); %
    
    %     Vertical: colors are NDE strength
%     Horizontal: colors are DE strength
% %         00|00  00|22  00|45  00|90 - yel
% %         22|00  22|22  22|45  22|90 - oran
% %         45|00  45|22  45|45  45|90 - red
% %         90|00  90|22  90|45  90|90 - purp
% %         -yel   -oran  -red   -purp
    
    zeroNDE     = di(1:4); % curve 1
    lowNDE      = di(5:8); % curve 2
    medNDE      = di(9:12); % curve 3
    highNDE     = di(13:16); % curve 4

    mn             = min(de_ref); 
    mx             = max(de_ref); 
    nzeroNDE       = (zeroNDE - mn)./(mx - mn); % curve 1
    nlowNDE        = (lowNDE - mn)./(mx - mn); % curve 2
    nmedNDE        = (medNDE - mn)./(mx - mn); % curve 3
    nhighNDE       = (highNDE - mn)./(mx - mn); % curve 4

    for curve = 1:size(curves,2)
        switch curves{curve}
            case 'zeroNDE'
                Data(2,:) = nzeroNDE;
            case 'lowNDE'
                Data(2,:) = nlowNDE;
            case 'medNDE'
                Data(2,:) = nmedNDE;
            case 'highNDE'
                Data(2,:) = nhighNDE;
        end
 
        [a, K, n, ~] = BMfitCRdata;
        predictions = 1:100;
            for c = 1:length(predictions)
                prd(curve,c) = a*[(c^n)/((c^n) + (K^n))]; % prediction
            end
    end

% Plot

subplot(1,3,o)
semilogx(predictions,prd(1,:),'color',[0.9290 0.6940 0.1250],'linestyle','-','linewidth',1.7,'HandleVisibility','off');  hold on
semilogx(predictions,prd(2,:),'color',[0.8500 0.3250 0.0980],'linestyle','-','linewidth',1.7,'HandleVisibility','off');  
semilogx(predictions,prd(3,:),'color',[0.6350 0.0780 0.1840],'linestyle','-','linewidth',1.7,'HandleVisibility','off');  
semilogx(predictions,prd(4,:),'color',[0.4940 0.1840 0.5560],'linestyle','-','linewidth',1.7,'HandleVisibility','off');  

semilogx(Data(1,:),nzeroNDE,'o','color',[0.9290 0.6940 0.1250],'linewidth',1.7,'markersize',5);
semilogx(Data(1,:),nlowNDE,'o','color',[0.8500 0.3250 0.0980],'linewidth',1.7,'markersize',5);
semilogx(Data(1,:),nmedNDE,'o','color',[0.6350 0.0780 0.1840],'linewidth',1.7,'markersize',5);
semilogx(Data(1,:),nhighNDE,'o','color',[0.4940 0.1840 0.5560],'linewidth',1.7,'markersize',5);

set(gca,'FontSize',12,'linewidth',1.5,'box','off',...
    'ylim',[0 max(nzeroNDE)],... %'XScale','log',...
    'xlim',[1 100],...
     'XTick',x);

    if o == 1
        occRange = abs([occValues(occGroups{o,1}(1)),occValues(occGroups{o,1}(end))]);
            title(sprintf('OccIndex [%.3g : %.3g]',occRange(1),occRange(2)),'FontSize',16);
            xlabel('Stimulus contrast'); ylabel(sprintf('Normalized Multi-unit \n response'));
            xticklabels({x});
  
    elseif o == 2 
        occRange = abs([occValues(occGroups{o,1}(1)),occValues(occGroups{o,1}(end))]);
            title(sprintf('OccIndex [%.3g : %.3g]',occRange(1),occRange(2)),'FontSize',16);
            %mCt = sum(Trls.mon(:,1,:),'all'); bCt = sum(Trls.bin(:,1,:),'all');
            %legend(sprintf('MON (%d trials)',mCt),sprintf('BIN (%d trials)',bCt),'location','northwest'); legend boxoff
            xlabel([]); ylabel([]); xticklabels([]);
    else
        occRange = abs([occValues(occGroups{o,1}(1)),occValues(occGroups{o,1}(end))]);
            title(sprintf('OccIndex [%.3g : %.3g]',occRange(1),occRange(2)),'FontSize',16);
            xlabel([]); ylabel([]); xticklabels([]);
    end
end

if flag_figsave == 1
    cd(strcat(figDir,'ocularity\'));
    saveas(gcf, strcat('CRF-dichoptic_',num2str(tw), '.svg'));
    sprintf("Figure saved");
else
    sprintf("Figure was not saved");
end

