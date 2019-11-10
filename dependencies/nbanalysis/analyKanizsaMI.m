% analyze laminar V1 data for kanizsa paper:

clear;
fdates = {'140724','140725','140807','140818','140731','140729'}; % '140820'

int    = 0;
saveon = 1;
kanpdfname = '/users/kaciedougherty/documents/analyzeddata/kanizsaMI';

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
    
    plotEVP([fnames.evp],brdrname,extension,el);
    title(gca,BRdatafile(1:6));
    subplot(1,4,1), title(gca,'single trial LFP');
    subplot(1,4,2), title(gca,'avg LFP');
    subplot(1,4,3), title(gca,'avg CSD');
    if ~exist(strcat(kanpdfname,'.pdf'))
        export_fig(kanpdfname,'-pdf','-nocrop');
    else
        export_fig(kanpdfname,'-pdf','-nocrop','-append');
    end
    %%
    if length(fnames.kanizsa) == 2
        % center file one:
        [f_trLFP,f_trTHETA,f_trGAMMA,f_onsets,f_conds,PRE,POST,Fs] = loadKanizsaFile([fnames.kanizsa{1}],brdrname,extension,el);
        [c_trLFP,c_trTHETA,c_trGAMMA,c_onsets,c_conds,PRE,POST,Fs] = loadKanizsaFile([fnames.kanizsa{2}],brdrname,extension,el);
    else
        BRdatafile = [fnames.kanizsa{1}];
        [trLFP,trTHETA,trGAMMA,onsets,conds,PRE,POST,Fs] = loadKanizsaFile([fnames.kanizsa{1}],brdrname,extension,el);
    end
    
    %% concatenate across in contour/in center files:
    if length(fnames.kanizsa) == 2
        conds    = [f_conds c_conds];
        onsets   = [f_onsets c_onsets];
        trLFP    = cat(3,f_trLFP,c_trLFP);
        trTHETA  = cat(3,f_trTHETA ,c_trTHETA);
        trGAMMA  = cat(3,f_trTHETA ,c_trGAMMA);
    end
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
        ANADATA(s).bothcond = MI;
    end
    
    if int == 1
        kMI   = fliplr(rot90(interp2(MI(:,:,1),10,'cubic'),3));
        conMI = fliplr(rot90(interp2(MI(:,:,2),10,'cubic'),3));
    else
        kMI   = fliplr(rot90(MI(:,:,1),3));
        conMI = fliplr(rot90(MI(:,:,2),3));
    end
    
    [~,~,sink] = getSink(BRdatafile(1:6));
    totchan = 24;
    depths     = [fliplr([1:sink-1]*100) 0 : -100: ((totchan-sink)*-100)];
    
    if c == 2
        ANADATA(s).depths = depths;
    end
    
    figure, set(gcf,'Color','w','Position',[1 1 1500 600]);
    subplot(1,2,1)
    imagesc(depths,depths,kMI),
    colormap('jet'),cb1 = colorbar;
    title(gca,strcat('kanizsa 1',' n=',num2str(length(find(conds == 1))),' sess',BRdatafile(1:6)));
    set(gca,'YDir','normal','XDir','reverse');
    xlabel('depth for theta (microns)');
    v = vline(0); h = hline(0); set(v,'LineWidth',2,'Color','w'); set(h,'LineWidth',2,'Color','w');
    ylabel('depth for gamma (microns)');axis square
    
    subplot(1,2,2)
    imagesc(depths,depths,conMI),cb2 = colorbar;
    title(gca,strcat('control 2',' n=',num2str(length(find(conds == 2))),' sess',BRdatafile(1:6)));
    set(gca,'YDir','normal','XDir','reverse');
    first = get(cb1,'Limits'); caxis([first]);
    xlabel('depth for theta (microns)');
    v = vline(0); h = hline(0); set(v,'LineWidth',2,'Color','w'); set(h,'LineWidth',2,'Color','w');
    ylabel('depth for gamma (microns)');axis square
    export_fig(kanpdfname,'-pdf','-nocrop','-append');
    
    if length(fnames.kanizsa) == 2
        %% CENTER ALONE
        clear MI kMI conMI
        for c = 1:2
            clear thesetrs hMI
            thesetrs = find(f_conds == c);
            tic
            [hMI] = runMILFP(f_trTHETA(:,:,thesetrs),f_trGAMMA(:,:,thesetrs),Fs);
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
        
        [~,~,sink] = getSink(BRdatafile(1:6));
        totchan = 24;
        depths     = [fliplr([1:sink-1]*100) 0 : -100: ((totchan-sink)*-100)];
        
        
        figure, set(gcf,'Color','w','Position',[1 1 1500 600]);
        subplot(1,2,1)
        imagesc(depths,depths,kMI),
        colormap('jet'),cb1 = colorbar;
        title(gca,strcat('CENTER kanizsa 1',' n=',num2str(length(find(f_conds == 1))),' sess',BRdatafile(1:6)));
        set(gca,'YDir','normal','XDir','reverse');
        xlabel('depth for theta (microns)');
        v = vline(0); h = hline(0); set(v,'LineWidth',2,'Color','w'); set(h,'LineWidth',2,'Color','w');
        ylabel('depth for gamma (microns)');axis square
        
        subplot(1,2,2)
        imagesc(depths,depths,conMI),cb2 = colorbar;
        title(gca,strcat('CENTER control 2',' n=',num2str(length(find(f_conds == 2))),' sess',BRdatafile(1:6)));
        set(gca,'YDir','normal','XDir','reverse');
        first = get(cb1,'Limits'); caxis([first]);
        xlabel('depth for theta (microns)');
        v = vline(0); h = hline(0); set(v,'LineWidth',2,'Color','w'); set(h,'LineWidth',2,'Color','w');
        ylabel('depth for gamma (microns)');axis square
        export_fig(kanpdfname,'-pdf','-nocrop','-append');
        
        %% CONTOUR ALONE
        clear MI kMI conMI
        for c = 1:2
            clear thesetrs hMI
            thesetrs = find(c_conds == c);
            tic
            [hMI] = runMILFP(c_trTHETA(:,:,thesetrs),c_trGAMMA(:,:,thesetrs),Fs);
            toc
            MI(:,:,c) = squeeze(hMI);
        end
        if c == 2
            ANADATA(s).contour = MI;
        end
        
        if int == 1
            kMI    = fliplr(rot90(interp2(MI(:,:,1),10,'cubic'),3));
            conMI  = fliplr(rot90(interp2(MI(:,:,2),10,'cubic'),3));
        else
            kMI    = fliplr(rot90(MI(:,:,1),3));
            conMI  = fliplr(rot90(MI(:,:,2),3));
        end
        
        [~,~,sink] = getSink(BRdatafile(1:6));
        totchan    = 24;
        depths     = [fliplr([1:sink-1]*100) 0 : -100: ((totchan-sink)*-100)];
        
        
        figure, set(gcf,'Color','w','Position',[1 1 1500 600]);
        subplot(1,2,1)
        imagesc(depths,depths,kMI),
        colormap('jet'),cb1 = colorbar;
        title(gca,strcat('CONTOUR kanizsa 1',' n=',num2str(length(find(c_conds == 1))),' sess',BRdatafile(1:6)));
        set(gca,'YDir','normal','XDir','reverse');
        xlabel('depth for theta (microns)');
        v = vline(0); h = hline(0); set(v,'LineWidth',2,'Color','w'); set(h,'LineWidth',2,'Color','w');
        ylabel('depth for gamma (microns)');axis square
        
        subplot(1,2,2)
        imagesc(depths,depths,conMI),cb2 = colorbar;
        title(gca,strcat('CONTOUR control 2',' n=',num2str(length(find(c_conds == 2))),' sess',BRdatafile(1:6)));
        set(gca,'YDir','normal','XDir','reverse');
        first = get(cb1,'Limits'); caxis([first]);
        xlabel('depth for theta (microns)');
        v = vline(0); h = hline(0); set(v,'LineWidth',2,'Color','w'); set(h,'LineWidth',2,'Color','w');
        ylabel('depth for gamma (microns)');axis square
        export_fig(kanpdfname,'-pdf','-nocrop','-append');
        
    end
end

%%

    clear tlab emparray  empdepths datadepths plotdepths mat
    
    totchan    = size(ANADATA(1).bothcond,1);
    cond       = size(ANADATA(1).bothcond,3);
    emparray   = nan(totchan*2,totchan*2,cond,length(ANADATA));
    
    empdepths    = [2400:-100:-2400];
    nsess = 0;
    for s = 1:length(ANADATA)
        
        
        [~,~,sink] = getSink(fdates{s});
        datadepths    = [fliplr([1:sink-1]*100) 0 : -100: ((totchan-sink)*-100)];
        
        [~,dataid]    = ismember(empdepths,datadepths);
        
       if isempty(ANADATA(s).center)
           
            emparray(find(dataid),find(dataid),:,s)   = ANADATA(s).bothcond;
       else
            emparray(find(dataid),find(dataid),:,s)   = ANADATA(s).center;
       end
        nsess = nsess + 1;
    end
    tlabel = 'center'; 

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
    
    

