%% find unique units analyze colorflicker and bwflicker files:
clear all
warning off

addpath('/volumes/drobo/users/kacie/code/processdata_code');
addpath('/volumes/drobo/users/kacie/code/nbanalysis');
addpath('/Users/kaciedougherty/Documents/Code/Spike Sorting with KiloSort/KiloSort Utils');

saveon       = 1;
use_evcodes  = 0;
savedrname   = '/Users/kaciedougherty/Documents/analyzeddata/LGNrespprop';
PRE          = -50;          % ms relative to stim onset
POST         = 1000;         % ms relative to stim onset
cmap         = lines;
colors       = fliplr(cmap(1:6,:));

subj = 'I';
[gdates,gfiles,stpath,kiloflick] = getgroupedData;
for fg = 1:length(gdates)
    files = gfiles{fg};
    files(find(cellfun('isempty',(strfind(files,'bwflicker'))))) = [];
    depth = sprintf('_d%02s',num2str(sum(~cellfun('isempty',strfind(gdates(1:fg-1),gdates{fg}))) + 1 ));
    gfiles{fg} = files;
    gdepths{fg} = depth;
    clear files depth
end
unqdates = unique(gdates,'stable');
%%
for days = 1:length(unqdates)
    
    fprintf('\n%s\n',char(unqdates(days)));
    ct = 0;
    for fg = 1:length(gdates)
        if ~strcmp(gdates{fg},unqdates(days))
            continue
        else
            
            if ~isempty(gfiles{fg}) & ~isempty(kiloflick{fg}.chan)
                for fl = 1
                    
                    clearvars -except ct WAVES c figH fg days unqdates fl gdates gfiles stpath kiloflick gdepths colors PRE POST savedrname use_evcodes saveon subj needfiles VISCHANS
                    BRdatafile   = char(strcat(gdates{fg},'_',subj,'_',gfiles{fg}(fl)));
                    brdrname     = sprintf('/volumes/Toshiba External/%s',BRdatafile(1:8));
                    mldrname     = brdrname;
                    kilofname    = strcat('/volumes/Toshiba External/flickerkilo/',BRdatafile(1:8),gdepths{fg},'_kilo_ss.mat');
                    
                    if ~exist(kilofname) | ~exist(strcat(mldrname,'/',BRdatafile,'.nev'))
                        needfiles{fg,fl} = BRdatafile;
                        continue
                    else
                        
                        load(kilofname);
                        clusterMap = sss.clusterMap;
                        chans      = kiloflick{fg}.chan;
                        [~,id]     = ismember(sss.clusterMap(:,2),chans);
                        chans      = chans(id(id>0));
                        clusters   = clusterMap(id>0,1);
                        
                        fn  = fieldnames(sss);
                        if strfind(BRdatafile,'bwflicker')
                            fd  = find(~cellfun(@isempty,strfind(fn, 'bwflicker')));
                        else strfind(BRdatafile,'colorflicker')
                            fd  = find(~cellfun(@isempty,strfind(fn, 'colorflicker')));
                        end
                        
                        SPKFL = sss.(fn{fd});
                        if fg ==1 || fg > 1 & ~strcmp(gdates{fg},gdates{fg-1})
                            figure(1),
                            close
                            figH = figure(1);
                            c = 1;
                        else
                            
                            figH = figure(1);
                            c = c +1;
                        end
                        
                        for ch = 1:length(clusters)
                            SPK = SPKFL.spikeTimes(SPKFL.spikeClusters == clusters(ch));
                            waveform = mean(SPKFL.spikeWaves(chans(ch),:,SPKFL.spikeClusters == clusters(ch)),3);
                            wavevar  = std(SPKFL.spikeWaves(chans(ch),:,SPKFL.spikeClusters == clusters(ch)),0,3)./(sqrt(length(find(SPKFL.spikeClusters == clusters(ch)))));
                            waves(:,ch) = waveform;
                            vwaves(:,ch) = wavevar;
                        end
                        
                        ct = ct + 1;
                        WAVES(days).file(ct).wave = waves;
                        WAVES(days).file(ct).wavevar = vwaves;
                        
                        clear waves vwaves;
                    end
                end
            end
        end
    end
end
close all
%%
for d = 1:length(WAVES)
    
    if ~isempty(WAVES(d).file)
        clear allwaves
        fls = length(WAVES(d).file);
        
        ct = 0;
        for f = 1:fls
            eachfl = WAVES(d).file(f).wave;
            for s = 1: size(eachfl,2)
                ct = ct + 1;
                allwaves(:,ct) = WAVES(d).file(f).wave(:,s);
            end
        end
        nunits = size(allwaves,2);
        if nunits > 1
            
            C = nchoosek([1:nunits],2);
            div = divisors(size(C,1));
            row = div(round(median([1:length(div)])));
            col = size(C,1)./row;
            figure,set(gcf,'Position',[1 1 1200 800]);
            ct = 0;
            for s = 1:size(C,1)
                ct = ct + 1;
                subplot(row,col,ct);
                plot(allwaves(:,C(s,1)),'LineWidth',2);
                hold on;
                plot(allwaves(:,C(s,2)),'LineWidth',2);
                hold on;
                axis tight;
                if ct == 1
                    title(gca,strcat('u1:',num2str(C(s,1)),' u2:',num2str(C(s,2))));
                else
                    title(gca,strcat('u1:',num2str(C(s,1)),' u2:',num2str(C(s,2))));
                end
            end
        end
    end
end
