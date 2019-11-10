% main.mn
% analyze our first months or DUAL UPROBE recordings
% 02/16- AVM, KD 3/7

clear all, close all

% updates to main_all (what makes it different from main):
% add unit dim to spk structure (before only 1 unit saved per chan)
% code that plots visual responses for each channel or subset of units
% divided waveform amplitudes by 4 to convert to microVolts
% plots csd


fname  = '121116_H_MLflash002';
saveon = 1; 
spath  = strcat('/volumes/drobo/data/neurophys/doughek/interelect/figures/',fname(1:8),'/'); if saveon == 1, mkdir(spath); end
nowaveplots = 1; % close waveform figures for units
pre = -100; post = 200; tvec = [pre:post];

% HOUSEKEEPING:
%
% INVENTORY OF ALL GOOD SESSIONS (n=42)
[~,~,~,~,~,alldates] = goodsessionlist('130303'); % enter one date to get out all session dates
monks     = [];
distances = [];
for i = 1:length(alldates)
    
    if str2num(alldates{i}) < 130411
        m = 'H';
    else
        m = 'B';
    end
    
    [~,~,~,~,dist] = goodsessionlist(alldates{i});
    distances = [distances dist]; % each element is interelectrode distance for alldates{i}
    monks     = strcat(monks,m);  % each element is subject (H or B) for alldates{i}
    
end
% remove sessions with unknown interleectrode distance for now:
distances(find(isnan(distances))) = [];
monks(find(isnan(distances))) = [];
alldates(find(isnan(distances))) = [];

%%
monk  = fname(8);
% LOAD LFP
lfp_filename = sprintf('/volumes/drobo/data/neurophys/doughek/interelect/stimfiles/%s/%s.lfp',monk,fname);
LFP = load(lfp_filename,'-MAT');
% processed LFP--filename is the same without "_ss."
% directory: /volumes/drobo/data/neurophys/doughek/interelect/restfiles/H
% and /volumes/drobo/data/neurophys/doughek/interelect/restfiles/B

% sort LFP:
[ID,ids] = sortElectrodeLabels(LFP.NeuralLabels);
lfpA = LFP.LFP(:,ids(ID(1),:));
lfpB = LFP.LFP(:,ids(ID(2),:)); clear LFP;

% % % LOAD MUAa
% % mua_filename = sprintf('/volumes/drobo/data/neurophys/doughek/interelect/stimfiles/%s/%s.mua',monk,fname);
% % MUAa = load(mua_filename,'-MAT');
% % % processed LFP--filename is the same without "_ss."
% % % directory: /volumes/drobo/data/neurophys/doughek/interelect/restfiles/H
% % % and /volumes/drobo/data/neurophys/doughek/interelect/restfiles/B


%%
% LOAD SPIKE SORTED DATA
%
filename = sprintf('/Lab Data/BOSS_sortedspikes/%s_ss',fname); % added "_ss" to make distinct from NEV files collected on recording day
ssNEV = openNEV(strcat(filename,'.nev'),'nomat','nosave');
% ssNEV.Data.Spikes
% timestamp of spike is paired with electrode label in field "Electrode".
% These electrode numbers aren't sorted .... code below to address that:
% also if there is only one "unit" for an eelctrode channel, it is MUA (some
% waveforms  look  better than others).
% If there's more than one, one is MUA and others single units

% get basic info about recorded data
neural = ~strcmp('E',{ssNEV.ElectrodesInfo.ConnectorBank}) & ~cellfun('isempty',{ssNEV.ElectrodesInfo.ConnectorBank}); % bank E is the BNCs on the front of the NSP
N.electrodes = length(neural);
N.neural = sum( neural);
N.analog = sum(~neural);

%get labels
NeuralLabels = {ssNEV.ElectrodesInfo(neural).ElectrodeLabel};
if strfind(NeuralLabels{1}','Channel')
    % load a file to get mapping between electrode ID and electrode label
    mfilename = '/Lab Data/BOSS_sortedspikes/121116_H_darkrest002_ss'; % added "_ss" to make distinct from NEV files collected on recording day
    mapNEV = openNEV(strcat(mfilename,'.nev'),'nomat','nosave');
    NeuralLabels = {mapNEV.ElectrodesInfo(~cellfun('isempty',{ssNEV.ElectrodesInfo.ElectrodeID})).ElectrodeLabel};
end
for e = 1:length(NeuralLabels)
    NeuralLabels{e} = NeuralLabels{e}';
end
NeuralInfo = ssNEV.ElectrodesInfo(neural);

[ID,ids] = sortElectrodeLabels(NeuralLabels);
% ID 1, 2, 3, 4 corresponds to banks A, B , C, D
% ids(ID(1),:) is electrode number that pairs with 'e(A)01','e(A)02', ...

%%
sdate = fname(1:6);
if monk == 'H'
    brdir = {'/volumes/drobo/data/neurophys/coxma'};
elseif monk == 'B'
    
    brdir = {'/volumes/drobo/data/neurophys/coxma';...
        '/volumes/drobo/data/neurophys/rig021'};
end

ex = exist(strcat(brdir{1}, '/', sdate,'_', monk));
if any(ex)
    drname = strcat(brdir{1}, '/', sdate,'_', monk,'/');
else
    drname = strcat(brdir{2}, '/', sdate,'_', monk,'/');
end

digNEV = openNEV(strcat(drname,'/',fname,'.nev'),'nomat','nosave');
EventCodes = digNEV.Data.SerialDigitalIO.UnparsedData-128;
EventTimes = digNEV.Data.SerialDigitalIO.TimeStampSec;

if isempty(EventCodes)
    % use photodiode signal to detect onset of flashes
    load(sprintf('/volumes/drobo/data/neurophys/doughek/interelect/stimfiles/%s/%s.bnc',monk,fname),'-MAT');
    thres  = mean(BNC(:,1),1) + (2*std(BNC(:,1),0,1));
    photoT = diff(BNC(:,1)< -thres) > 0;
    onsetT = find(photoT == 1);
    evon = onsetT(diff(onsetT)>1000);
    %                     if length(onsetT) > 10
    %                         pid = pid + 1;
    %                         randT = randi(length(onsetT)-3);
    %                         time = onsetT(randT) - 200 : onsetT(randT + 2) + 200;
    %                         figure,
    %                         set(gcf,'Color','w');
    %                         plot(time,BNC(time,1),'k');
    %                         vline(onsetT(randT:randT+2),'r');
    %                         set(gca,'LineWidth',2,'TickDir','out','Box','off');
    %                         xlabel('time (ms)');
    %                         randT = randi(length(onsetT)-3);
    %                         time = onsetT(randT) - 200 : onsetT(randT + 2) + 200;
    %                     end
    
else
    evon = floor(EventTimes(EventCodes == 23).*1000); % stim on times (ms)
end

%%
for probe = 1:length(ID)
    
    for e = 1:length(ids)
        
        % find spikes for chan
        chan = ids(ID(probe),e);
        chanspks=find(ssNEV.Data.Spikes.Electrode == chan);
        
        %find out how many units there are for chan
        numunits=unique(ssNEV.Data.Spikes.Unit(chanspks));
        numunits = numunits(numunits>0);
        if any(numunits>6)
            fprintf('ignoring units numbered greater than 6\n');
            numunits = numunits(numunits<=6);
        end
        cunit = 0;
        
        for unit = numunits
            unitspks=find(ssNEV.Data.Spikes.Unit == unit);
            if any(unitspks)
                cunit = cunit + 1;
                chanunitspks=intersect(chanspks,unitspks);
                
                tstamps = ssNEV.Data.Spikes.TimeStamp(chanunitspks)./(ssNEV.MetaTags.SampleRes./1000);
                wvform  = ssNEV.Data.Spikes.Waveform(:,chanunitspks)./4; % divide by 4 if input argument "uV" is not used in openNEV call
                
                spks{chan}(cunit).ts = tstamps;
                spks{chan}(cunit).wf = wvform;
                
                
                figure
                % plot mean spike waveform and 95% confidence interval
                mnwf=mean(spks{chan}(cunit).wf,2);
                plot(mnwf,'k')
                hold on
                %compute 95% confidence interval (=1.96xSEM)
                spkstd=std(double(spks{chan}(cunit).wf),0,2); %note: var needs double as input
                spksem=spkstd./sqrt(size(spks{chan}(cunit).wf,2));
                pspkci=spkstd;%1.96*spksem; % PLOT STD
                plot(mnwf+pspkci,'--k'),plot(mnwf-pspkci,'--k'),
                axis tight
                xlabel('time (samples)')
                ylabel('amplitude (uV)')
                chanlabel=chan;
                unitlabel=unit;
                figurelabel=sprintf('chan=%i unit=%i',chanlabel,unitlabel);
                title(figurelabel)
                if nowaveplots == 1
                    close 
                end
                try
                    %compute trough-to-peak (T2P):
                    wfmin = min(mnwf); trough=find(mnwf==wfmin); % trough
                    %usually two peaks: one pre and one post trough (here we want the latter)
                    premax  = max(mnwf(1:trough)); prepeak = find(mnwf==premax);
                    postmax = max(mnwf(trough:end)); postpeak = find(mnwf==postmax);
                    t2p = (postpeak-trough)/30; % T2P in ms (assuming 30kHz sampling)
                    %compute peak amplitude asymmetry (PAA):
                    paa = (postmax-premax)/(postmax+premax);
                    %compute half width
                    hmax=wfmin/2; % define half max (i.e., min) of spike deflection
                    lmst=min(find(mnwf(1:trough) <= hmax));
                    rmst=trough+(min(find(mnwf(trough:end) >= hmax)));
                    fwhm=(rmst-lmst)/30; %Half Width in ms
                    %create 3D vector
                    spks{chan}(cunit).spk_id = [t2p paa fwhm];
                end
            end
        end
        
    end
end

% PC vs IN (see Harris Neuron paper, SF3)
%
% CHECK FOR CLUSTERS WITHIN 3D PLOT OF:
% Y: PEAK AMPLITUDE ASYMMETRY
% X: TROUGH-TO-PEAK (ms)
% Z: HALF WIDTH (ms)

kmat=[];chrec = []; unitrec = []; 
for chan = 1:length(spks)
    for unit = 1:size(spks{chan},2)
        try
            idmat=[spks{chan}(unit).spk_id(1),spks{chan}(unit).spk_id(3),spks{chan}(unit).spk_id(2)];
            kmat=vertcat(kmat,idmat);
            chrec = vertcat(chrec,chan);  % keep track of which ch, unit data points belong to
            unitrec = vertcat(unitrec,unit); 
        end
    end
end

% run k-mans clustering with k=2 clusters
[G,C] = kmeans(kmat,2);
% plot result
clr = lines(2);
figure, hold on
scatter3(kmat(:,1), kmat(:,2), kmat(:,3), 100, clr(G,:), 'Marker','.')
hold off
view(3), axis vis3d, box on, rotate3d on
xlabel('x'), ylabel('y'), zlabel('z')
title(gca,fname,'interpreter','none');
saveas(gcf,'PIN_kmeanscluster','eps2c'); 

% plot average waveform for each cluster
unqG = unique(G);    figure, set(gcf,'Color','w'); figN = get(gcf,'Number'); 
for i = 1:length(unqG)
    
    it = 0; 
    chans = chrec(find(G == unqG(i))); 
    units = unitrec(find(G == unqG(i))); 
 
    for ch = 1:length(chans)
            it = it + 1; 
            avgwf(:,it) = mean(spks{chans(ch)}(units(ch)).wf,2);   
            if abs(min(avgwf(:,it))) < abs(max(avgwf(:,it)))
                avgwf(:,it) = avgwf(:,it).*-1; 
            end
            
            %figure, plot(mean(avgwf,2),'Color',clr(i,:)); 
            
    end
    
    spkstd=std(double(avgwf),0,2); %note: var needs double as input
    spksem=spkstd./sqrt(size(avgwf,2));
    pspkci=spkstd;%1.96*spksem; % PLOT STD

    figure(figN)
    subplot(length(unqG),1,i)
    plot(mean(avgwf,2)); 
    hold on; 
    plot(mean(avgwf,2)+pspkci,'--k'),plot(mean(avgwf,2)-pspkci,'--k'),
    set(gca,'Box','off'); 
    title(gca,strcat(fname,' group ',num2str(i)),'interpreter','none'); 
    clear avgwf
end
saveas(gcf,'avgwaveform_2kgroups','eps2c'); 

% is G 1 interneuron or pyramidal group?
% time to peak, half width, amplitude assymetry
if mean(kmat(G==1,2))<mean(kmat(G==2,2)) && mean(kmat(G==1,3))>mean(kmat(G==2,3))
    g1  = 'IN';g2 = 'P';
elseif mean(kmat(G==1,2))>mean(kmat(G==2,2)) && mean(kmat(G==1,3))<mean(kmat(G==2,3))
    g1 = 'P'; g2 = 'IN';
else
    error('check clusters to see which is P and IN\n');
end

% add cell type label in spks variable

for i = 1:length(unqG)
    
    it = 0; 
    chans = chrec(find(G == unqG(i))); 
    units = unitrec(find(G == unqG(i))); 
 
    for ch = 1:length(chans)
        if i == 1
           spks{chans(ch)}(units(ch)).type  =   g1; 
        else
            spks{chans(ch)}(units(ch)).type =   g2; 
        end
    end
end

%%
% check every neuron's waveform + response
it = 0 ; nunits = 3;
for chans = 1:length(spks)
    
    if ~isempty(spks{chans})
        for unit = 1:size(spks{chans},2)
            if mod(it,nunits) == 0
                it = 0;
            end
            it = it + 1;
            
            if it == 1  
                figure, set(gcf,'Color','w','Position',[1 1 700 500],'PaperPositionMode','auto');
                figN = get(gcf,'Number');
            end
            
            % waveform:
            subplot(nunits,3,(it*3-2))
            % plot mean spike waveform and 95% confidence interval
            mnwf=mean(spks{chans}(unit).wf,2);
            plot(mnwf,'k')
            hold on
            %compute 95% confidence interval (=1.96xSEM)
            spkstd=std(double(spks{chans}(unit).wf),0,2); %note: var needs double as input
            spksem=spkstd./sqrt(size(spks{chans}(unit).wf,2));
            pspkci=spkstd;%1.96*spksem; % PLOT STD
            plot(mnwf+pspkci,'--k'),plot(mnwf-pspkci,'--k'),
            axis tight
            set(gca,'TickDir','out','Box','off');
            xlabel('time (samples)')
            ylabel('amplitude (uV)')
            if chans > 24
                chanlabel=strcat('eB',num2str(chans-24));
            else
                chanlabel=strcat('eA',num2str(chans));
            end
            unitlabel=unit;
            figurelabel=sprintf('chan=%s unit=%i',chanlabel,unitlabel);
            title(figurelabel)
            
            % raster plot showing visual response
            subplot(nunits,3,(it*3-2)+1)
            pre = -50; post = 200; tvec = [pre:post];
            counts = zeros(length(tvec),length(evon));
            for tr = 1:length(evon)
                refwin = evon(tr) + pre:evon(tr)+post;
                sts = tvec(ismember(refwin,spks{chans}(unit).ts));
                plot(sts,repmat(tr,length(sts),1),'o','Color',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerSize',1);
                hold on;
                
                counts(ismember(refwin,spks{chans}(unit).ts),tr) =  counts(ismember(refwin,spks{chans}(unit).ts),tr) + 1;
            end
            ylabel('trial');
            set(gca,'Box','off','TickDir','out');
            ylim([0 length(evon)+1]); xlim([tvec(1) tvec(end)]);
            set(gca,'TickDir','out','Box','off');
            title(gca,fname,'interpreter','none');
            
            % calculate spike density
            subplot(nunits,3,(it*3-2)+2)
            binsz = 10; %st dev of gaussian kernel (ms)
            [convspk] = calcSpikeDensity(counts,binsz);
            sem = std(convspk,0,2)./sqrt(size(convspk,2)); int = 20; 
            h = plot(tvec,mean(convspk,2));
            hold on; 
            errorbar(tvec(1:int:end),mean(convspk(1:int:end,:),2),sem(1:int:end),'LineStyle','none','Color',get(h,'Color')); 
            title(gca,sprintf('spk density std dev %u ms',binsz));
            set(gca,'TickDir','out','Box','off'); xlim([tvec(1) tvec(end)]);
            hold on;
            
            if it == 3 & saveon == 1
                saveas(gcf,strcat(spath,'unitvisualresp_',num2str(figN)),'eps2c');
            end
            
        end
        
    end
end
%% Like Figure 1 in Luczak, Bartho, and Harris (2013)

% trigger response to stimulus onset for nunits 
nunits = 4;  
chans = datasample(find(~cellfun('isempty',spks)),nunits,'replace',false); 
unit = 1; 
figure, set(gcf,'Color','w','position',[1 1 500 600],'PaperPositionMode','auto');
for i = 1:nunits
    subplot(nunits,2,i*2-1)
    counts = zeros(length(tvec),length(evon));
    tvec = pre:post; 
    for tr = 1:length(evon)
        refwin = evon(tr) + pre:evon(tr)+post;
        sts = tvec(ismember(refwin,spks{chans(i)}(unit).ts));
        plot(sts,repmat(tr,length(sts),1),'o','Color',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerSize',1);
        hold on;
        
        counts(ismember(refwin,spks{chans(i)}(unit).ts),tr) =  counts(ismember(refwin,spks{chans(i)}(unit).ts),tr) + 1;
    end
    if chans(i) > 24
        title(gca,sprintf('eB%u',chans(i)-24));
    else
    title(gca,strcat(fname,sprintf('eA%u',chans(i))),'interpreter','none');
    end
    ylabel('trial');
    set(gca,'Box','off','TickDir','out');
    ylim([0 length(evon)+1]); xlim([tvec(1) tvec(end)]);
    hold on;
    figN = get(gcf,'Number');
    
    % calculate spike density 
    binsz = 10; %st dev of gaussian kernel (ms)
    [convspk] = calcSpikeDensity(counts,binsz); 
    subplot(nunits,2,(i*2-1)+1), plot(tvec,mean(convspk,2));
    title(gca,sprintf('spk density, std dev = %d ms',binsz)); 
    xlim([tvec(1) tvec(end)]);
    set(gca,'Box','off','TickDir','out');
    
end
xlabel('t (ms)');
if saveon == 1
    saveas(gcf,strcat(spath,'rndunits_visualresp'),'eps2c'); 
end
% 
% % nunits shuffled between trials (different trial for each neuron) 
% clear trials;
% trials = randperm([length(evon)]);
% trials = trials(1:length(find(~cellfun('isempty',spks))));
% figure, set(gcf,'Color','w','position',[1 1 400 400]); 
% tvec = pre:post;
% it = 0; 
% for i = 1:length(spks)
%     if ~isempty(spks{i})
%         it = it + 1; 
%     tr = trials(it);
%     refwin = evon(tr) + pre:evon(tr)+post;
%     sts = tvec(ismember(refwin,spks{i}(unit).ts));
%     plot(sts,repmat(i,length(sts),1),'o','Color',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerSize',1);
%     hold on;
%     set(gca,'Box','off','TickDir','out');
%    end
% end
% ylabel('neuron'); xlabel('t (ms)'); 
% ylim([0 length(find(~cellfun('isempty',spks)))+1]); xlim([tvec(1) tvec(end)]); 
% title(gca,'neurons shuffled between trials'); 

% all neurons over snippet of time (neurons by time)
figure, set(gcf,'Color','w','position',[1 1 300 450],'PaperPositionMode','auto'); 
ntr = 6; %how many light flashes to include in t snippet
rndtr = randi([1 length(evon)-ntr]); 
extratime = (1000*ntr) + 200; 
refwin = [evon(rndtr) + pre:(evon(rndtr) + extratime)]; 
clear tvec; tvec = [pre:extratime-pre]; 
it = 0; 
for i = 1:length(spks)
    if ~isempty(spks{i})
        for unit = 1:size(spks{i},2)
        it  = it + 1; 
    sts = tvec(ismember(refwin,spks{i}(unit).ts));
    plot(sts,repmat(it,length(sts),1),'o','Color',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerSize',1);
    hold on;
    set(gca,'Box','off','TickDir','out');
        end
    end
end
ylabel('neuron'); xlabel('t (ms)'); 
xlim([tvec(1) tvec(end)]); 
title(gca,strcat(fname,'  example of actual data')); 
if saveon == 1
    saveas(gcf,strcat(spath,'actualdata_allunits'),'eps2c'); 
end


%% laminar analysis 
%% plot CSD to gauge depth alignment

for tr = 1:length(evon)
    refwin = [evon(tr) + pre:evon(tr) + post];
    csdA(:,:,tr) = mod_iCSD(lfpA(refwin,:)')';
    csdB(:,:,tr) = mod_iCSD(lfpB(refwin,:)')';
end

mcsdA = nanmean(csdA,3);
mcsdB = nanmean(csdB,3);

figure
set(gcf,'Color',[1 1 1],'Position',[1 1 700 800],'PaperPositionMode','auto');

subplot(1,2,1)
[fcsd] = filterCSD(mcsdA');
nn = nan(10,size(fcsd,2));
pd_fcsd = [nn ;fcsd; nn]; % pad with nans for lost channels
was = 1:size(lfpA,2);     % n contacts

olddepths = was;
ticks     = linspace(1,size(pd_fcsd,1),length(was));

pr = 1;  %set lims for colorbar
mn = ((min(min(fcsd))).*pr) ;
mx = ((max(max(fcsd))).*pr) ;
tvec = pre:post;
imagesc(tvec,1:size(pd_fcsd,1),pd_fcsd,[mn -mn]); cb = colorbar;
cm =  size(colormap('jet'),1);
dmap = ((-mn)-mn)/cm;
caxis([mn-dmap -mn]); %make sure only nans go white


map = flipud(colormap('jet'));
colormap([1 1 1;map]);
ylim(cb,[mn -mn]);  %set lims for colorbar and add in white for nans
set(cb,'YDir','reverse');
set(gca,'YTick',ticks(2:2:end-1),'YTicklabel',olddepths(2:2:end-1),'FontSize',14); set(gca,'TickDir','out');
xlabel('t (ms) from flash');
set(get(cb,'ylabel'),'string','nA/mm^3','FontSize',14);

t = title(gca,sprintf('%s CSD A',fname),'interpreter','none');
set(t,'FontWeight','bold');
v = vline(0); set(v,'Color',[0.5 0.5 0.5],'LineWidth',1);
h = text(post+80,csd(1,1)-15,sprintf('n = %d trs',length(evon))); set(h,'Color','k','FontSize',14);
ylabel('channel'); 

subplot(1,2,2)
[fcsd] = filterCSD(mcsdB');
nn = nan(10,size(fcsd,2));
pd_fcsd = [nn ;fcsd; nn]; % pad with nans for lost channels
was = 1:size(lfpA,2); % n contacts

olddepths = was;
ticks     = linspace(1,size(pd_fcsd,1),length(was));

pr = 1;  %set lims for colorbar
mn = ((min(min(fcsd))).*pr) ;
mx = ((max(max(fcsd))).*pr) ;
tvec = pre:post;
imagesc(tvec,1:size(pd_fcsd,1),pd_fcsd,[mn -mn]); cb = colorbar;
cm =  size(colormap('jet'),1);
dmap = ((-mn)-mn)/cm;
caxis([mn-dmap -mn]); %make sure only nans go white


map = flipud(colormap('jet'));
colormap([1 1 1;map]);
ylim(cb,[mn -mn]);  %set lims for colorbar and add in white for nans
set(cb,'YDir','reverse');
set(gca,'YTick',ticks(2:2:end-1),'YTicklabel',olddepths(2:2:end-1),'FontSize',14); set(gca,'TickDir','out');
xlabel('t (ms) from flash');
set(get(cb,'ylabel'),'string','nA/mm^3','FontSize',14);

t = title(gca,sprintf('%s CSD B',fname),'interpreter','none');
set(t,'FontWeight','bold');
v = vline(0); set(v,'Color',[0.5 0.5 0.5],'LineWidth',1);
h = text(post+80,csd(1,1)-15,sprintf('n = %d trs',length(evon))); set(h,'Color','k','FontSize',14);
if saveon == 1
    saveas(gcf,strcat(spath,'csd'),'eps2c'); 
end
% 

[~, Asink,Bsink,rank] = goodsessionlist(sdate); 
fprintf('\nassigned sink A:%u sink B:%u rank: %u\n',Asink,Bsink,rank)

%%
% guess laminar compartments for each electrode
supdep = [1:-0.1:0.3]; grandep = [0.2:-.1:0]; infdep = [-.1:-.1:-.5]; 
alldep = [1:-.1:-.5];
supra(1).chan = Asink-10:Asink-3; supra(1).chan = supra(1).chan(supra(1).chan>=1 & supra(1).chan<=24); 
gran(1).chan  = Asink-2:Asink;
infra(1).chan = Asink+1:Asink+5; infra(1).chan = infra(1).chan(infra(1).chan>=1 & infra(1).chan<=24); 

supra(2).chan = Bsink-10:Bsink-3; supra(2).chan = supra(2).chan(supra(2).chan>=1 & supra(2).chan<=24); 
gran(2).chan  = Bsink-2:Bsink;
infra(2).chan = Bsink+1:Bsink+5; infra(2).chan = infra(2).chan(infra(2).chan>=1 & infra(2).chan<=24); 

%%
% Sakata and Harris define sparseness with response probability measure
% (probability that neuron fires 1 time in response to stimulus event)
% look at this as function of cell type and layer
clear RPROB
evpre = 30; evpost = 150; % define event window

thesechans = [supra(1).chan gran(1).chan infra(1).chan]; 
comspks = {spks{thesechans}};
it = 0;
for i = 1:length(comspks)
    if isempty(comspks{i})
        continue
    else
        for unit = 1:size(comspks{i},2)
            it = it + 1;
            counts = zeros(length(tvec),length(evon));
            for tr = 1:length(evon)
                refwin = evon(tr) + evpre:evon(tr)+evpost;
                counts(ismember(refwin,comspks{i}(unit).ts),tr) =  counts(ismember(refwin,comspks{i}(unit).ts),tr) + 1;
            end
            
            % probability firing 1 spike?
            rprob(it) = sum(sum(counts,1)>0)./size(counts,2);
            depth(it) = alldep(i);
      
            if strcmp(comspks{i}(unit).type,'P')
                type(it)  = 0;
            elseif strcmp(comspks{i}(unit).type,'IN')
                type(it)  = 1;
            else
                type(it) = nan;
            end
        end
    end
end


RPROB(1).prob  = rprob; 
RPROB(1).depth = depth;  
RPROB(1).type  = type; clear rprob depth it comspks type

thesechans = [supra(2).chan gran(2).chan infra(2).chan]+24; 
comspks = {spks{thesechans}};
it = 0; 
for i = 1:length(comspks)
    if isempty(comspks{i})
        continue
    else
        for unit = 1:size(comspks{i},2)
            it = it + 1; 
            counts = zeros(length(tvec),length(evon));
            for tr = 1:length(evon)
                refwin = evon(tr) + evpre:evon(tr)+evpost;
                counts(ismember(refwin,comspks{i}(unit).ts),tr) =  counts(ismember(refwin,comspks{i}(unit).ts),tr) + 1;
            end
            
            % probability firing 1 spike?
            rprob(it) = sum(sum(counts,1)>0)./size(counts,2);
            depth(it) = alldep(i);
            if strcmp(comspks{i}(unit).type,'P')
                type(it)  = 0;
            elseif strcmp(comspks{i}(unit).type,'IN')
                type(it)  = 0;
            else
                type(it) = nan;
            end
            
        end
    end 
end

RPROB(2).prob  = rprob;
RPROB(2).depth = depth;
RPROB(2).type  = type; clear rprob depth it comspks type
%%
mxunit = 6; % max number of units on one channel
Pprob  = nan(size(RPROB,2),mxunit,length(alldep)); % electrode by Nunits by depth
INprob = nan(size(RPROB,2),mxunit,length(alldep));
for i = 1:size(RPROB,2)
    if any(RPROB(i).type == 0)
  unqdep = unique(RPROB(i).depth(RPROB(i).type == 0)); 
  for dep = 1:length(unqdep)
      
      did = find(alldep == unqdep(dep)); 
      Pprob(i,1:length(find(RPROB(i).depth == unqdep(dep))),did) =  RPROB(i).prob(RPROB(i).depth == unqdep(dep)); 
        
  end   
    end
  
  clear unqdep dep 
   if any(RPROB(i).type == 1)
  unqdep = unique(RPROB(i).depth(RPROB(i).type == 1)); 
  for dep = 1:length(unqdep)
      
      did = find(alldep == unqdep(dep)); 
      INprob(i,1:length(find(RPROB(i).depth == unqdep(dep))),did) =  RPROB(i).prob(RPROB(i).depth == unqdep(dep)); 
        
  end  
   end
end

avgelect = 1; %average across electrode probes?
if avgelect == 1
    figure, set(gcf,'Color','w','Position',[1 1 500 650]);
    subplot(1,3,[1:2])
    clear pbars inbars bars
    rPprob = reshape(Pprob,[size(Pprob,1)*size(Pprob,2) size(Pprob,3)]);
    rINprob = reshape(INprob,[size(INprob,1)*size(INprob,2) size(INprob,3)]);
    pbars = squeeze(nanmean(rPprob,1));  epbars = nanstd(rPprob,0,1)./sqrt(sum(~isnan(rPprob),1));
    inbars = squeeze(nanmean(rINprob,1)); einbars = nanstd(rINprob,0,1)./sqrt(sum(~isnan(rINprob),1));
    bars = cat(2,pbars',inbars'); ebars = cat(2,epbars',einbars');
    bar(alldep,bars);
    hold on;
    errorbar(alldep,bars(:,1),ebars(:,1),'LineStyle','none','Color','k'); hold on;
    errorbar(alldep,bars(:,2),ebars(:,2),'LineStyle','none','Color','k');
    set(gca,'TickDir','out','View',[90 90],'XDir','reverse','FontSize',14,'LineWidth',2,'Box','off');
    xlim([alldep(end)-.1 alldep(1)+.1]);
    xlabel('cortical depth (mm)'); ylabel('response probability');
    title(gca,fname,'interpreter','none');
    
    subplot(1,3,3)
    clear bars; 
    bars = cat(1,sum(~isnan(rPprob),1),sum(~isnan(rINprob),1)); 
    bar(alldep',bars'); 
    set(gca,'TickDir','out','View',[90 90],'XDir','reverse','FontSize',14,'LineWidth',2,'Box','off');
    xlim([alldep(end)-.1 alldep(1)+.1]);
    ylabel('n units');
    
    if saveon == 1
        saveas(gcf,strcat(spath,'avg_alldepths_sparseness'),'eps2c');
    end

else
    figure, set(gcf,'Color','w','Position',[1 1 700 700]);
    subplot(1,2,1)
    clear pbars inbars bars
    pbars = squeeze(nanmean(Pprob(1,:,:),2));
    inbars = squeeze(nanmean(INprob(1,:,:),2));
    bars = cat(2,pbars,inbars); ebars = cat(2,epbars',einbars');
    bar(alldep,bars);
    hold on;
    errorbar(alldep,bars(:,1),ebars(:,1),'LineStyle','none','Color','k'); hold on;
    errorbar(alldep,bars(:,2),ebars(:,2),'LineStyle','none','Color','k');
    set(gca,'TickDir','out','View',[90 90],'XDir','reverse','FontSize',14,'LineWidth',2,'Box','off');
    xlim([alldep(end)-.1 alldep(1)+.1]);
    xlabel('cortical depth (mm)'); ylabel('response probability');
    title(gca,strcat(fname,'electrode A'),'interpreter','none'); 
    
    subplot(1,2,2)
    clear pbars inbars bars
    pbars = squeeze(nanmean(Pprob(2,:,:),2));
    inbars = squeeze(nanmean(INprob(2,:,:),2));
    bars = cat(2,pbars,inbars); ebars = cat(2,epbars',einbars');
    bar(alldep,bars);
    hold on;
    errorbar(alldep,bars(:,1),ebars(:,1),'LineStyle','none','Color','k'); hold on;
    errorbar(alldep,bars(:,2),ebars(:,2),'LineStyle','none','Color','k');
    set(gca,'TickDir','out','View',[90 90],'XDir','reverse','FontSize',14,'LineWidth',2,'Box','off');
    xlim([alldep(end)-.1 alldep(1)+.1]);
    xlabel('cortical depth (mm)'); ylabel('response probability');
    title(gca,'electrode B'); 
         if saveon == 1
        saveas(gcf,strcat(spath,'alldepths_sparseness'),'eps2c');
    end   
end

% plot average response probability upper vs. deep with cell type as group
clear pbars inbars bars
Pupper  = squeeze(nanmean(rPprob(:,alldep>=0),1)); 
Pdeep = squeeze(nanmean(rPprob(:,alldep<0),1)); 
INupper  = squeeze(nanmean(rINprob(:,alldep>=0),1)); 
INdeep = squeeze(nanmean(rINprob(:,alldep<0),1)); 

pbars  = [nanmean(Pupper); nanmean(Pdeep)];  
epbars = [nanstd(Pupper,0,2)./sqrt(sum(~isnan(Pupper))); nanstd(Pdeep,0,2)./sqrt(sum(~isnan(Pdeep)))];

inbars = [squeeze(nanmean(INupper)); squeeze(nanmean(INdeep))]; 
einbars = [nanstd(INupper,0,2)./sqrt(sum(~isnan(INupper))); nanstd(INdeep,0,2)./sqrt(sum(~isnan(INdeep)))];
bars = cat(2,pbars,inbars); ebars = cat(2,epbars,einbars);

figure, set(gcf,'Color','w','Position',[1 1 500 300]);
subplot(1,2,1)
bar([1 2],bars);
hold on;
errorbar([1 2],bars(:,1),ebars(:,1),'LineStyle','none','Color','k'); hold on;
errorbar([1 2],bars(:,2),ebars(:,2),'LineStyle','none','Color','k');
set(gca,'TickDir','out','View',[90 90],'XDir','normal','FontSize',14,'LineWidth',2,'Box','off','XTickLabel',{'upper'; 'deep'});
xlim([.5 2.5]);
ylabel('response probability');title(gca,fname,'interpreter','none');

subplot(1,2,2)
clear bars
bars = cat(2,[sum(~isnan(Pupper)); sum(~isnan(Pdeep))],[sum(~isnan(INupper)); sum(~isnan(INdeep))]); 
bar(bars); 
set(gca,'TickDir','out','View',[90 90],'XDir','normal','FontSize',14,'LineWidth',2,'Box','off','XTickLabel',{'upper'; 'deep'});
xlim([.5 2.5]); ylabel('n units'); 
    if saveon == 1
        saveas(gcf,strcat(spath,'avg_uppervlower_sparseness'),'eps2c');
    end