% analyze drifting gratings files:
% set path
clear all
close all force
warning off

if ~ispc && ~exist('/users/kaciedougherty/documents/code')
    addpath('/volumes/drobo/users/kacie/code/processdata_code');
    addpath('/volumes/drobo/users/kacie/code/nbanalysis');
    addpath(genpath('/volumes/drobo/lab software/neurophys analysis'));
elseif ~ispc
    addpath(genpath('/users/kaciedougherty/documents/code/nbanalysis'))
    addpath('/users/kaciedougherty/documents/neurophysdata')
    addpath(genpath('/users/kaciedougherty/documents/code/BLACKROCK'))
else
    addpath('/users/MLab/documents');
    addpath('/users/MLab/documents/mlanalysisonline');
    addpath('/users/MLab/documents/utils');
end

%% set filename and parameters
unitid = '180405_I_cinterocdrft001';

% [unit] = getUnitInfo(unitid)
% theseelectrodes = {unit.channel};
% if ~isfield(unit,'flipeye'), unit.flipeye = 0, end
% if ~isfield(unit,'refresh'), unit.refresh = 85, end

theseelectrodes = {'eD18'};
unit.channel = 'eD18';
unit.refresh = 85;
unit.num     = 0;
unit.flipeye = 0;

PRE         = 0;            % ms relative to stim onset
POST        = 1000;         % ms relative to stim onset
BRdatafile  = unitid; %unit.cinteroc;

saveon      = 0;
chn         = 1;            %  channel number
un          = 1;            %  unit number
chooseunqc  = [];           %  subselect DE contrast levels // leave empty to use all levels in file
choosenunqc = [0 1];        %  subselect NDE contrast levels // leave empty to use all levels in file
colid       = [1 12 2 4];
use_evcodes = 1;
analycyc    = 1;
aligncodes  = 0;
UNITS       = unit;
el          = theseelectrodes{chn};
spiking     = 1;
csdonly     = 0;
crfs        = 1;
dofit       = 1;
binoccrf    = 1;
psth        = 0 ;
raster      = 1;

figtitle = strcat(BRdatafile(1:9),BRdatafile(end-2:end),'_ch',theseelectrodes{1}(3:4),'_u',num2str(un));

unq_unit = strcat(BRdatafile(1:8),'_',theseelectrodes{1},'_un',num2str(unit(1).num));

if strfind(BRdatafile,'pmk')
    paradigm = 'pmk';
else
    paradigm   = BRdatafile(10:end-3); paradigm = paradigm(1:4);
end

if ispc
    brdrname = strcat('\\129.59.230.179\CerebrusData\',BRdatafile(1:8));
    mldrname = strcat('\\129.59.230.116\MLData\',BRdatafile(1:8));
    
elseif exist('/users/kaciedougherty/documents/code')
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

if ~isempty(strfind(BRdatafile,'ori'))
    ext = '.gRFORIDRFTGrating_di';
elseif ~isempty(strfind(BRdatafile,'rfsf'))
    ext = '.gRFSFDRFTGrating_di';
elseif ~isempty(strfind(BRdatafile,'tfsf'))
    ext = '.gTFSFDRFTGrating_di';
elseif  ~isempty(strfind(BRdatafile,'size'))
    ext = '.gRFSIZEDRFTGrating_di';
elseif  ~isempty(strfind(BRdatafile,'cinteroc'))
    ext = '.gCINTEROCDRFTGrating_di';
elseif ~isempty(strfind(BRdatafile,'cpatch'))
    ext = '.gCPATCHDRFTGrating_di';
elseif ~isempty(strfind(BRdatafile,'pmk'))
    ext = '.gPMKDRFTGrating_di';
elseif ~isempty(strfind(BRdatafile,'cone'))
    ext = '.gCONEDRFTGrating_di';
elseif ~isempty(strfind(BRdatafile,'color'))
    ext = '.gCOLORFLICKERDRFTGrating_di';
end

badobs = [1];
flag_RewardedTrialsOnly = true;
grating = readgDRFTGrating([mldrname filesep BRdatafile ext]);

if strfind(BRdatafile,'cinteroc') & unit.flipeye == 1
    hcontrast = grating.contrast;
    grating.contrast = grating.fixedc;
    grating.fixedc = hcontrast;
    
end

if unit.refresh ~= 85
    h_temporalfreq =  grating.temporal_freq;
    grating.temporal_freq = [];
    grating.temporal_freq = round( h_temporalfreq.*(unit.refresh./85));
    
end

if strcmp(paradigm,'cone') && ~isfield(grating,'path') && strcmp(BRdatafile,'161214')
    % forgot to add path field initially in grating text file (12/14/16)
    currentdir = pwd;
    cd(brdrname);
    
    finfo = dir(strcat(BRdatafile,'*','GRATINGRECORD','*','.mat'));
    [~,sid] = sort([finfo(:).datenum],'ascend');
    hGRATINGRECORD = [];
    for fl = 1:length(finfo)
        
        load(finfo(sid(fl)).name);
        
        hGRATINGRECORD = [hGRATINGRECORD GRATINGRECORD];
        clear GRATINGRECORD;
    end
    GRATINGRECORD = hGRATINGRECORD; clear hGRATINGRECORD;
    grating.path = [];
    for tr = 1:max(grating.trial)
        grating.path = [grating.path GRATINGRECORD(tr).path GRATINGRECORD(tr).path];
    end
    
end
%% load digital codes and neural data:
filename = fullfile(brdrname,BRdatafile);
if strcmp(BRdatafile(1:6),'151214')
    filename(end-19:end-14) = '151314';
end
% check if file exist and load NEV
if exist(strcat(filename,'.nev'),'file') == 2;
    NEV = openNEV(strcat(filename,'.nev'),'read','overwrite');
    labels = {NEV.ElectrodesInfo.ElectrodeLabel};
    for i = 1:length(labels)
        if strfind(labels{i}',el)
            idxlab = i;
        end
    end
    pin = NEV.ElectrodesInfo(idxlab).ElectrodeID;
else
    error('the following file does not exist\n%s.nev',filename);
end
% get event codes from NEV
EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventTimes = floor(NEV.Data.SerialDigitalIO.TimeStampSec .* 1000); %ms, to match 1kHz
EventSampels = NEV.Data.SerialDigitalIO.TimeStamp;
[pEvC, pEvT] = parsEventCodesML(EventCodes,EventSampels);


kilo = strcat(brdrname,'/',BRdatafile(1:8),'_p01_kilo_ss.mat'); 
load(kilo,'-MAT')


% if ~ispc
%     
%     ex = exist(strcat('/volumes/drobo/data/neurophys/doughek/spikesorted_binoc/',BRdatafile,'_ss'));
%     if any(ex)
%         clear NEV; % load spike NEV
%         filename = strcat('/volumes/drobo/data/neurophys/doughek/spikesorted_binoc/',BRdatafile,'_ss');
%         if strcmp(BRdatafile(1:6),'151214')
%             filename = strcat('/volumes/drobo/data/neurophys/doughek/spikesorted_binoc/','151314_I_cinteroc001','_ss');
%             
%         end
%         NEV = openNEV(strcat(filename,'.nev'),'nomat','nosave');
%         elabel = NEV.ElectrodesInfo(pin).ElectrodeLabel;
%     else
%         elabel =  el;
%         clear NEV;
%        
%         if exist(strcat(brdrname,'/',BRdatafile,'_ss'))
%         filename = strcat(brdrname,'/',BRdatafile,'_ss');
%         else
%                 filename = strcat(brdrname,'/',BRdatafile);  
%         end
%         NEV = openNEV(strcat(filename,'.nev'),'nomat','nosave');
%         elabel = NEV.ElectrodesInfo(pin).ElectrodeLabel;
%     end
% else
%     
%     ex = exist(strcat('c:/users/mlab/documents/spikesorted_binoc/',BRdatafile,'_ss.nev'));
%     if any(ex)
%         clear NEV;
%         nevname = strcat('c:\users\mlab\documents\spikesorted_binoc\',BRdatafile,'_ss.nev');
%         NEV = openNEV(nevname,'nomat','nosave');
%         elabel = NEV.ElectrodesInfo(pin).ElectrodeLabel;
%     else
%         elabel =  el;
%     end
%     
% end


% sort/pick trials [before iterating unit]
if strfind(BRdatafile,'cone')
    stimfeatures = {...
        'tilt'...
        'sf'...
        'contrast'...
        'fixedc'...
        'diameter'...
        'eye'...
        'oridist'...
        'phase'...
        'temporal_freq'...
        'pathw'};
    
    grating.path([strcmp(grating.pathw,'LM')])  = 1;
    grating.path([strcmp(grating.pathw,'S')])   = 2;
    grating.path([strcmp(grating.pathw,'LMo')]) = 3;
    grating.path([strcmp(grating.pathw,'ach')]) = 4;
    grating.pathw = []; 
    grating.pathw = grating.path; grating.path = []; 
else
    stimfeatures = {...
        'tilt'...
        'sf'...
        'contrast'...
        'fixedc'...
        'diameter'...
        'eye'...
        'oridist'...
        'phase'...
        'temporal_freq'};
end


clear(stimfeatures{:})

extension = 'ns6';

[STIM,spkTPs,all_cyc,realtr] = sortTrialData(BRdatafile,brdrname,pEvC,pEvT,grating,stimfeatures,flag_RewardedTrialsOnly,badobs,extension,use_evcodes);

if ~isempty(all_cyc)
    photo_yes = 1;
else
    photo_yes = 0;
end

%% which electrode contact?
% get pin id
if spiking == 1
    if ~isempty(el)
        if size(elabel,1) > 1
            elabel = elabel';
        end
        eidx = find(cell2mat(cellfun(@(x) ~isempty(strfind(x',elabel)),{NEV.ElectrodesInfo.ElectrodeLabel},'UniformOutput',0)));
        if isempty(eidx)
            error('no %s',elabel)
        end
    else
        % run one channel with sorted data:
        eidx = unique(NEV.Data.Spikes.Electrode);
        eidx = eidx(1);
    end
    
    % get spikes for unit of interest:
    eI    =  NEV.Data.Spikes.Electrode == eidx;
    unit  = UNITS(chn).num(un);
    if unit > 0
        I = eI &  NEV.Data.Spikes.Unit == unit;
    else
        I = eI;
    end
    SPK   = double(NEV.Data.Spikes.TimeStamp(I)); % in samples
end

wvform = double(NEV.Data.Spikes.Waveform(:,I))./4; 
figure,set(gcf,'Color','w','Position',[1 1 500 300]); 
plot(mean(wvform,2),'Color','k','LineWidth',3);
hold on; 
errbar = std(wvform,0,2)./sqrt(size(wvform,2)); 
  shadedErrorBar([1:size(wvform,1)],mean(wvform,2),errbar,{'-','color',[0 0 0]});
set(gca,'Box','off','TickDir','out','FontSize',16); 
xlim([1 size(wvform,1)]); 

%% trigger data
fs     = double(NEV.MetaTags.SampleRes);
if photo_yes == 1
WPRE   = PRE - ((1000/grating.temporal_freq(1))./2); % OR +
WPOST  = POST - ((1000/grating.temporal_freq(1))./2);
else
    WPRE = PRE; 
    WPOST = POST; 
end

pre    = WPRE.*(fs/1000);
post   = WPOST.*(fs/1000);

svec   = floor((PRE.*(fs/1000))): floor((POST.*(fs/1000)));
tvec   = svec./(fs/1000);

if csdonly == 1
    [r_spkTPs,r_all_cyc,r_STIM] = removeAnaNaNTrials(spkTPs,all_cyc,STIM);
    tf      = nanunique(STIM.temporal_freq);
    per    =  floor(1000./tf).*(fs./1000);
else
    [spktr,spkcyc] = triggerDRFTdata(SPK,STIM,spkTPs,all_cyc,fs,pre,post,photo_yes);
    [r_spktr,r_spkcyc,r_spktintr,r_trls,r_spktincyc,r_cycn,r_cyctn,r_STIM] = removeNaNTrials(spktr,spkcyc,STIM);
end

[sdf, sua, tm] = spk2sdf(SPK,fs);
[STIM,realtr]  = sortStimandTimeData(grating,pEvC,pEvT);
triggerpts     = double(STIM.onsets);
sdftr          = squeeze(trigData(sdf',triggerpts./30,-PRE,POST));
trls           = STIM.contrast == 1 & STIM.fixedc == 1; 
resp           = sdftr(:,trls);
% 
% nfft = 512;
% win  = 512;
% noverlap = win./2;
% for tr = 1:size(resp,2)
%     [Pxx, F] = new_psd(resp(:,tr),nfft, fs, win, noverlap);
%     power(:,tr) = Pxx;
% end
% 
% figure,
% plot(F(f<30),mean(power(f<30,:),2));

    for tr = 1:size(resp,2)
        
        spkdens         = resp(:,tr);
        NFFT            = 2^nextpow2(length(spkdens));       % Next power of 2 from length of y
        Y               = fft(spkdens,NFFT)/length(spkdens); % keep units
        pvec(:,tr)      = abs(Y(1:NFFT/2+1,:));
    end
        
        f               = (1000./1)/2*linspace(0,1,NFFT/2+1);
  figure, plot(f(f<30),mean(pvec(f<30,:),2))      
                                                     

%% plot data

if strcmp(paradigm,'colo') || strcmp(paradigm,'cone')
    figsize = [500 500 500 800];
else
    figsize = [500 500 1000 800];
end

global figH

if csdonly ~= 1
    
    switch paradigm
        
        case 'cint'
            
            figure,set(gcf,'Color','w'),
            text(0.1,0.7,BRdatafile,'interpreter','none','FontSize',20,'FontName','arial');
            hold on;
            text(0.1,0.4,unq_unit,'interpreter','none','FontSize',20,'FontName','arial');
            hold on;
            text(0.1,0.1,strcat(' run on ',date),'interpreter','none','FontSize',20,'FontName','arial');
            set(gca,'XColor','w','YColor','w')
            
            
            unqc  = unique(STIM.contrast(STIM.contrast>=0));
            if ~isempty(chooseunqc)
            unqc  = unqc(ismember(unqc,chooseunqc));
            end
            
            nunqc = unique(STIM.fixedc(STIM.fixedc>=0));
            nunqc = nunqc(ismember(nunqc,choosenunqc));
            
            tf    = nanunique(STIM.temporal_freq);
            ct = 0; psdct = 0; 
            for c  = 1:length(unqc)
                for fc = 1:length(nunqc)
                    
                    thesetrs = find(r_STIM.contrast == unqc(c) & r_STIM.fixedc == nunqc(fc));
                           
                    % exclude trials with 1 or 0 spikes: 
                    thesetrs = thesetrs(sum(r_spktr(:,thesetrs),1) > 1); 
                    h_thesetrs = zeros(length(r_STIM.contrast),1); 
                    h_thesetrs(thesetrs) = 1; 
                    clear thesetrs; thesetrs = h_thesetrs; 
                  
                    ct = ct + 1;
                    % raster:
                    fh=findall(0,'type','figure');
                    if isempty(fh),
                        figH.N = 1;
                        figH.subplot = ct;
                        figH.row    = length(unqc);
                        figH.col    = length(nunqc);
                        rasterN = figH.N;
                    else
                        if exist('rasterN')
                            figH.N = rasterN;
                        else
                            figH.N = fh(1).Number + 1;
                            rasterN = figH.N;
                        end
                        figH.subplot = ct;
                        figH.row    = length(unqc);
                        figH.col    = length(nunqc);
                        
                    end
                    shufthesetrs = find(thesetrs); 
                    shufthesetrs = shufthesetrs(randperm(length(shufthesetrs))); 
                    plotSpikeRaster(logical(r_spktr(:,shufthesetrs)'),'PlotType','scatter');
                    setFigure(gcf,gca,figsize,1,tvec,PRE,STIM)
                    title(gca,strcat(' DE ',num2str(unqc(c)),' NDE ',num2str(nunqc(fc))),'FontSize',12)
                    xlim([1 length([pre:post])]);
                    
                    
                    % psth:
                    [tedges,resp,binsize] = myPSTH(length([pre:post]),PRE,thesetrs,r_spktintr,r_trls,fs);
                    fh=findall(0,'type','figure');
                    if isempty(fh),
                        figP = 1;
                    else
                        figP = fh(1).Number + 1;
                    end
                    if ct == 1,
                        figure(figP), psthN = figP;
                    else
                        figure(psthN),
                    end
                    subplot(length(unqc),length(nunqc),ct)
                    ph = bar(tedges(1:end-1),mean(resp(1:end-1,:),2),'histc');
                    setFigure(gcf,gca,figsize,2,tvec,PRE,STIM,tedges);
                    title(gca,strcat(' DE ',num2str(unqc(c)),' NDE ',num2str(nunqc(fc))),'FontSize',12)
                    xlim([PRE POST]); set(ph,'FaceColor','k');
                    ylim([0 90]); hold on;
              
                 if size(resp,2) > 1
                     psdct = psdct + 1; 
                   
                    [zF1, f,freqrange,frid,pvec,meanpower,stdpower,powF1] = myzF1(resp,binsize,tf);
                    
                    fh=findall(0,'type','figure');
                    if isempty(fh),
                        figPS = 1;
                    else
                        figPS = fh(1).Number + 1;
                    end
                    if psdct == 1,
                        figure(figPS), psdN = figPS;
                    else
                        figure(psdN),
                    end
                    subplot(length(unqc),length(nunqc),ct),
                    plot(f(f>0 & f<30),mean(pvec(f>0 & f<30,:),2),'-o','LineWidth',1.5); h = vline(tf); set(h,'LineWidth',2);
                    set(gca,'Box','off','TickDir','out','FontSize',16);
                    title(gca,strcat(' DE ',num2str(unqc(c)),' NDE ',num2str(nunqc(fc))),'FontSize',12);
                    setFigure(gcf,gca,figsize,0); if unqc(c) == unqc(end), xlabel('freq (Hz)'); end
                    if unqc(c) == unqc(floor(length(unqc)./2)),
                        ylabel('spks/s'),
                    end
                    
                    
                    zF1s(c,fc)      = nanmean(zF1,2);
                    stdzF1s(c,fc)   = nanstd(zF1,0,2)./(sqrt(size(zF1,2)));
                    trzF1s{c,fc}    = zF1;
                    
                    spkrTF(c,fc)    = nanmean(pvec(frid,:),2);
                    stdspkrTF(c,fc) = nanstd(pvec(frid,:),0,2)./sqrt(size(pvec,2));
                    trspkrTF{c,fc}  = pvec(frid,:);
                    
                    F1s(c,fc)       = nanmean(pvec(frid,:)./(pvec(1,:)),2);
                    stdF1s(c,fc)    = nanstd(pvec(frid,:)./(pvec(1,:)),0,2)./(sqrt(size(pvec,2)));
                    trF1s{c,fc}     = pvec(frid,:)./(pvec(1,:));
                    
                    F0s(c,fc)       = nanmean(mean(resp,1),2);
                    stdF0s(c,fc)    = nanstd(mean(resp,1),0,2)./(sqrt(size(resp,2)));
                    trF0s{c,fc}     = nanmean(resp,1);
                 else
                        
                    zF1s(c,fc)      = nan;
                    stdzF1s(c,fc)   = nan;
                    trzF1s{c,fc}    = nan;
                    
                    spkrTF(c,fc)    = nan;
                    stdspkrTF(c,fc) = nan;
                    trspkrTF{c,fc}  = nan;
                    
                    F1s(c,fc)       = nan;
                    stdF1s(c,fc)    = nan;
                    trF1s{c,fc}     = nan;
                    
                    F0s(c,fc)       = nan;
                    stdF0s(c,fc)    = nan;
                    trF0s{c,fc}     = nan;
                     
                 end
                end
            end
            
            figure, set(gcf,'Color','w','Position',[1 1 800 350],'PaperPositionMode','auto');
            subplot(1,4,1)
            for fc = 1:length(nunqc)
                h = plot(unqc,zF1s(:,fc),'-o','LineWidth',2); rgbcol = get(h,'Color');
                hold on;
                errorbar(unqc,zF1s(:,fc),stdzF1s(:,fc),'LineStyle','none','Color',rgbcol);
                hold on;
            end
            xlim([0 1]); set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('zF1'); xlabel('contrast DE (%)');
            
            subplot(1,4,2)
            for fc = 1:length(nunqc)
                h = plot(unqc,spkrTF(:,fc),'-o','LineWidth',2); rgbcol = get(h,'Color');
                hold on;
                errorbar(unqc,spkrTF(:,fc),stdspkrTF(:,fc),'LineStyle','none','Color',rgbcol);
                hold on;
            end
            xlim([0 1]); set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('spks/s at TF'); xlabel('contrast DE (%)');
            
            subplot(1,4,3)
            for fc = 1:length(nunqc)
                h = plot(unqc,F1s(:,fc),'-o','LineWidth',2); rgbcol = get(h,'Color');
                hold on;
                errorbar(unqc,F1s(:,fc),stdF1s(:,fc),'LineStyle','none','Color',rgbcol);
                hold on;
            end
            xlim([0 1]); set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('F1'); xlabel('contrast DE (%)');
            
            subplot(1,4,4)
            for fc = 1:length(nunqc)
                h = plot(unqc,F0s(:,fc),'-o','LineWidth',2); rgbcol = get(h,'Color');
                hold on;
                errorbar(unqc,F0s(:,fc),stdF0s(:,fc),'LineStyle','none','Color',rgbcol);
                hold on;
            end
            xlim([0 1]); set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('F0 (spk/s)'); xlabel('contrast DE (%)');
            
            if dofit == 1
                figure, set(gcf,'Color','w','Position',[1 1 800 350],'PaperPositionMode','auto');
                subplot(1,4,1)
                for fc = 1:length(nunqc)
                    h = plot(unqc.*100,zF1s(:,fc),'o','LineWidth',2); rgbcol = get(h,'Color'); set(h,'MarkerFaceColor',rgbcol);
                    hold on;
                    errorbar(unqc.*100,zF1s(:,fc),stdzF1s(:,fc),'LineStyle','none','Color',rgbcol);
                    hold on;
                    [prd,xprd] = runCRFFit(zF1s(:,fc),unqc);
                    hold on;
                    plot(xprd,prd,'-','LineWidth',2,'Color',rgbcol);
                    clear h prd xprd
                end
                xlim([0 100]); set(gca,'Box','off','TickDir','out','FontSize',16);
                ylabel('zF1'); xlabel('contrast DE (%)');
              
                subplot(1,4,2)
                for fc = 1:length(nunqc)
                    h = plot(unqc.*100,spkrTF(:,fc),'o','LineWidth',2); rgbcol = get(h,'Color'); set(h,'MarkerFaceColor',rgbcol);
                    hold on;
                    errorbar(unqc.*100,spkrTF(:,fc),stdspkrTF(:,fc),'LineStyle','none','Color',rgbcol);
                    hold on;
                    [prd,xprd] = runCRFFit(spkrTF(:,fc),unqc);
                    hold on;
                    plot(xprd,prd,'-','LineWidth',2,'Color',rgbcol);
                    clear h prd xprd
                end
                xlim([0 100]); set(gca,'Box','off','TickDir','out','FontSize',16);
                ylabel('spk/s'); xlabel('contrast DE (%)');
                
                subplot(1,4,3)
                for fc = 1:length(nunqc)
                    h = plot(unqc.*100,F1s(:,fc),'o','LineWidth',2); rgbcol = get(h,'Color'); set(h,'MarkerFaceColor',rgbcol);
                    hold on;
                    errorbar(unqc.*100,F1s(:,fc),stdF1s(:,fc),'LineStyle','none','Color',rgbcol);
                    hold on;
                    [prd,xprd] = runCRFFit(F1s(:,fc),unqc);
                    hold on;
                    plot(xprd,prd,'-','LineWidth',2,'Color',rgbcol);
                    clear h prd xprd
                end
                xlim([0 100]); set(gca,'Box','off','TickDir','out','FontSize',16);
                ylabel('F1'); xlabel('contrast DE (%)');
                
                
                subplot(1,4,4)
                for fc = 1:length(nunqc)
                    h = plot(unqc.*100,F0s(:,fc),'o','LineWidth',2); rgbcol = get(h,'Color'); set(h,'MarkerFaceColor',rgbcol);
                    hold on;
                    errorbar(unqc.*100,F0s(:,fc),stdF0s(:,fc),'LineStyle','none','Color',rgbcol);
                    hold on;
                    [prd,xprd] = runCRFFit(F0s(:,fc),unqc);
                    hold on;
                    plot(xprd,prd,'-','LineWidth',2,'Color',rgbcol);
                    clear h prd xprd
                end
                xlim([0 100]); set(gca,'Box','off','TickDir','out','FontSize',16);
                ylabel('F0'); xlabel('contrast DE (%)');
            end
            
            figure, set(gcf,'Color','w','Position',[1 1 250 300],'PaperPositionMode','auto');
            for fc = 1:length(nunqc)
                h = plot(unqc.*100,F0s(:,fc),'o','LineWidth',2); rgbcol = get(h,'Color'); set(h,'MarkerFaceColor',rgbcol);
                hold on;
                errorbar(unqc.*100,F0s(:,fc),stdF0s(:,fc),'LineStyle','none','Color',rgbcol);
                hold on;
                [prd,xprd] = runCRFFit(F0s(:,fc),unqc);
                hold on;
                plot(xprd,prd,'-','LineWidth',2,'Color',rgbcol);
                clear h prd xprd
            end
            xlim([0 100]); set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('F0'); xlabel('contrast DE (%)');
            saveas(gcf,strcat(savedrname,unq_unit,'_crfF0'),'eps2c'); 
            
            figure, set(gcf,'Color','w','Position',[1 1 250 300],'PaperPositionMode','auto');
            
                for fc = 1:length(nunqc)
                    h = plot(unqc.*100,zF1s(:,fc),'o','LineWidth',2); rgbcol = get(h,'Color'); set(h,'MarkerFaceColor',rgbcol);
                    hold on;
                    errorbar(unqc.*100,zF1s(:,fc),stdzF1s(:,fc),'LineStyle','none','Color',rgbcol);
                    hold on;
                    [prd,xprd] = runCRFFit(zF1s(:,fc),unqc);
                    hold on;
                    plot(xprd,prd,'-','LineWidth',2,'Color',rgbcol);
                    clear h prd xprd
                end
                xlim([0 100]); set(gca,'Box','off','TickDir','out','FontSize',16);
                ylabel('zF1'); xlabel('contrast DE (%)');
               saveas(gcf,strcat(savedrname,unq_unit,'_crfzF1'),'eps2c'); 
         jkjk
            
            if saveon == 1
                fh=findall(0,'type','figure','Name','');
                
                for f  = length(fh):-1:1
                    get(figure(fh(f).Number))
                    if isempty(fh(1).Name)
                        export_fig(strcat(savedrname,unq_unit),'-pdf','-nocrop','-append')
                    end
                end
            end
            
            
            [anovapvals,ttestpvals] = runCRFstats(unqc,nunqc,trzF1s);
            fprintf('\n anova zF1 p = %f\n',anovapvals(1));
            if any(ttestpvals<0.05)
                fprintf(' ttest zF1 contrast %f p = %f\n',unqc(find(ttestpvals<.05)),ttestpvals(find(ttestpvals<.05)));
            end
            clear anovapvals ttestpvals
            
            [anovapvals,ttestpvals] = runCRFstats(unqc,nunqc,trspkrTF);
            fprintf('\n anova s @ TF p = %f\n',anovapvals(1));
            if any(ttestpvals<0.05)
                fprintf(' ttest s @ TF contrast %f p = %f\n',unqc(find(ttestpvals<.05)),ttestpvals(find(ttestpvals<.05)));
            end
            clear anovapvals ttestpvals
            
            [anovapvals,ttestpvals] = runCRFstats(unqc,nunqc,trF1s);
            fprintf('\n anova F1 p = %f\n',anovapvals(1));
            if any(ttestpvals<0.05)
                fprintf(' ttest F1 contrast %f p = %f\n',unqc(find(ttestpvals<.05)),ttestpvals(find(ttestpvals<.05)));
            end
            clear anovapvals ttestpvals
            
            
            [anovapvals,ttestpvals] = runCRFstats(unqc,nunqc,trF0s);
            fprintf('\n anova F0 p = %f\n',anovapvals(1));
            if any(ttestpvals<0.05)
                fprintf(' ttest F0 contrast %f p = %f\n',unqc(find(ttestpvals<.05)),ttestpvals(find(ttestpvals<.05)));
            end
            clear anovapvals ttestpvals
            
        case 'colo'
            
            figure,set(gcf,'Color','w'),
            text(0.1,0.7,BRdatafile,'interpreter','none','FontSize',20,'FontName','arial');
            hold on;
            text(0.1,0.4,unq_unit,'interpreter','none','FontSize',20,'FontName','arial');
            hold on;
            text(0.1,0.1,strcat(' run on ',date),'interpreter','none','FontSize',20,'FontName','arial');
            set(gca,'XColor','w','YColor','w')
            
            colors = nanunique(STIM.contrast);
            tf     = nanunique(STIM.temporal_freq);
            ct = 0;
            for c = 1:length(colors)
                ct = ct + 1;
                thesetrs = r_STIM.contrast == colors(c);
                
                % raster:
                fh=findall(0,'type','figure');
                if isempty(fh),
                    figH.N    = 1;
                    figH.row  = length(colors);
                    figH.col  = 1;
                    figH.subplot = ct;
                    rasterN = figH.N;
                else
                    if exist('rasterN')
                        figH.N = rasterN;
                    else
                        figH.N = fh(1).Number + 1;
                        
                    end
                    figH.row  = length(colors);
                    figH.col  = 1;
                    figH.subplot = ct;
                    rasterN = figH.N;
                end
                
                plotSpikeRaster(logical(r_spktr(:,thesetrs)'),'PlotType','scatter');
                setFigure(gcf,gca,figsize,1,tvec,PRE,STIM)
                title(gca,strcat(' color ',num2str(colors(c))),'FontSize',12)
                xlim([1 length([pre:post])]);
                
     
                
                % psth:
                fh=findall(0,'type','figure');
                if isempty(fh),
                    psthN = 1 ;
                else
                    
                    psthN = fh(1).Number + 1;
                end
                
                [tedges,resp] = myPSTH(length([pre:post]),PRE,thesetrs,r_spktintr,r_trls,fs);
                
                figure(psthN),
                subplot(length(colors),1,ct)
                ph = bar(tedges(1:end-1),mean(resp(1:end-1,:),2),'histc');
                setFigure(gcf,gca,figsize,2,tvec,PRE,STIM,tedges);
                title(gca,strcat('color',num2str(colors(c))));
                xlim([PRE POST]); set(ph,'FaceColor','k'); ylim([0 70]); 
                
                % on or off?
                if colors(c) == -4
                    
                    md    = ceil(size(resp,1)./2);
                    
                    means = [mean(mean(resp(1:md,:),2),1); mean(mean(resp(md+1:end,:),2),1)];
                    
                    if means(1) > means(2)
                        
                        cell = 'off';
                        m_bsl   = mean(mean(resp(md+1:end,:),1),2);
                        std_bsl = nanstd(mean(resp(md+1:end,:),1),0,2);
                        thr     = (3*m_bsl);
                        faker   = mean(resp(1:md,:));
                        
                    else
                        
                        cell = 'on';
                        m_bsl   = mean(mean(resp(1:md,:),1),2);
                        std_bsl = nanstd(mean(resp(1:md,:),1),0,2);
                        thr     = (5*m_bsl);
                        faker   = mean(resp(md+1:end,:),2);
                        onlat   = tedges(find(faker>thr,1,'first')); % Jiang et al. (2015): P (25.2, 2.0) M  (19.6, 2.2)
                        peaklat    = tedges(find(faker == max(faker)));
                    end
                    
                end
                
                trresp{c} = resp;
                clear thesetrs;
            end
            
            figure,set(gcf,'Color','w','Position',[1 1 1000 300],'PaperPositionMode','auto');
            for c = 1:length(colors)
                subplot(1,length(colors),c)
                md    = floor(size(trresp{c},1)./2);
                means = [mean(mean(trresp{c}(1:md-1,:),2),1); mean(mean(trresp{c}(md:end,:),2),1)];
                vars  = [std(mean(trresp{c}(1:md-1,:),1),0,2); std(mean(trresp{c}(md:end,:),1),0,2)];
                [~,p] = ttest2(mean(trresp{c}(1:md-1,:),1),mean(trresp{c}(md:end,:),1));
                h = bar([1 2],means); set(h,'FaceColor',getColor(5));
                hold on;
                errorbar([1 2],means,vars,'LineStyle','none','Color', getColor(5));
                set(gca,'Box','off','TickDir','out','FontSize',16)
                title(gca,sprintf('p = %0.4f',p));
                if colors(c) == -1
                    set(gca,'XTickLabel',{'-LM';'+LM'});
                elseif colors(c) == -2
                    set(gca,'XTickLabel',{'S-';'S+'});
                elseif colors(c) == -3
                    set(gca,'XTickLabel',{'+L-M';'+M-L'});
                else
                    set(gca,'XTickLabel',{'black';'white'});
                end
                if c == 1
                    ylabel('spks/sec');
                end
                xlim([0 3]);   
            end
            
            if saveon == 1
                fh=findall(0,'type','figure');
                
                for f  = length(fh):-1:1
                    get(figure(fh(f).Number))
                    export_fig(strcat(savedrname,unq_unit),'-pdf','-nocrop','-append')
                    
                end
            end
            
        case 'rfor'
            
            figure,set(gcf,'Color','w'),
            text(0.1,0.7,BRdatafile,'interpreter','none','FontSize',20,'FontName','arial');
            hold on;
            text(0.1,0.4,unq_unit,'interpreter','none','FontSize',20,'FontName','arial');
            hold on;
            text(0.1,0.1,strcat(' run on ',date),'interpreter','none','FontSize',20,'FontName','arial');
            set(gca,'XColor','w','YColor','w')
            
            
            ori = nanunique(r_STIM.tilt);
            tf  = nanunique(r_STIM.temporal_freq);
            
            for o = 1:length(ori)
                
                thesetrs = r_STIM.tilt == ori(o);
                
                %raster:
                fh=findall(0,'type','figure');
                if isempty(fh),
                    figH.N       = 1;
                    figH.subplot = o;
                    figH.row     = floor(length(ori)./2);
                    figH.col     = 2;
                    rasterN      = figH.N;
                else
                    if exist('rasterN')
                        figH.N   = rasterN;
                    else
                        figH.N   = fh(1).Number + 1;
                        rasterN  = figH.N;
                    end
                    figH.subplot = o;
                    figH.row     = floor(length(ori)./2);
                    figH.col     = 2;
                    
                end
                
                plotSpikeRaster(logical(r_spktr(:,thesetrs)'),'PlotType','scatter');
                setFigure(gcf,gca,figsize,1,tvec,PRE,STIM)
                title(gca,strcat(' ori ',num2str(ori(o))),'FontSize',12)
                xlim([1 length([pre:post])]);
                
                % psth:
                [tedges,resp,binsize] = myPSTH(length([pre:post]),PRE,thesetrs,r_spktintr,r_trls,fs);
                fh=findall(0,'type','figure');
                if isempty(fh),
                    figP = 1;
                else
                    figP = fh(1).Number + 1;
                end
                if o == 1,
                    figure(figP), psthN = figP;
                else
                    figure(psthN),
                end
                subplot(floor(length(ori)./2),2,o)
                ph = bar(tedges(1:end-1),mean(resp(1:end-1,:),2),'histc');
                setFigure(gcf,gca,figsize,2,tvec,PRE,STIM,tedges);
                title(gca,strcat(' ori ',num2str(ori(o))),'FontSize',12)
                xlim([PRE POST]); set(ph,'FaceColor','k')
                ylim([0 90]); hold on;
                
                [zF1, f,freqrange,frid,pvec,meanpower,stdpower,powF1] =myzF1(resp,binsize,tf);
                
                fh=findall(0,'type','figure');
                if isempty(fh),
                    figPS = 1;
                else
                    figPS = fh(1).Number + 1;
                end
                if o == 1,
                    figure(figPS), psdN = figPS;
                else
                    figure(psdN),
                end
                subplot(floor(length(ori)./2),2,o), plot(f(f>0 & f<30),pvec(f>0 & f<30),'-o'); vline(tf);
                set(gca,'Box','off','TickDir','out','FontSize',16);
                title(gca,strcat('ORI',num2str(ori(o))));
                setFigure(gcf,gca,figsize,0)
                
                zF1s(o)      = mean(zF1,2);
                stdzF1s(o)   = std(zF1,0,2)./(sqrt(size(zF1,2)));
                trzF1s{o}    = zF1;
                
                spkrTF(o)    = mean(pvec(frid,:),2);
                stdspkrTF(o) = std(pvec(frid,:),0,2)./sqrt(size(pvec,2));
                trspkrTF{o}  = pvec(frid,:);
                
                F1s(o)       = mean(pvec(frid,:)./(pvec(1,:)),2);
                stdF1s(o)    = std(pvec(frid,:)./(pvec(1,:)),0,2)./(sqrt(size(pvec,2)));
                trF1s{o}     = pvec(frid,:)./(pvec(1,:));
                
                F0s(o)       = mean(mean(resp,1),2);
                stdF0s(o)    = std(mean(resp,1),0,2)./(sqrt(size(resp,1)));
                trF0s{o}     = mean(resp,1);
                
                realORI(o)  = ori(o);
                clear thesetrs;
            end
            
            % Ringach et al. (2002)
            
            figure,set(gcf,'Color','w','Position',[1 1 1000 400],'PaperPositionMode','auto'),
            
            subplot(1,3,1),
            h = plot(realORI,F1s,'-o','LineWidth',2); rgbcol = get(h,'Color'); hold on;
            errorbar(realORI,F1s,stdF1s,'LineStyle','none','Color',rgbcol);
            clear CV,CV = 1 - abs(sum(F1s.*exp(2i*deg2rad(ori)))./(sum(F1s))); % 0 perfectly tuned, 1 not at all tuned
            [~,mx] = max(F1s); [~,mn] = min(F1s); vline(ori([mn mx]));
            set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('F1');title(gca,sprintf('max: %u, min: %u, CV: %f',floor(ori(mn)),floor(ori(mx)),CV),'fontweight','normal','fontsize',14);
            xlim([ori(1) ori(end)]);
            fprintf('\nF1: max: %u, min: %u, CV: %f\n',floor(ori(mn)),floor(ori(mx)),CV);
            
            subplot(1,3,2),
            plot(realORI,F0s,'-o','LineWidth',2); hold on;
            errorbar(realORI,F0s,stdF0s,'LineStyle','none','Color',rgbcol);
            clear CV,CV = 1 - abs(sum(F0s.*exp(2i*deg2rad(ori)))./(sum(F0s))); % 0 perfectly tuned, 1 not at all tuned
            [~,mx] = max(F0s); [~,mn] = min(F0s); vline(ori([mn mx]));
            set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('F0');xlabel('ori (deg)');title(gca,sprintf('max: %u, min: %u, CV: %f',floor(ori(mn)),floor(ori(mx)),CV),'fontweight','normal','fontsize',14);
            xlim([ori(1) ori(end)]);
            fprintf('\nF0: max: %u, min: %u, CV: %f\n',floor(ori(mn)),floor(ori(mx)),CV);
            
            subplot(1,3,3),plot(realORI,zF1s,'-o','LineWidth',2); hold on;
            errorbar(realORI,zF1s,stdzF1s,'LineStyle','none','Color',rgbcol);
            clear CV,CV = 1 - abs(sum(zF1s.*exp(2i*deg2rad(ori)))./(sum(zF1s))); % 0 perfectly tuned, 1 not at all tuned
            [~,mx] = max(zF1s); [~,mn] = min(zF1s); vline(ori([mn mx]));
            set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('zF1'); title(gca,sprintf('max: %u, min: %u, CV: %f',floor(ori(mn)),floor(ori(mx)),CV),'fontweight','normal','fontsize',14);
            xlim([ori(1) ori(end)]);
            fprintf('\nzF1: max: %u, min: %u, CV: %f\n',floor(ori(mn)),floor(ori(mx)),CV);
            
            if dofit == 1
                
                figure, set(gcf,'Color','w','Position',[500 500 800 800],'PaperPositionMode','auto')
                subplot(3,2,1)
                polar([deg2rad(ori) (deg2rad(ori)+pi) deg2rad(ori(1))],[zF1s zF1s zF1s(1)]);
                
                subplot(3,2,2)
                y = [zF1s zF1s]';
                x = [ori (ori+180)]';
                f = fit(x,y,'smoothingspline');
                plot(f,x,y); hold on
                h = plot(x,y,'o'); rgbcol = get(h,'Color'); set(h,'MarkerFaceColor',rgbcol); hold on
                errorbar(x,y,[stdzF1s stdzF1s],'LineStyle','none','Color',rgbcol);
                axis tight; axis square; box off; legend('off')
                plot([180 180],ylim,'k:');
                set(gca,'TickDir','out','FontSize',16); xlabel('ori (deg)'); ylabel('zF1');
                
                subplot(3,2,3)
                polar([deg2rad(ori) (deg2rad(ori)+(pi)) deg2rad(ori(1))],[F1s F1s F1s(1)]);
                
                subplot(3,2,4)
                y = [F1s F1s]';
                x = [ori (ori+180)]';
                f = fit(x,y,'smoothingspline');
                plot(f,x,y); hold on
                h = plot(x,y,'o'); rgbcol = get(h,'Color'); set(h,'MarkerFaceColor',rgbcol); hold on
                errorbar(x,y,[stdF1s stdF1s],'LineStyle','none','Color',rgbcol);
                axis tight; axis square; box off; legend('off')
                plot([180 180],ylim,'k:');
                set(gca,'TickDir','out','FontSize',16); xlabel('ori (deg)'); ylabel('F1');
                
                subplot(3,2,5)
                polar([deg2rad(ori) (deg2rad(ori)+pi) deg2rad(ori(1))],[F0s F0s F0s(1)]);
                
                subplot(3,2,6)
                y = [F0s F0s]';
                x = [ori (ori+180)]';
                f = fit(x,y,'smoothingspline');
                plot(f,x,y); hold on
                h = plot(x,y,'o'); rgbcol = get(h,'Color'); set(h,'MarkerFaceColor',rgbcol); hold on
                errorbar(x,y,[stdF0s stdF0s],'LineStyle','none','Color',rgbcol);
                axis tight; axis square; box off; legend('off')
                plot([180 180],ylim,'k:');
                set(gca,'TickDir','out','FontSize',16); xlabel('ori (deg)'); ylabel('F0');
                
            end
            
            if saveon == 1
                fh=findall(0,'type','figure');
                
                for f  = length(fh):-1:1
                    get(figure(fh(f).Number))
                    export_fig(strcat(savedrname,unq_unit),'-pdf','-nocrop','-append')
                    
                end
            end
            
        case 'rfsf'
            
            figure,set(gcf,'Color','w'),
            text(0.1,0.7,BRdatafile,'interpreter','none','FontSize',20,'FontName','arial');
            hold on;
            text(0.1,0.4,unq_unit,'interpreter','none','FontSize',20,'FontName','arial');
            hold on;
            text(0.1,0.1,strcat(' run on ',date),'interpreter','none','FontSize',20,'FontName','arial');
            set(gca,'XColor','w','YColor','w')
            
            
            sf  = 1./nanunique(r_STIM.sf);
            tf  = nanunique(r_STIM.temporal_freq);
            
            for s = 1:length(sf)
                
                thesetrs = r_STIM.sf == 1./sf(s);
                
                %raster:
                fh=findall(0,'type','figure');
                if isempty(fh),
                    figH.N = 1;
                    figH.subplot = s;
                    figH.row    = ceil(length(sf)./2);
                    figH.col    = 2;
                    rasterN = figH.N;
                else
                    if exist('rasterN')
                        figH.N = rasterN;
                    else
                        figH.N = fh(1).Number + 1;
                        rasterN = figH.N;
                    end
                    figH.subplot = s;
                    figH.row    = ceil(length(sf)./2);
                    figH.col    = 2;
                    
                end
                
                plotSpikeRaster(logical(r_spktr(:,thesetrs)'),'PlotType','scatter');
                setFigure(gcf,gca,figsize,1,tvec,PRE,STIM)
                title(gca,strcat(' sf ',num2str(sf(s))),'FontSize',12)
                xlim([1 length([pre:post])]);
                
                % psth:
                [tedges,resp,binsize] = myPSTH(length([pre:post]),PRE,thesetrs,r_spktintr,r_trls,fs);
                fh=findall(0,'type','figure');
                if isempty(fh),
                    figP = 1;
                else
                    figP = fh(1).Number + 1;
                end
                if s == 1,
                    figure(figP), psthN = figP;
                else
                    figure(psthN),
                end
                
                subplot(ceil(length(sf)./2),2,s)
                ph = bar(tedges(1:end-1),mean(resp(1:end-1,:),2),'histc');
                setFigure(gcf,gca,figsize,2,tvec,PRE,STIM,tedges);
                title(gca,strcat(' sf ',num2str(sf(s))),'FontSize',12);
                xlim([PRE POST]); set(ph,'FaceColor','k');
                hold on;
                
                [zF1, f,freqrange,frid,pvec,meanpower,stdpower,powF1] = myzF1(resp,binsize,tf);
                
                fh=findall(0,'type','figure');
                if isempty(fh),
                    figPS = 1;
                else
                    figPS = fh(1).Number + 1;
                end
                if s == 1,
                    figure(figPS), psdN = figPS;
                else
                    figure(psdN),
                end
                subplot(ceil(length(sf)./2),2,s), plot(f(f>0 & f<30),pvec(f>0 & f<30),'-o'); vline(tf);
                set(gca,'Box','off','TickDir','out','FontSize',16);
                title(gca,strcat('SF ',num2str(sf(s))));
                setFigure(gcf,gca,figsize,0);
                if s == length(sf), xlabel('sf (cyc/deg)'), end
                if s == floor(length(sf)./2), ylabel('zF1'), end
                
                zF1s(s)      = mean(zF1,2);
                stdzF1s(s)   = std(zF1,0,2)./(sqrt(size(zF1,2)));
                trzF1s{s}    = zF1;
                
                spkrTF(s)    = mean(pvec(frid,:),2);
                stdspkrTF(s) = std(pvec(frid,:),0,2)./sqrt(size(pvec,2));
                trspkrTF{s}  = pvec(frid,:);
                
                F1s(s)       = mean(pvec(frid,:)./(pvec(1,:)),2);
                stdF1s(s)    = std(pvec(frid,:)./(pvec(1,:)),0,2)./(sqrt(size(pvec,2)));
                trF1s{s}     = pvec(frid,:)./(pvec(1,:));
                
                F0s(s)       = mean(mean(resp,1),2);
                stdF0s(s)    = std(mean(resp,1),0,2)./(sqrt(size(resp,1)));
                trF0s{s}     = mean(resp,1);
                
                realSF(s)    = sf(s);
                
                clear thesetrs;
            end
            
            figure,set(gcf,'Color','w','Position',[1 1 1000 400],'PaperPositionMode','auto'),
            
            subplot(1,3,1),
            h = plot(realSF,F1s,'-o','LineWidth',2); rgbcol = get(h,'Color'); hold on;
            errorbar(realSF,F1s,stdF1s,'LineStyle','none','Color',rgbcol);
            [~,mx] = max(F1s); [~,mn] = min(F1s); vline(realSF([mn mx]));
            set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('F1');title(gca,sprintf('max: %0.2f, min: %0.2f',realSF(mx),realSF(mn)),'fontweight','normal','fontsize',14);
            fprintf('\nF1 max: %0.2f, min: %0.2f\n',realSF(mx),realSF(mn))
            
            subplot(1,3,2),
            plot(realSF,F0s,'-o','LineWidth',2); hold on;
            errorbar(realSF,F0s,stdF0s,'LineStyle','none','Color',rgbcol);
            [~,mx] = max(F0s); [~,mn] = min(F0s); vline(realSF([mn mx]));
            set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('F0');xlabel('sf(cyc/deg)'); title(gca,sprintf('max: %0.2f, min: %0.2f',realSF(mx),realSF(mn)),'fontweight','normal','fontsize',14);
            fprintf('\nF0 max: %0.2f, min: %0.2f\n',realSF(mx),realSF(mn))
            
            subplot(1,3,3),plot(realSF,zF1s,'-o','LineWidth',2); hold on;
            errorbar(realSF,zF1s,stdzF1s,'LineStyle','none','Color',rgbcol);
            [~,mx] = max(zF1s); [~,mn] = min(zF1s); vline(realSF([mn mx]));
            set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('zF1'); title(gca,sprintf('max: %0.2f, min: %0.2f',realSF(mx),realSF(mn)),'fontweight','normal','fontsize',14);
            fprintf('\nzF1 max: %0.2f, min: %0.2f\n',realSF(mx),realSF(mn))
            
            if dofit == 1
                
                figure, set(gcf,'Color','w','Position',[500 500 1000 300],'PaperPositionMode','auto')
                subplot(1,3,1)
                y = [zF1s]';
                x = [realSF]';
                f = fit(x,y,'smoothingspline'); [~,mn] = min(zF1s); [~,mx] = max(zF1s);
                plot(f,x,y); hold on
                h = plot(x,y,'o'); rgbcol = get(h,'Color'); set(h,'MarkerFaceColor',rgbcol); hold on
                errorbar(x,y,[stdzF1s],'LineStyle','none','Color',rgbcol);
                axis tight; axis square; box off; legend('off'); vline(realSF([mn mx]));
                set(gca,'TickDir','out','FontSize',16); xlabel('sf (cyc/deg)'); ylabel('zF1');
                title(gca,sprintf('max: %0.2f, min: %0.2f',realSF(mx),realSF(mn)),'fontweight','normal','fontsize',14);
                
                
                subplot(1,3,2)
                clear x y
                y = [F1s]';
                x = [realSF]';[~,mn] = min(F1s); [~,mx] = max(F1s);
                f = fit(x,y,'smoothingspline');
                plot(f,x,y); hold on
                h = plot(x,y,'o'); rgbcol = get(h,'Color'); set(h,'MarkerFaceColor',rgbcol); hold on
                errorbar(x,y,[stdF1s],'LineStyle','none','Color',rgbcol);
                axis tight; axis square; box off; legend('off'); vline(realSF([mn mx]));
                set(gca,'TickDir','out','FontSize',16); xlabel('sf (cyc/deg)'); ylabel('F1');
                title(gca,sprintf('max: %0.2f, min: %0.2f',realSF(mx),realSF(mn)),'fontweight','normal','fontsize',14);
                
                subplot(1,3,3)
                clear x y
                y = [F0s]';
                x = [realSF]';
                f = fit(x,y,'smoothingspline');
                plot(f,x,y); hold on; [~,mn] = min(F0s); [~,mx] = max(F0s);
                h = plot(x,y,'o'); rgbcol = get(h,'Color'); set(h,'MarkerFaceColor',rgbcol); hold on
                errorbar(x,y,[stdF0s],'LineStyle','none','Color',rgbcol);
                axis tight; axis square; box off; legend('off'); vline(realSF([mn mx]));
                set(gca,'TickDir','out','FontSize',16); xlabel('sf (cyc/deg)'); ylabel('F0');
                title(gca,sprintf('max: %0.2f, min: %0.2f',realSF(mx),realSF(mn)),'fontweight','normal','fontsize',14);
                
            end
            
            if saveon == 1
                fh=findall(0,'type','figure');
                
                for f  = length(fh):-1:1
                    get(figure(fh(f).Number))
                    export_fig(strcat(savedrname,unq_unit),'-pdf','-nocrop','-append')
                    
                end
            end
            
        case 'cone'
            
            figure,set(gcf,'Color','w'),
            text(0.1,0.7,BRdatafile,'interpreter','none','FontSize',20,'FontName','arial');
            hold on;
            text(0.1,0.4,unq_unit,'interpreter','none','FontSize',20,'FontName','arial');
            hold on;
            text(0.1,0.1,strcat(' run on ',date),'interpreter','none','FontSize',20,'FontName','arial');
            set(gca,'XColor','w','YColor','w')
            
            path = nanunique(STIM.pathw);
            
            ct = 0;
            for c = 1:length(path)
                ct = ct + 1;
                thesetrs = r_STIM.pathw == path(c);
                tf       = STIM.temporal_freq(find(thesetrs,1,'first'));
                % raster:
                fh=findall(0,'type','figure');
                if isempty(fh),
                    figH.N    = 1;
                    figH.row  = length(path);
                    figH.col  = 1;
                    figH.subplot = ct;
                    rasterN = figH.N;
                else
                    if exist('rasterN')
                        figH.N = rasterN;
                    else
                        figH.N = fh(1).Number + 1;
                        
                    end
                    figH.row  = length(path);
                    figH.col  = 1;
                    figH.subplot = ct;
                    rasterN = figH.N;
                end
                
                plotSpikeRaster(logical(r_spktr(:,thesetrs)'),'PlotType','scatter');
                setFigure(gcf,gca,figsize,1,tvec,PRE,STIM)
                title(gca,strcat(' path ',num2str(path(c))),'FontSize',12)
                xlim([1 length([pre:post])]);
                
                % psth:
                fh=findall(0,'type','figure');
                if isempty(fh),
                    psthN = 1 ;
                else
                    
                    psthN = fh(1).Number + 1;
                end
                
                [tedges,resp,binsize] = myPSTH(length([pre:post]),PRE,thesetrs,r_spktintr,r_trls,fs);
                figure(psthN),
                subplot(length(path),1,ct)
                ph = bar(tedges(1:end-1),mean(resp(1:end-1,:),2),'histc');
                setFigure(gcf,gca,figsize,2,tvec,PRE,STIM,tedges);
                title(gca,strcat('path ',num2str(path(c))));
                xlim([PRE POST]); set(ph,'FaceColor','k');
                
                [zF1, f,freqrange,frid,pvec,meanpower,stdpower,powF1] = myzF1(resp,binsize,tf);
                
                fh=findall(0,'type','figure');
                if isempty(fh),
                    figPS = 1;
                else
                    figPS = fh(1).Number + 1;
                end
                if c == 1,
                    figure(figPS), psdN = figPS;
                else
                    figure(psdN),
                end
                subplot(length(path),1,c), plot(f(f>0 & f<30),pvec(f>0 & f<30),'-o'); vline(tf);
                set(gca,'Box','off','TickDir','out','FontSize',16);
                title(gca,strcat('path ',num2str(path(c))));
                setFigure(gcf,gca,figsize,0);
                
                zF1s(c)      = mean(zF1,2);
                stdzF1s(c)   = std(zF1,0,2)./(sqrt(size(zF1,2)));
                trzF1s{c}    = zF1;
                
                spkrTF(c)    = mean(pvec(frid,:),2);
                stdspkrTF(c) = std(pvec(frid,:),0,2)./sqrt(size(pvec,2));
                trspkrTF{c}  = pvec(frid,:);
                
                F1s(c)       = mean(pvec(frid,:)./(pvec(1,:)),2);
                stdF1s(c)    = std(pvec(frid,:)./(pvec(1,:)),0,2)./(sqrt(size(pvec,2)));
                trF1s{c}     = pvec(frid,:)./(pvec(1,:));
                
                F0s(c)       = mean(mean(resp,1),2);
                stdF0s(c)    = std(mean(resp,1),0,2)./(sqrt(size(resp,1)));
                trF0s{c}     = mean(resp,1);
                
                realPATH(c)    = path(c);
                
                clear thesetrs;
            end
            
            figure,set(gcf,'Color','w','Position',[1 1 1000 400],'PaperPositionMode','auto'),
            
            subplot(1,3,1),
            h = bar(realPATH,F1s,'FaceColor',[0.8 0.2 0.2]); rgbcol = get(h,'FaceColor'); hold on;
            errorbar(realPATH,F1s,stdF1s,'LineStyle','none','Color',rgbcol);
            [~,mx] = max(F1s); [~,mn] = min(F1s);
            set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('F1');title(gca,sprintf('max: %0.2f, min: %0.2f',realPATH(mx),realPATH(mn)),'fontweight','normal','fontsize',14);
            fprintf('\nF1 max: %0.2f, min: %0.2f\n',realPATH(mx),realPATH(mn))
            
            subplot(1,3,2),
            bar(realPATH,F0s,'FaceColor',rgbcol); hold on;
            errorbar(realPATH,F0s,stdF0s,'LineStyle','none','Color',rgbcol);
            [~,mx] = max(F0s); [~,mn] = min(F0s);
            set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('F0');xlabel('path'); title(gca,sprintf('max: %0.2f, min: %0.2f',realPATH(mx),realPATH(mn)),'fontweight','normal','fontsize',14);
            fprintf('\nF0 max: %0.2f, min: %0.2f\n',realPATH(mx),realPATH(mn))
            
            subplot(1,3,3),
            bar(realPATH,zF1s,'FaceColor',rgbcol); hold on;
            errorbar(realPATH,zF1s,stdzF1s,'LineStyle','none','Color',rgbcol);
            [~,mx] = max(zF1s); [~,mn] = min(zF1s);
            set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('zF1'); title(gca,sprintf('max: %0.2f, min: %0.2f',realPATH(mx),realPATH(mn)),'fontweight','normal','fontsize',14);
            fprintf('\nzF1 max: %0.2f, min: %0.2f\n',realPATH(mx),realPATH(mn))
            
        case 'disp'
            
            
            figure,set(gcf,'Color','w'),
            text(0.1,0.7,BRdatafile,'interpreter','none','FontSize',20,'FontName','arial');
            hold on;
            text(0.1,0.4,unq_unit,'interpreter','none','FontSize',20,'FontName','arial');
            hold on;
            text(0.1,0.1,strcat(' run on ',date),'interpreter','none','FontSize',20,'FontName','arial');
            set(gca,'XColor','w','YColor','w')
            
            
            dis  = nanunique(r_STIM.disparity);
            tf   = nanunique(r_STIM.temporal_freq);
            
            for d = 1:length(dis)
                
                thesetrs = r_STIM.disparity == dis(d);
                
                %raster:
                fh=findall(0,'type','figure');
                if isempty(fh),
                    figH.N = 1;
                    figH.subplot = d;
                    figH.row    = floor(length(dis)./2);
                    figH.col    = 2;
                    rasterN = figH.N;
                else
                    if exist('rasterN')
                        figH.N = rasterN;
                    else
                        figH.N = fh(1).Number + 1;
                        rasterN = figH.N;
                    end
                    figH.subplot = d;
                    figH.row    = floor(length(dis)./2);
                    figH.col    = 2;
                    
                end
                
                plotSpikeRaster(logical(r_spktr(:,thesetrs)'),'PlotType','scatter');
                setFigure(gcf,gca,figsize,1,tvec,PRE,STIM)
                title(gca,strcat(' disparity ',num2str(dis(d))),'FontSize',12)
                xlim([1 length([pre:post])]);
                
                % psth:
                [tedges,resp,binsize] = myPSTH(length([pre:post]),PRE,thesetrs,r_spktintr,r_trls,fs);
                fh=findall(0,'type','figure');
                if isempty(fh),
                    figP = 1;
                else
                    figP = fh(1).Number + 1;
                end
                if d == 1,
                    figure(figP), psthN = figP;
                else
                    figure(psthN),
                end
                
                subplot(floor(length(dis)./2),2,d)
                ph = bar(tedges(1:end-1),mean(resp(1:end-1,:),2),'histc');
                setFigure(gcf,gca,figsize,2,tvec,PRE,STIM,tedges);
                title(gca,strcat(' dis ',num2str(dis(d))),'FontSize',12);
                xlim([PRE POST]); set(ph,'FaceColor','k');
                hold on;
                
                [zF1, f,freqrange,frid,pvec,meanpower,stdpower,powF1] = myzF1(resp,binsize,tf);
                
                fh=findall(0,'type','figure');
                if isempty(fh),
                    figPS = 1;
                else
                    figPS = fh(1).Number + 1;
                end
                if d == 1,
                    figure(figPS), psdN = figPS;
                else
                    figure(psdN),
                end
                subplot(floor(length(dis)./2),2,d), plot(f(f>0 & f<30),pvec(f>0 & f<30),'-o'); vline(tf);
                set(gca,'Box','off','TickDir','out','FontSize',16);
                title(gca,strcat('DISPARITY ',num2str(dis(d))));
                setFigure(gcf,gca,figsize,0);
                if d == length(dis), xlabel('disparity (deg)'), end
                if d == floor(length(dis)./2), ylabel('zF1'), end
                
                zF1s(d)      = mean(zF1,2);
                stdzF1s(d)   = std(zF1,0,2)./(sqrt(size(zF1,2)));
                trzF1s{d}    = zF1;
                
                spkrTF(d)    = mean(pvec(frid,:),2);
                stdspkrTF(d) = std(pvec(frid,:),0,2)./sqrt(size(pvec,2));
                trspkrTF{d}  = pvec(frid,:);
                
                F1s(d)       = mean(pvec(frid,:)./(pvec(1,:)),2);
                stdF1s(d)    = std(pvec(frid,:)./(pvec(1,:)),0,2)./(sqrt(size(pvec,2)));
                trF1s{d}     = pvec(frid,:)./(pvec(1,:));
                
                F0s(d)       = mean(mean(resp,1),2);
                stdF0s(d)    = std(mean(resp,1),0,2)./(sqrt(size(resp,1)));
                trF0s{d}     = mean(resp,1);
                
                realDIS(d)    = dis(d);
                
                clear thesetrs;
            end
            
            figure,set(gcf,'Color','w','Position',[1 1 1000 400],'PaperPositionMode','auto'),
            
            subplot(1,3,1),
            h = plot(realDIS,F1s,'-o','LineWidth',2); rgbcol = get(h,'Color'); hold on;
            errorbar(realDIS,F1s,stdF1s,'LineStyle','none','Color',rgbcol);
            [~,mx] = max(F1s); [~,mn] = min(F1s); vline(dis([mn mx]));
            set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('F1');title(gca,sprintf('max: %u, min: %u',floor(dis(mn)),floor(dis(mx))),'fontweight','normal','fontsize',14);
            xlim([dis(1) dis(end)]);
            fprintf('\nF1: max: %u, min: %u\n',floor(dis(mn)),floor(dis(mx)));
            
            subplot(1,3,2),
            plot(realDIS,F0s,'-o','LineWidth',2); hold on;
            errorbar(realDIS,F0s,stdF0s,'LineStyle','none','Color',rgbcol);
            [~,mx] = max(F0s); [~,mn] = min(F0s); vline(dis([mn mx]));
            set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('F0');xlabel('disparity (deg)');title(gca,sprintf('max: %u, min: %u',floor(dis(mn)),floor(dis(mx))),'fontweight','normal','fontsize',14);
            xlim([dis(1) dis(end)]);
            fprintf('\nF0: max: %u, min: %u\n',floor(dis(mn)),floor(dis(mx)));
            
            subplot(1,3,3),plot(realDIS,zF1s,'-o','LineWidth',2); hold on;
            errorbar(realDIS,zF1s,stdzF1s,'LineStyle','none','Color',rgbcol);
            [~,mx] = max(zF1s); [~,mn] = min(zF1s); vline(dis([mn mx]));
            set(gca,'Box','off','TickDir','out','FontSize',16);
            ylabel('zF1'); title(gca,sprintf('max: %u, min: %u',floor(dis(mn)),floor(dis(mx))),'fontweight','normal','fontsize',14);
            xlim([dis(1) dis(end)]);
            fprintf('\nzF1: max: %u, min: %u\n',floor(dis(mn)),floor(dis(mx)));
            
    end
end
%%
if csdonly == 1
    
    switch paradigm
        
        case 'cint'
            
            % load LFP
            
            load(strcat(brdrname,'/',BRdatafile,'.lfp'),'-MAT'); 
     
            [ID,ids] = sortElectrodeLabels(NeuralLabels);
            chans = [1:24]; 
            lfp = LFP(:,ids(ID(2),:)); 
            lfp = lfp(:,chans); 
            csd = mod_iCSD(lfp')'; 
            csdchans = chans(2:end-1); 
            [trneural,cycneural] = triggerAnalogData(csd,PRE,POST,r_spkTPs,r_all_cyc,fs,per);
            
         
            unqc  = unique(STIM.contrast(STIM.contrast>=0));
            nunqc = unique(STIM.fixedc(STIM.fixedc>=0));
            tf    = nanunique(STIM.temporal_freq);
            ct    = 0;
            
            for c  = 1:length(unqc)
                for fc = 1:length(nunqc)
                    
                    thesetrs = r_STIM.contrast == unqc(c) & r_STIM.fixedc == nunqc(fc);
                    
                    thesedata = trneural(:,:,thesetrs);
                    mcsd   = mean(thesedata,3); %- repmat(mean(mean(trneural([1:abs(PRE)],:,thesetrs),1),3),[size(trneural,1) 1 1]);
                    stdcsd = std(thesedata,0,3);
                    
                    figure,set(gcf,'Color','w','PaperPositionMode','auto','Position',[2146           8         306         727]);
                    imagesc([WPRE:WPOST],csdchans,filterCSD(mcsd'))
                    cmap = colormap('jet');
                    colormap(flipud(cmap));
                    caxis([-1200 1200]); colorbar
                    set(gca,'TickDir','out','FontSize',15,'Box','off');
                    
                    rescal = rsEqChans(mcsd,2000);
                    figure, set(gcf,'Color','w','PaperPositionMode','auto','Position',[1865         -82         302         855]);
                    for ch = 1:size(mcsd,2)
                        plot([WPRE:WPOST],rescal(:,ch),'LineWidth',2,'Color',[.8 .2 .2])
                        hold on;
                        text(WPRE-300,rescal(1,ch),strcat('ch',num2str(csdchans(ch))));
                        hold on;
                    end
                    set(gca,'YTick',[],'YColor','w','TickDir','out','Box','off'); vline(0);
                    xlabel('t (ms)'); title(gca,strcat('DE ',num2str(unqc(c)),' NDE ',num2str(nunqc(fc))));
                    xlim([PRE POST]);
                    
                    figure,set(gcf,'Color','w','Position',[  1446         -38         700         800]); ct = 0;
                    for ch = 1:size(thesedata,2)
                        ct = ct + 1;
                        subplot(ceil(size(thesedata,2)./3),3,ct)
                        plot([WPRE:WPOST],mean(thesedata(:,ch,:),3));
                        set(gca,'Box','off','TickDir','out','FontSize',15);
                        xlim([WPRE WPOST]); h = title(gca,strcat(ch,num2str(csdchans(ch)))); set(h,'FontSize',8);
                        vline(0);
                    end
                    
                    figure, set(gcf,'Color','w','Position',[50 50 700 300]); 
                    h = plot([WPRE:WPOST],mean(thesedata(:,csdchans==15,:),3),'LineWidth',2); rgbcol = get(h,'Color'); 
                    hold on; 
                    errorbar([WPRE:50:WPOST],mean(thesedata(1:50:end,csdchans==15,:),3),std(thesedata(1:50:end,csdchans==15,:),0,3)./sqrt(size(thesedata,3)),'LineStyle','none','Color',rgbcol,'LineWidth',2); 
                    set(gca,'Box','off','TickDir','out','FontSize',15);
                    xlim([WPRE WPOST]);
                    vline(0); title(gca,'channel 15'); xlabel('t(ms)'); ylabel('nA/mm^3'); 
                    
                    for ch = 1:size(thesedata,2)
                        [pow,f] = myPSD(squeeze(thesedata(:,ch,:)));
                        frid            = find(abs(f-tf) == min(abs(f - tf)));
                        patTF(ch) = mean(pow(frid,:),2);
                        stdatTF(ch) = std(pow(frid,:),0,2)./size(pow,2);
                    end
                    figure,set(gcf,'Color','w','Position', [2726         -35         300         800]);
                    h = bar(csdchans,patTF); set(h,'FaceColor',getColor(6)); hold on, errorbar(csdchans,patTF,stdatTF,'LineStyle','none','Color',getColor(6));
                    set(gca,'TickDir','out','FontSize',15,'Box','off','XDir','reverse');
                    view(90,-90); xlim([csdchans(1)-1 csdchans(end)+1]); xlabel('channel'); ylabel('pow at TF (nA/mm^3');
                    
                    mresp(:,c,fc)   = mean(mcsd(abs(WPRE):end,:),1).*-1; 
                    stdresp(:,c,fc) = std(mean(thesedata(abs(WPRE):end,:,:),1),0,3)./(sqrt(size(thesedata,3))); 
                    
                end
            end  
           
           for chan = 1:size(mresp,1) 
            figure,set(gcf,'Color','w'); 
            for fc = 1:size(mresp,3)
            h = plot(unqc,mresp(chan,:,fc),'-o','LineWidth',2); 
            rgbcol = get(h,'Color'); 
            hold on; 
            errorbar(unqc,mresp(chan,:,fc),squeeze(stdresp(chan,:,fc)),'LineStyle','none','Color',rgbcol); 
            end
            set(gca,'Box','off','TickDir','out','FontSize',16); ylabel('-1*nA/mm3'), xlabel('contrast DE (%)'); 
            title(gca,strcat('chan',num2str(csdchans(chan)))); xlim([0 1]); 
           end
            
    end
end
%%
                    
%                     sizes = cellfun(@(C) size(C,2), r_all_cyc); 
%                     
%                     for cy = 1:min(sizes)
%                            clear mcsd
%                            mcsd = squeeze(mean(cycneural(:,:,cy,thesetrs),4)); 
%                        figure,set(gcf,'Color','w','PaperPositionMode','auto','Position',[1 1 300 700]);
%                     imagesc([0:floor(per./(fs/1000))],csdchans,filterCSD(mcsd'))
%                     cmap = colormap('jet');
%                     colormap(flipud(cmap));
%                     caxis([-1200 1200]); colorbar
%                     set(gca,'TickDir','out','FontSize',15,'Box','off'); 
%                         
%                     end
                    
            
% %%
% if analycyc == 1
%     clearvars -global
%     clear rasterN tvec svec psthN figH figP figPS
%     global figH
%     svec    = 0: per;
%     tvec    = svec./(fs/1000);
%     figsize = [1 1 500 800];
%     switch paradigm
%         
%         case 'cint'
%             
%             unqc  = unique(STIM.contrast(STIM.contrast>=0));
%             nunqc = unique(STIM.fixedc(STIM.fixedc>=0));
%             tf    = STIM.temporal_freq(1);
%             ct = 0;
%             for c  = 1:length(unqc)
%                 for fc = 1:length(nunqc)
%                     
%                     thesetrs = r_STIM.contrast == unqc(c) & r_STIM.fixedc == nunqc(fc);
%                     ct = ct + 1;
%                     % raster:
%                     fh=findall(0,'type','figure');
%                     if isempty(fh),
%                         figH.N = 1;
%                         figH.subplot = ct;
%                         figH.row    = length(unqc);
%                         figH.col    = length(nunqc);
%                         rasterN = figH.N;
%                     else
%                         if exist('rasterN')
%                             figH.N = rasterN;
%                         else
%                             figH.N = fh(1).Number + 1;
%                             rasterN = figH.N;
%                         end
%                         figH.subplot = ct;
%                         figH.row    = length(unqc);
%                         figH.col    = length(nunqc);
%                         
%                     end
%                     
%                     cyctrs = reshape([r_spkcyc(:,:,thesetrs)],[size(r_spkcyc,1) length(find(thesetrs))*size(r_spkcyc,2)]);
%                     
%                     plotSpikeRaster(logical(cyctrs)','PlotType','scatter');
%                     setFigure(gcf,gca,figsize,1,tvec,0,STIM)
%                     title(gca,strcat(' DE ',num2str(unqc(c)),' NDE ',num2str(nunqc(fc))),'FontSize',12)
%                     xlim([0 svec(end)+1]);
%                     
%                     % psth:
%                     triallen = length(svec);
%                     [tedges,resp,binsize] = myCYCPSTH(triallen,thesetrs,r_spkcyc,fs);
%                     fh=findall(0,'type','figure');
%                     if isempty(fh),
%                         figP = 1;
%                     else
%                         figP = fh(1).Number + 1;
%                     end
%                     if ct == 1,
%                         figure(figP), psthN = figP;
%                     else
%                         figure(psthN),
%                     end
%                     subplot(length(unqc),length(nunqc),ct)
%                     ph = bar(tedges(1:end-1),mean(resp(1:end-1,:),2),'histc');
%                     setFigure(gcf,gca,figsize,2,tvec,0,STIM,tedges);
%                     title(gca,strcat(' DE ',num2str(unqc(c)),' NDE ',num2str(nunqc(fc))),'FontSize',12)
%                     xlim([tvec(1) tvec(end)]); ylim([0 100]);
%                     set(ph,'FaceColor','k')
%                     hold on;
%                     
%                     [zF1, f,freqrange,frid,pvec,meanpower,stdpower,powF1] = myzF1(resp,binsize,tf);
%                     
%                     zF1s(c,fc)      = mean(zF1,2);
%                     stdzF1s(c,fc)   = std(zF1,0,2)./(sqrt(size(zF1,2)));
%                     trzF1s{c,fc}    = zF1;
%                     
%                     spkrTF(c,fc)    = mean(pvec(frid,:),2);
%                     stdspkrTF(c,fc) = std(pvec(frid,:),0,2)./sqrt(size(pvec,2));
%                     trspkrTF{c,fc}  = pvec(frid,:);
%                     
%                     F1s(c,fc)       = mean(pvec(frid,:)./(pvec(1,:)),2);
%                     stdF1s(c,fc)    = std(pvec(frid,:)./(pvec(1,:)),0,2)./(sqrt(size(pvec,2)));
%                     trF1s{c,fc}     = pvec(frid,:)./(pvec(1,:));
%                     
%                     F0s(c,fc)       = mean(mean(resp,1),2);
%                     stdF0s(c,fc)    = std(mean(resp,1),0,2)./(sqrt(size(resp,1)));
%                     trF0s{c,fc}     = mean(resp,1);
%                     
%                 end
%             end
%             
%             
%             [pvals] = runCRFstats(unqc,nunqc,trzF1s); fprintf('\n zF1 p = %f\n',pvals(1));clear pvals
%             [pvals] = runCRFstats(unqc,nunqc,trspkrTF); fprintf('\n s @ TF p = %f\n',pvals(1));clear pvals
%             [pvals] = runCRFstats(unqc,nunqc,trF1s); fprintf('\n F1 p = %f\n',pvals(1));clear pvals
%             [pvals] = runCRFstats(unqc,nunqc,trF0s); fprintf('\n F0 p = %f\n',pvals(1));clear pvals
%             
%             
%             figure, set(gcf,'Color','w','Position',[1 1 800 350],'PaperPositionMode','auto');
%             subplot(1,4,1)
%             for fc = 1:length(nunqc)
%                 plot(unqc,zF1s(:,fc),'-o','LineWidth',2);
%                 hold on;
%                 errorbar(unqc,zF1s(:,fc),stdzF1s(:,fc),'LineStyle','none');
%                 hold on;
%             end
%             xlim([0 1]); set(gca,'Box','off','TickDir','out','FontSize',16);
%             ylabel('zF1'); xlabel('contrast DE (%)');
%             title(gca,'cycles');
%             
%             subplot(1,4,2)
%             for fc = 1:length(nunqc)
%                 plot(unqc,spkrTF(:,fc),'-o','LineWidth',2);
%                 hold on;
%                 errorbar(unqc,spkrTF(:,fc),stdspkrTF(:,fc),'LineStyle','none');
%                 hold on;
%             end
%             xlim([0 1]); set(gca,'Box','off','TickDir','out','FontSize',16);
%             ylabel('spks/s at TF'); xlabel('contrast DE (%)');
%             
%             subplot(1,4,3)
%             for fc = 1:length(nunqc)
%                 plot(unqc,F1s(:,fc),'-o','LineWidth',2);
%                 hold on;
%                 errorbar(unqc,F1s(:,fc),stdF1s(:,fc),'LineStyle','none');
%                 hold on;
%             end
%             xlim([0 1]); set(gca,'Box','off','TickDir','out','FontSize',16);
%             ylabel('F1'); xlabel('contrast DE (%)');
%             
%             subplot(1,4,4)
%             for fc = 1:length(nunqc)
%                 plot(unqc,F0s(:,fc),'-o','LineWidth',2);
%                 hold on;
%                 errorbar(unqc,F0s(:,fc),stdF0s(:,fc),'LineStyle','none');
%                 hold on;
%             end
%             xlim([0 1]); set(gca,'Box','off','TickDir','out','FontSize',16);
%             ylabel('F0 (spk/s)'); xlabel('contrast DE (%)');
%             
%             
%             
%             
%             
%     end
% end
% 
% 
% 
% %     %%
% %
% %             if crfs
% %
% %                 % spike rate:
% %                 for c = 1:length(unqc)
% %                     for fc = 1:length(nunqc)
% %
% %                         spkrate{c,fc} = double(sum(spkdata{c,fc}(abs(pre):end,:),1))./...
% %                             length(tvec)./(double(NEV.MetaTags.SampleRes)/1000)...%divided by window length at 1kHz sampling rate
% %                             .*(1000./(length(tvec)./(double(NEV.MetaTags.SampleRes)/1000))); % scale window length to 1 sec to get spks/sec;
% %
% %                         mspk(c,fc)   = nanmean(spkrate{c,fc});
% %                         semspk(c,fc) = std(spkrate{c,fc},0,2)./sqrt(size(spkrate{c,fc},2));
% %
% %                     end
% %                 end
% %
% %                 clear Data prd
% %                 figure,
% %
% %                 for fc = 1:length(nunqc)
% %                     global Data
% %                     Data = [unqc.*100 mspk(:,fc)]';
% %                     [Gr Gc q s b fbest] = fitCRdata;
% %                     prcon = [0:100];
% %                     for pc = 1:length(prcon)
% %                         prd(pc) = Gr*[prcon(pc)^(q+s)]/[prcon(pc)^q + Gc^q]+b; % prediction
% %                     end
% %                     plot(prcon,prd,'Color',getColor(colid(fc)),'LineWidth',2);
% %                     hold on;
% %                     plot(unqc*100,mspk(:,fc),'o','Color',getColor(colid(fc)),'MarkerFaceColor',getColor(colid(fc)),'LineWidth',2);
% %                     hold on;
% %                     errorbar(unqc*100,mspk(:,fc),semspk(:,fc),'LineStyle','none','Color',getColor(colid(fc)));
% %                     hold on;
% %                 end
% %                 setFigure(gcf,gca,[500 500 300 500]);
% %                 set(gca,'LineWidth',2,'XLim',[0 100],'FontSize',16,'XScale','log');
% %                 hold on
% %                 xlabel('% contrast DE'); ylabel('spks/s');
% %                 title(gca,strcat(BRdatafile,'chan=',elabel,' unit=',num2str(unit)),'interpreter','none','FontSize',10);
% %                 lc = strcat(num2str(nunqc*100),'%ND');
% %                 legend(gca,lc,'Location','best');
% %
% %                 h1 = gcf;
% %                 h2 = figure;
% %                 objects = allchild(h1);
% %                 copyobj(get(h1,'children'),h2);
% %                 set(gcf,'Color','w','PaperPositionMode','auto','Position',[500 500 300 500]);
% %                 set(gca,'XScale','lin');
% %
% %                 % anova (2x1): Compare mean responses across monocular condition against same aross dichoptic condition
% %                 y = reshape(mspk,size(mspk,1)*size(mspk,2),1);
% %                 g1 = [repmat(1,size(mspk,1),1); repmat(2,size(mspk,1),1)];
% %                 anovan(y,{g1})
% %
% %                 %anova (5x2): Compare responses at each contrast level
% %                 %and for each ocular conditions
% %                 y = []; g1 = []; g2= [];
% %                 for c = 1:length(unqc)
% %                     for fc = 1:length(nunqc)
% %                         y  = [y spkrate{c,fc}];
% %                         g1 = [g1; repmat(fc,length(spkrate{c,fc}),1)]; % ocular condition
% %                         g2 = [g2; repmat(c,length(spkrate{c,fc}),1)];  % contrast level
% %                     end
% %                 end
% %                 anovan(y,{g1,g2})
% %
% %                 % spike rate, normalize between 0 and 1:
% %
% %                 for c = 1:length(unqc)
% %                     for fc = 1:length(nunqc)
% %                         mns(c,fc) = min(spkrate{c,fc});
% %                         mxs(c,fc) = max(spkrate{c,fc});
% %                     end
% %                 end
% %                 mn = min(min(mns));
% %                 mx = max(max(mxs));
% %                 clear normspk;
% %                 for c = 1:length(unqc)
% %                     for fc = 1:length(nunqc)
% %                         normspk{c,fc} = (spkrate{c,fc} - mn)./(mx - mn);
% %                         mnormspk(c,fc) = nanmean(normspk{c,fc});
% %                         semnormspk(c,fc) = std(normspk{c,fc},0,2)./sqrt(size(normspk{c,fc},2));
% %                     end
% %                 end
% %
% %                 clear Data prd
% %                 figure,
% %                 for fc = 1:length(nunqc)
% %                     global Data
% %
% %                     Data = [unqc.*100 mnormspk(:,fc)]';
% %                     [Gr Gc q s b fbest] = fitCRdata;
% %                     prcon = [0:100];
% %                     for pc = 1:length(prcon)
% %                         prd(pc) = Gr*[prcon(pc)^(q+s)]/[prcon(pc)^q + Gc^q]+b; % prediction
% %                     end
% %                     plot(prcon,prd,'Color',getColor(colid(fc)),'LineWidth',2);
% %                     hold on;
% %                     plot(unqc*100,mnormspk(:,fc),'o','Color',getColor(colid(fc)),'MarkerFaceColor',getColor(colid(fc)),'LineWidth',2);
% %                     hold on;
% %                     errorbar(unqc*100,mnormspk(:,fc),semnormspk(:,fc),'LineStyle','none','Color',getColor(colid(fc)));
% %                     hold on;
% %                 end
% %                 set(gcf,'Color','w','PaperPositionMode','auto','Position',[500 500 300 500]);
% %                 set(gca,'Box','off','TickDir','out','LineWidth',2,'XLim',[0 100],'FontSize',16,'XScale','log');
% %                 hold on
% %                 xlabel('% contrast DE'); ylabel('normalized spiking response');
% %                 title(gca,strcat(BRdatafile,'chan=',elabel,' unit=',num2str(unit)),'interpreter','none','FontSize',10);
% %                 lc = strcat(num2str(nunqc*100),'%ND');
% %                 legend(gca,lc,'Location','best');
% %
% %                 h1 = gcf;
% %                 h2 = figure;
% %                 objects = allchild(h1);
% %                 copyobj(get(h1,'children'),h2);
% %                 set(gcf,'Color','w','PaperPositionMode','auto','Position',[500 500 300 500]);
% %                 set(gca,'XScale','lin');
% %
% %             end
% %
% %             if photo_yes
% %
% %                 % first convolve signal and look at whole trial:
% %                 for c = 1:length(unqc)
% %                     for fc = 1:length(nunqc)
% %
% %                         binsize       = 35; % ms
% %                         hspkdens      = calcSpikeDensity(spkdata{c,fc},binsize*(double(NEV.MetaTags.SampleRes)/1000));
% %                         spkdens{c,fc} = hspkdens;
% %
% %                     end
% %                 end
% %
% %                 % normalize across all trials:
% %
% %                 it = 0;
%                 for fc = 1:length(nunqc)
%                     for c = 1:length(unqc)
%                         it = it + 1;
%                         cmin(it) = min(min(spkdens{c,fc}));
%                         cmax(it) = max(max(spkdens{c,fc}));
%                     end
%                 end
%                 maxn = max(cmax); minn = (min(cmin));
%                 for fc = 1:length(nunqc)
%                     for c = 1:length(unqc)
%                         for tr = 1:size(spkdens{c,fc},2)
%                             nspkdens{c,fc}(:,tr) = (spkdens{c,fc}(:,tr) - minn)./(maxn-minn);
%                         end
%
%                     end
%                 end
%
%                 % whole trial:
%                 figure, set(gcf,'Color','w','Position',[1 1 1050 700],'PaperPositionMode','auto');
%                 sp = 0 ;
%                 for c = 1:length(unqc)
%                     sp = sp + 1;
%                     for fc = 1:length(nunqc)
%
%                         subplot(length(unqc),1,sp)
%                         if fc == 1
%                         plot(tvec,mean(nspkdens{c,fc},2),'LineWidth',2,'Color',[1 0 0]);
%                         else
%                              plot(tvec,mean(nspkdens{c,fc},2),'LineWidth',2,'Color',[0 0 1]);
%                         end
%                         hold on;
%                         legend(gca,sprintf('N %u',round(nunqc(fc).*100)),'Location','best','FontSize',8);
%                         title(gca,sprintf('D %u',round(unqc(c).*100)));
%                         set(gca,'Box','off','TickDir','out','FontSize',10);
%                         legend boxoff
%                         xlabel('t (ms)'); ylim([0.2 0.6]);
%                     end
%                     hold on;
%                 end
%
%                 % then do above analyses on individual cycles:
%                 for c = 1:length(unqc)
%                     for fc = 1:length(nunqc)
%
%                         spkvec = spkcyc{c,fc};
%                         for cy = 1:size(spkvec,2)
%
%                             hspkdenscy(:,:,cy) = calcSpikeDensity(squeeze(spkvec(:,cy,:)),binsize*(double(NEV.MetaTags.SampleRes)/1000));
%
%                         end
%
%                         spkdenscy{c,fc} = permute(hspkdenscy,[1 3 2]);
%
%                     end
%                 end
%
%                 % normalize across all trials:
%                 it = 0;
%                 for fc = 1:length(nunqc)
%                     for c = 1:length(unqc)
%                         it = it + 1;
%                         cmin(it) = min(min(nanmean(spkdenscy{c,fc},2)));
%                         cmax(it) = max(max(nanmean(spkdenscy{c,fc},2)));
%                     end
%                 end
%                 maxn = max(cmax); minn = (min(cmin));
%                 for fc = 1:length(nunqc)
%                     for c = 1:length(unqc)
%                         for tr = 1:size(spkdens{c,fc},2)
%                             for cy = 1:size(spkdenscy{c,fc},2)
%                                 nspkdenscy{c,fc}(:,cy,tr) = (spkdenscy{c,fc}(:,cy,tr) - minn)./(maxn-minn);
%                             end
%                         end
%                     end
%                 end
%
%                 % single cycle:
%                 figure, set(gcf,'Color','w','Position',[1 1 1050 700],'PaperPositionMode','auto'); sp = 0;
%                 for c = 1:length(unqc)
%                     for fc = 1:length(nunqc)
%                         sp = sp + 1;
%                         subplot(length(unqc),length(nunqc),sp)
%                         plot(tcyc,nanmean(nanmean(nspkdenscy{c,fc},2),3),'LineWidth',2);
%                         hold on;
%                         legend(gca,sprintf('N %u',round(nunqc(fc).*100)),'Location','best','FontSize',8);
%                         title(gca,sprintf('D %u',round(unqc(c).*100)));
%                         set(gca,'Box','off','TickDir','out','FontSize',10);
%                         legend boxoff
%                         xlabel('t (ms)');
%                         axis tight
%                     end
%                 end
%
%
%             end
%
%
%         case 'rfor'
%
%             oris  = unique(STIM.tilt(STIM.tilt>=0));
%             tfs   = unique(STIM.temporal_freq(STIM.temporal_freq>=0));
%
%
%
%         case 'rfsf'
%
%             sfs  = unique(STIM.sf(STIM.sf>=0));
%
%     end
%
%
%
% if csd == 1
%
%     % compute CSD
%     CSD = mod_iCSD(LFP')';
%
%     % time window variables:
%     pre  = PRE.*double((lfpFs/1000));
%     post = POST*double((lfpFs/1000));
%
%     svec = pre:post;
%     tvec = svec./double((lfpFs/1000));
%
%     if strcmp(BRdatafile,'161004_E_pmk001')
%         load('/Users/kaciedougherty/Documents/neurophysdata/161004_E/161004_E_pmk001_driftingGRATINGRECORD0001.mat','-MAT');
%         oldGratingrecord = GRATINGRECORD; clear GRATINGRECORD;
%
%         load('/Users/kaciedougherty/Documents/neurophysdata/161004_E/161004_E_pmk001_driftingGRATINGRECORD0046.mat','-MAT');
%         newGratingrecord = GRATINGRECORD; clear GRATINGRECORD;
%         GRATINGRECORD = cat(2,oldGratingrecord,newGratingrecord);
%         pmktr = [GRATINGRECORD(realtr).path] ;
%     end
%
%     if  strcmp(BRdatafile,'161003_E_pmk002')
%
%         load('/Users/kaciedougherty/Documents/neurophysdata/161003_E/161003_E_pmk002_driftingGRATINGRECORD0001.mat','-MAT');
%         oldGratingrecord = GRATINGRECORD; clear GRATINGRECORD;
%
%         load('/Users/kaciedougherty/Documents/neurophysdata/161003_E/161003_E_pmk002_driftingGRATINGRECORD0046.mat','-MAT');
%         newGratingrecord = GRATINGRECORD; clear GRATINGRECORD;
%         GRATINGRECORD = cat(2,oldGratingrecord,newGratingrecord);
%         pmktr = [GRATINGRECORD(realtr).path] ;
%
%     end
%
%     if  strcmp(BRdatafile,'161003_E_pmk003')
%
%         load('/Users/kaciedougherty/Documents/neurophysdata/161003_E/161003_E_pmk003_driftingGRATINGRECORD0001.mat','-MAT');
%         oldGratingrecord = GRATINGRECORD; clear GRATINGRECORD;
%
%         load('/Users/kaciedougherty/Documents/neurophysdata/161003_E/161003_E_pmk003_driftingGRATINGRECORD0046.mat','-MAT');
%         newGratingrecord = GRATINGRECORD; clear GRATINGRECORD;
%         GRATINGRECORD = cat(2,oldGratingrecord,newGratingrecord);
%         pmktr = [GRATINGRECORD(realtr).path] ;
%
%     end
%
%     if lfpFs == 1000, oldspkTPs = spkTPs; spkTPs = ceil(spkTPs./(double(NEV.MetaTags.SampleRes)./lfpFs)); end
%     for r = 1:length(spkTPs)
%         if spkTPs(r,1) == 0
%             continue
%         else
%
%             refwin = spkTPs(r,1) + pre : spkTPs(r,1) + post;
%             if pmktr(r) == 1
%                 pcsd(:,:,r) = CSD(refwin,:);
%             elseif pmktr(r) == 2
%                 mcsd(:,:,r) = CSD(refwin,:);
%             else
%                 kcsd(:,:,r) = CSD(refwin,:);
%             end
%
%         end
%
%     end
%     %%
%     scfact = 5000;
%     figure, set(gcf,'Color','w','PaperPositionMode','auto','Position',[1 1 800 800]);
%     chans = 1:22;
%     subplot(1,3,1)
%     mpcsd = rsEqChans(nanmean(pcsd,3),scfact);
%     for chan = 1:length(chans)
%         plot(tvec,mpcsd(:,chans(chan)),'k');
%         hold on;
%         text(pre-175,mpcsd(1,chan),num2str(chans(chan)));
%         hold on;
%     end
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mpcsd(:,chans(end))) max(mpcsd(:,chans(1)))]);
%
%     subplot(1,3,2)
%     mmcsd = rsEqChans(nanmean(mcsd,3),scfact);
%     for chan = 1:length(chans)
%         plot(tvec,mmcsd(:,chans(chan)),'k');
%         hold on;
%         text(pre-175,mmcsd(1,chan),num2str(chans(chan)));
%         hold on;
%     end
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mmcsd(:,chans(end))) max(mmcsd(:,chans(1)))]);
%
%     subplot(1,3,3)
%     mkcsd = rsEqChans(nanmean(kcsd,3),scfact);
%     for chan = 1:length(chans)
%         plot(tvec,mkcsd(:,chans(chan)),'k');
%         hold on;
%         text(pre-175,mkcsd(1,chan),num2str(chans(chan)));
%         hold on;
%     end
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%
%     %%
%     figure, set(gcf,'Color','w','PaperPositionMode','auto','Position',[1 1 800 800]);
%     subplot(1,3,1)
%     plot(tvec,mpcsd,'r');
%     hold on;
%     plot(tvec,mmcsd,'k');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     subplot(1,3,2);
%     plot(tvec,mpcsd,'r');
%     hold on;
%     plot(tvec,mkcsd,'b');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     subplot(1,3,3);
%     plot(tvec,mmcsd,'k');
%     hold on;
%     plot(tvec,mkcsd,'b');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%     %%
%     clear mpcsd mkcsd mmcsd
%     scfact = 4000;
%     mpcsd = rsEqChans(nanmean(pcsd,3),scfact);
%     mmcsd = rsEqChans(nanmean(mcsd,3),scfact);
%     mkcsd = rsEqChans(nanmean(kcsd,3),scfact);
%     chans = [1:4];
%     figure, set(gcf,'Color','w','PaperPositionMode','auto','Position',[1 1 800 800]);
%     subplot(1,3,1)
%     plot(tvec,mpcsd(:,chans),'r');
%     hold on;
%     plot(tvec,mmcsd(:,chans),'k');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     subplot(1,3,2);
%     plot(tvec,mpcsd(:,chans),'r');
%     hold on;
%     plot(tvec,mkcsd(:,chans),'b');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     subplot(1,3,3);
%     plot(tvec,mmcsd(:,chans),'k');
%     hold on;
%     plot(tvec,mkcsd(:,chans),'b');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     chans = [5:9];
%     figure, set(gcf,'Color','w','PaperPositionMode','auto','Position',[1 1 800 800]);
%     subplot(1,3,1)
%     plot(tvec,mpcsd(:,chans),'r');
%     hold on;
%     plot(tvec,mmcsd(:,chans),'k');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     subplot(1,3,2);
%     plot(tvec,mpcsd(:,chans),'r');
%     hold on;
%     plot(tvec,mkcsd(:,chans),'b');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     subplot(1,3,3);
%     plot(tvec,mmcsd(:,chans),'k');
%     hold on;
%     plot(tvec,mkcsd(:,chans),'b');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     chans = [10:14];
%     figure, set(gcf,'Color','w','PaperPositionMode','auto','Position',[1 1 800 800]);
%     subplot(1,3,1)
%     plot(tvec,mpcsd(:,chans),'r');
%     hold on;
%     plot(tvec,mmcsd(:,chans),'k');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     subplot(1,3,2);
%     plot(tvec,mpcsd(:,chans),'r');
%     hold on;
%     plot(tvec,mkcsd(:,chans),'b');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     subplot(1,3,3);
%     plot(tvec,mmcsd(:,chans),'k');
%     hold on;
%     plot(tvec,mkcsd(:,chans),'b');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     chans = [15:19];
%     figure, set(gcf,'Color','w','PaperPositionMode','auto','Position',[1 1 800 800]);
%     subplot(1,3,1)
%     plot(tvec,mpcsd(:,chans),'r');
%     hold on;
%     plot(tvec,mmcsd(:,chans),'k');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     subplot(1,3,2);
%     plot(tvec,mpcsd(:,chans),'r');
%     hold on;
%     plot(tvec,mkcsd(:,chans),'b');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     subplot(1,3,3);
%     plot(tvec,mmcsd(:,chans),'k');
%     hold on;
%     plot(tvec,mkcsd(:,chans),'b');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     chans = [20:22];
%     figure, set(gcf,'Color','w','PaperPositionMode','auto','Position',[1 1 800 800]);
%     subplot(1,3,1)
%     plot(tvec,mpcsd(:,chans),'r');
%     hold on;
%     plot(tvec,mmcsd(:,chans),'k');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     subplot(1,3,2);
%     plot(tvec,mpcsd(:,chans),'r');
%     hold on;
%     plot(tvec,mkcsd(:,chans),'b');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%
%     subplot(1,3,3);
%     plot(tvec,mmcsd(:,chans),'k');
%     hold on;
%     plot(tvec,mkcsd(:,chans),'b');
%     set(gca,'Box','off','TickDir','out','YTick',[],'YColor','w');
%     xlim([pre post]); ylim([min(mkcsd(:,chans(end))) max(mkcsd(:,chans(1)))]);
%     v =  vline(0); set(v,'Color','k','LineWidth',1);
%     %%
% end
