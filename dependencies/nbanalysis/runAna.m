function [STIM,sdftr,time] = runAna(BRdatafile,clustN)    

% extension     = 'nev';
% el            = 'eD01';
% 
% PRE           = -300;            % ms relative to stim onset
% POST          = 1000;         % ms relative to stim onset
% RASTER        = 1;
% photo_yes     = 0;
% 
% brdrname      = sprintf('/users/kaciedougherty/documents/neurophysdata/%s',BRdatafile(1:8));
% mldrname      = brdrname;
% refresh       = 85;
% flipeye       = 0;
% dofit         = 0;

%% load text file for condition info
    
%     if ~isempty(strfind(BRdatafile,'ori'))
%         ext = '.gRFORIDRFTGrating_di';
%     elseif ~isempty(strfind(BRdatafile,'disparity'))
%         ext = '.gDISPARITYDRFTGrating_di';
%     elseif ~isempty(strfind(BRdatafile,'rfsf'))
%         ext = '.gRFSFDRFTGrating_di';
%     elseif ~isempty(strfind(BRdatafile,'tfsf'))
%         ext = '.gTFSFDRFTGrating_di';
%     elseif  ~isempty(strfind(BRdatafile,'size'))
%         ext = '.gRFSIZEDRFTGrating_di';
%     elseif  ~isempty(strfind(BRdatafile,'cinteroc'))
%         ext = '.gCINTEROCDRFTGrating_di';
%     elseif ~isempty(strfind(BRdatafile,'cpatch'))
%         ext = '.gCPATCHDRFTGrating_di';
%     elseif ~isempty(strfind(BRdatafile,'pmk'))
%         ext = '.gPMKDRFTGrating_di';
%     elseif ~isempty(strfind(BRdatafile,'cone'))
% %         ext = '.gCONEDRFTGrating_di';
% %     elseif ~isempty(strfind(BRdatafile,'color'))
% %         ext = '.gCOLORFLICKERDRFTGrating_di';
% %     elseif ~isempty(strfind(BRdatafile,'bw'))
% %         ext = '.gBWFLICKERDRFTGrating_di';
% %     end
% %
% %
% %     if strfind(BRdatafile,'pmk')
% %         paradigm = 'pmk';
% %     else
% %         paradigm   = BRdatafile(10:end-3); paradigm = paradigm(1:4);
% %     end
% %
% %     badobs                   = [1];
% %     flag_RewardedTrialsOnly  = true;
% %     grating                  = readgDRFTGrating([mldrname filesep BRdatafile ext]);
% %
% %     %% file dependent switches/mods
% %     if strfind(BRdatafile,'cinteroc') & flipeye == 1
% %         hcontrast        = grating.contrast;
% %         grating.contrast = grating.fixedc;
% %         grating.fixedc   = hcontrast;
% %
% %     end
% 
% %     if  refresh ~= 85
% %         h_temporalfreq =  grating.temporal_freq;
% %         grating.temporal_freq = [];
% %         grating.temporal_freq = round( h_temporalfreq.*(refresh./85));
% %     end
% %
% % if strcmp(paradigm,'cone') && ~isfield(grating,'path') && strcmp(BRdatafile,'161214')
% %     % forgot to add path field initially in grating text file (12/14/16)
% %     currentdir = pwd;
% %     cd(brdrname);
% %
% %     finfo = dir(strcat(BRdatafile,'*','GRATINGRECORD','*','.mat'));
% %     [~,sid] = sort([finfo(:).datenum],'ascend');
% %     hGRATINGRECORD = [];
% %     for fl = 1:length(finfo)
% %
% %         load(finfo(sid(fl)).name);
% %
% %         hGRATINGRECORD = [hGRATINGRECORD GRATINGRECORD];
% %         clear GRATINGRECORD;
% %     end
% %     GRATINGRECORD = hGRATINGRECORD; clear hGRATINGRECORD;
% %     grating.path = [];
% %     for tr = 1:max(grating.trial)
% %         grating.path = [grating.path GRATINGRECORD(tr).path GRATINGRECORD(tr).path];
% %     end
% % end
% %%
% % sort/pick trials [before iterating unit]
% if any(strfind(BRdatafile,'cone')) | any(strfind(BRdatafile,'color'))
%     stimfeatures = {...
%         'tilt'...
%         'sf'...
%         'contrast'...
%         'fixedc'...
%         'diameter'...
%         'eye'...
%         'oridist'...
%         'phase'...
%         'temporal_freq'...
%         'pathw'};
%     grating.path = [];
%     grating.path([strcmp(grating.pathw,'LM')])  = 1;
%     grating.path([strcmp(grating.pathw,'S')])   = 2;
%     grating.path([strcmp(grating.pathw,'LMo')]) = 3;
%     grating.path([strcmp(grating.pathw,'ach')]) = 4;
% else
%     stimfeatures = {...
%         'tilt'...
%         'sf'...
%         'contrast'...
%         'fixedc'...
%         'diameter'...
%         'eye'...
%         'oridist'...
%         'phase'...
%         'temporal_freq'};
% end
% 
% clear(stimfeatures{:})
    %% load digital codes and neural data:
%     filename = fullfile(brdrname,BRdatafile);
%     
%     % check if file exist and load NEV
%     if exist(strcat(filename,'.nev'),'file') == 2;
%         NEV = openNEV(strcat(filename,'.nev'),'read','overwrite');
%         labels = {NEV.ElectrodesInfo.ElectrodeLabel};
%         for i = 1:length(labels)
%             if strfind(labels{i}',el)
%                 idxlab = i;
%             end
%         end
%         pin = NEV.ElectrodesInfo(idxlab).ElectrodeID;
%     else
%         error('the following file does not exist\n%s.nev',filename);
%     end
%     
%     % get event codes from NEV
%     EventCodes   = NEV.Data.SerialDigitalIO.UnparsedData - 128;
%     EventTimes   = floor(NEV.Data.SerialDigitalIO.TimeStampSec .* 1000); %ms, to match 1kHz
%     EventSampels = NEV.Data.SerialDigitalIO.TimeStamp;
%     [pEvC, pEvT] = parsEventCodesML(EventCodes,EventSampels);
%     [STIM,realtr] = sortStimandTimeData(grating,pEvC,pEvT);
%     fs            = double(NEV.MetaTags.SampleRes);
%     if photo_yes
%         [STIM]    = photodiodeTrigger(filename,STIM,refresh);
%     end
%     %%
%         nevname = strcat('/Users/kaciedougherty/Documents/neurophysdata/kiloout/',BRdatafile,'_kilo_ss.mat');
%         load(nevname,'ss');
%         spk = ss; clear ss; 
%         kchans  = 1; 
%         kclust  = find(spk.clusterMap(:,1) == clustN); 
%         SPK     = spk.spikeTimes(spk.spikeClusters == clustN); 
%         wvform  = spk.spikeWaves(kchans,:,kclust); 
%                
%     %% get spikes for unit of interest:
% 
%         figure,set(gcf,'Color','w','Position',[1 1 500 300]);
%         plot(wvform,'Color','k','LineWidth',3);
%         %hold on;
%         %errbar = std(wvform,0,2)./sqrt(size(wvform,2));
%         %shadedErrorBar([1:size(wvform,1)],mean(wvform,2),errbar,{'-','color',[0 0 0]});
%         set(gca,'Box','off','TickDir','out','FontSize',16);
%         axis tight; ylabel('microV'); xlabel('samples');
%         title(gca,el);
% 
%     %% setup time vectors
%     
%     if photo_yes
%         WPRE   = PRE - ((1000/grating.temporal_freq(1))./2); % OR +
%         WPOST  = POST - ((1000/grating.temporal_freq(1))./2);
%     else
%         WPRE   = PRE;
%         WPOST  = POST;
%     end
%     
%     pre        = floor(WPRE);
%     post       = floor(WPOST);
%     svec       = floor((PRE.*(fs/1000))): floor((POST.*(fs/1000)));
%     tvec       = svec./(fs/1000);
%     time       = [PRE:POST];
%     %% trigger data neural data
%     
%     [sdf, sua, tm] = spk2sdf(SPK,fs);

% if photo_yes
%     triggerpts    = double(STIM.onsets_p);
%     
% else
%     triggerpts    = double(STIM.onsets);
% end
% if any((triggerpts./30 + post) > size(sdf,2))
%     keep = (triggerpts./30 + post) < size(sdf,2);
%     fn = fieldnames(STIM);
%     for f= 1:length(fn)
%         STIM.(fn{f}) = STIM.(fn{f})(keep);
%     end
%     triggerpts(triggerpts./30 > size(sdf,2)) = [];
% end
% sdftr              = squeeze(trigData(sdf',floor(triggerpts./30),-pre,post));


    
    if RASTER  == 1
        [spktr] = triggerPointData(floor(pre*30),floor(30*post),triggerpts,pin,NEV,SPK);
        spktr   = squeeze(spktr);
    else
        spktr   = [];
    end
    
    %%
    STIM.current_chan  = el;

        STIM.current_clust = kclust; 
    
    switch paradigm
        case 'rfor'
            %rforidrft
            anaOriDrft(STIM,sdftr,time,spktr)
            
        case 'rfsf'
            %rfsfdrft
            anaSFDrft(STIM,sdftr,time,spktr)
            
        case 'cone'
            %conedrft
            anaConeDrft(STIM,sdftr,time,spktr)
            
        case 'coni'
            %coneinterocdrft
            anaCinterocDrft(STIM,sdftr,time,spktr,dofit)
            
        case 'cint'
            %cinterocdrft
            
            anaCinterocDrft(STIM,sdftr,time,spktr,dofit)
            
        case 'rfsi'
            %rfsizedrft
            anaSizeDrft(STIM,sdftr,time,spktr)
            
        case 'colo'
            %colorflicker
            anaColorFlicker(STIM,sdftr,time,spktr)
            
        case 'bwfl'
            anaBWFlicker(STIM,sdftr,time,spktr)
            
        case 'disp'
            anaDisparityDrft(STIM,sdftr,time,spktr)
            
    end
    
