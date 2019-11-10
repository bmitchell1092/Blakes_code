function [trLFP,trSDF,wholesdf,wholelfp,trCSD,TTL,STIM] = getInjectionData(BRdatafile,el,chans,pre,post)
STIM            = []; 
TTL             = []; 
extension       = 'ns2';
NEV             = openNEV(strcat(BRdatafile,'.nev'),'read','overwrite','uV');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec .* 1000); %ms, to match 1kHz

if contains(BRdatafile,'evp')
    triggerpoints   = EventTimes(EventCodes == 23 | EventCodes == 25 | EventCodes == 27 | EventCodes == 29| EventCodes == 31);
else
    
    [pEvC, pEvT]    = parsEventCodesML(EventCodes,EventTimes);
    ext             = '.gRFORIGrating_di';
    grating         = readgGrating([BRdatafile ext]);
    [STIM,realtr]   = sortStimandTimeData(grating,pEvC,pEvT,'stim');
    triggerpoints   = STIM.onsets; 
end

LFP             = getLFP(BRdatafile,extension,el,'ascending'); LFP = LFP(:,chans); wholelfp = LFP;
trLFP           = trigData(LFP,triggerpoints,pre,post); 

CSD             = mod_iCSD(LFP')'; CSD = [nan(size(CSD,1),1) CSD nan(size(CSD,1), 1)]; 
trCSD           = trigData(CSD,triggerpoints,pre,post); 
trSDF           = nan(length(-pre:post),length(chans),length(triggerpoints));
for ch = 1:length(chans)
    clear elabel SPK eidx I sdf
    elabel      = sprintf('%s%02u',el,chans(ch));
    eidx        = find(cell2mat(cellfun(@(x) contains(x',elabel),{NEV.ElectrodesInfo.ElectrodeLabel},'UniformOutput',0)));
    
    if isempty(eidx)
        error('no %s',elabel)
    end
    
    I               =  NEV.Data.Spikes.Electrode == eidx;
    Fs              = double(NEV.MetaTags.SampleRes);
    SPK             = double(NEV.Data.Spikes.TimeStamp(I)); % in samples
    sdf             = spk2sdf(SPK,Fs);
    trSDF(:,ch,:)   = trigData(sdf',triggerpoints,pre,post);
    wholesdf(1:length(sdf),ch) = sdf; 
end

NS_Header     = openNSx(strcat(BRdatafile,'.',extension),'noread');
analog        = find(strcmp({NS_Header.ElectrodesInfo.ConnectorBank},'E'));
AnalogLabels  = {NS_Header.ElectrodesInfo(analog).Label};
injection_ttl = analog(contains(AnalogLabels,'ainp8')); 
electrode     = sprintf('c:%u',injection_ttl);
TTL           = openNSx(strcat(BRdatafile,'.',extension),electrode,'read','uV');
TTL           = TTL.Data; 








