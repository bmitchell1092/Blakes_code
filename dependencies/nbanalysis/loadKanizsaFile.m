function [trLFP,trTHETA,trGAMMA,onsets,conds,PRE,POST,Fs] = loadKanizsaFile(BRdatafile,brdrname,extension,el)

%%
PRE = 300; POST = 1000; 
[LFP, EventCodes, EventTimes]= getLFP(strcat(brdrname,'/',BRdatafile),extension,el);

%% low pass filter LFP at 256 Hz
Fs = 1000;
lpc = 256;
nyq =Fs/2;
lWn = lpc/nyq;
[bwb,bwa] = butter(4,lWn,'low');
LFP= filtfilt(bwb,bwa,LFP);

%% filter data in theta and gamma bands
lpc = [3 6];
lWn = lpc/nyq;
[bwb,bwa] = butter(4,lWn,'bandpass');
thetaLFP= filtfilt(bwb,bwa,LFP);

lpc = [25 50];
lWn = lpc/nyq;
[bwb,bwa] = butter(4,lWn,'bandpass');
gammaLFP= filtfilt(bwb,bwa,LFP);
%% divide data into trials
% different event code system:
% 1- fix on
% 2- fix off
% 3 - stimulus off
% 4-fix starts
% 5-no initial fixation
% 6-broke initial fixation
% 7-broke fixation during stim pres
% 8-broke fix between sitm presentations
%10-correct
%11-reward
% condition number = condn - 20 <---coded when stimulus comes on
% [ look for fix on --1//condition number//%fix off (2)//correct
% (10)//reward (11)]
% load MUA, LFP, and BHV

bhvfile = strcat(brdrname,'/',BRdatafile,'.mat'); 
load(bhvfile,'-mat');

% parse event codes
[pEvC, pEvT] = parsEventCodesML(EventCodes,EventTimes);
% check that tr numbers match (DEV: make better)
if length(pEvC) ~= length(BHV.TrialError)
    ecst = EventCodes(1);
    ecen = EventCodes(end);
    if ecst ~= 9
        fprintf('\nfirst trial is bad, correcting...\n')
        BHV.TrialError(1) = [];
        BHV.ConditionNumber(1) = [];
    elseif ecen ~= 18
        fprintf('\nlast trial is bad, correcting...\n')
        BHV.TrialError(end) = [];
        BHV.ConditionNumber(end) = [];
    elseif strcmp(datafile,'140806_B_qkanizsa_inducer001')
        % BR was turned off before ML, manually checked 1/28/2014
        ntr = length(pEvC);
        BHV.TrialError(length(pEvC)+1:end) = [];
        BHV.ConditionNumber(length(pEvC)+1:end) = [];
    else
        fprintf('\ntrls do not match, skipping file...\n')
        sct = sct + 1;
        skipfile{sct} = datafile;
        skipmessage{sct} = 'issue syncing behavoral events';
    end
end

% determin conditions of interest
allcondnum = unique(BHV.ConditionNumber(BHV.TrialError == 0));
stimnames = {'C1' 'C2' 'C3' 'C4' 'C5' 'C6' 'C7'};
stimcolors = [ 1 0 0; 0 0 1];
conditonnumbers = []; cct = 0;
for c = 1:length(allcondnum)
    condition = allcondnum(c);
    condstr = BHV.TaskObject{condition,2};
    for s = 1:length(stimnames)
        if ~isempty(strfind(condstr,stimnames{s}))
            cct = cct+1;
            conditonnumbers(cct) = condition;
        end
    end
end


% get triggerpoints for each condition
conds = [];
onsets = []; 
for c = 1:length(conditonnumbers)
    condition = conditonnumbers(c);
    ct = 0; condstimon=[];
    maxtr = length(BHV.TrialError);
    for tr = 1:maxtr;
        if BHV.TrialError(tr) == 0 ...
                && BHV.ConditionNumber(tr) == condition;
            ct = ct + 1;
            
            tr_codes = pEvC{tr};
            tr_times = pEvT{tr};
            
            condstimon(ct) = tr_times(tr_codes == 20 + condition); % Schmid Lab Event Code Scheme, see timing file
        end
    end
    onsets = [onsets condstimon]; 
    conds  = [conds repmat(c,length(condstimon),1)']; 
end
%%
%% cut into trials : 

for tr = 1:length(onsets)
    
    refwin        =  floor([onsets(tr)./30 + PRE : onsets(tr)./30 + POST]);
    trLFP(:,:,tr) = LFP(refwin,:);
    trTHETA(:,:,tr) = thetaLFP(refwin,:);
    trGAMMA(:,:,tr) = gammaLFP(refwin,:);
    
end