%% BMcinteroc
% script to: select an input file, load in .nev,
% pull out stim onset, load in .ns2 and .ns6 and sync event codes to the
% raw neural data. Each stim onset with animal fixation for the duration of presentation is a trial. 
% Generate LFP, aMUA, and CSD --triggered to a reference
% window. Average across trials and baseline corrected when necessary. 
% Several plotting options for data visualization at the end
clear
cd('C:\Users\bmitc\')

%% Establish directories and set path

if strcmp(getenv('USER'),'maierav')                                      %retrieves environment variable 'USER' 
    npmkdir  = '/Users/alex 1/Desktop/LAB/Brock/OLD/NPMK-4.5.3.0/NPMK/'; %directory for Alex's machine
    nbanalysisdir   = '/Users/alex 1/Desktop/LAB/bootcamp/nbanalysis/';  %directory for Alex's machine
    datadir  = '/Users/alex 1/Desktop/LAB/';                             %directory for the stored data
else
    npmkdir  = '/users/bmitc/Documents/MATLAB/NPMK/';                    %neural processing matlab kit (NPMK)
    nbanalysisdir   = '/users/bmitc/Documents/MATLAB/nbanalysis/';       %directory with various tools for opening, loading, and processing 
    datadir  = '/users/bmitc/Box Sync/DATA/';                            %this is my Vanderbilt Box sync 
    %datadir = 'users/bmitc/Documents/MATLAB/data/';
end

addpath(genpath(npmkdir))
addpath(genpath(nbanalysisdir))
addpath(genpath(datadir))

BRdatafile = '160505_E_mcosinteroc002'; 
filename = [datadir BRdatafile];

%% Define layers by contact (informed by datalogs)

[supra,gran,infra] = layerLog(BRdatafile);

%% Define stimulus patterns

patterns   = {'rforidrft','rfsfdrft','posdisparitydrft','disparitydrft','cinterocdrft','coneinterocdrft','conedrft', ...
                'colorflicker','bwflicker','rfori','rfsize','cinteroc','color','rfsf','mcosinteroc','dotmapping'}; 
for p = 1:length(patterns) 
    
    pattern = patterns{p};    %pattern is individual elements of the above array 'patterns'
    
    if any(strfind(BRdatafile,pattern))         %if BRdatafile contains any string the same as pattern
        startlog = strfind(BRdatafile,pattern); 
        if ~isequal(BRdatafile(startlog:end-3),pattern), continue
        else
        match = patterns{p};        % finding which file type the BRdatafile is
        end
    end
    
end

if isequal(match,'dotmapping')
    ext = '.gDotsXY_di';
elseif isequal(BRdatafile, '190415_B_cinteroc002') || isequal(BRdatafile, '190321_B_cinteroc001')...  % specific files that are named cinteric but actually cinterocDRFT
        || isequal(BRdatafile,'161003_E_cinteroc002') || isequal(BRdatafile,'190210_B_cinteroc001')
    ext = ['.g' upper(match) 'DRFTGrating_di']; % defining an extension for DRFT files
else
    ext = ['.g' upper(match) 'Grating_di']; % defining an extension for grating files
end

if contains(ext,'DRFT')                     
      grating     = readgDRFTGrating([filename ext]); % from nbanalysis 
elseif contains(ext,'Dots')
      grating     = readgDotsXY([filename ext]);
else
      grating     = readgGrating([filename ext]);
end

%% Load event times and codes

NEV             = openNEV([filename '.nev'],'noread','overwrite');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;          
EventSamples    = NEV.Data.SerialDigitalIO.TimeStamp;                   %Events in samples 
EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec.*1000);   %floor rounds to nearest integer and then convert event to ms 
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSamples);            %sorts codes, samps or times into trials

%The following is a structure that only contains trials where the animal
%did not break fixation, and includes grating info and stimulus onsets. 

STIM            = sortStimandTimeData(grating,pEvC,pEvT,'stim'); 
STIM.onsetsdown = floor(STIM.onsets./30);
STIM.laminae.supra = supra;
STIM.laminae.gran = gran;
STIM.laminae.infra = infra;

%% Load LFP with NS2 file
clear ext
ext = 'ns2';

%retrieve sort direction and bank info
NS_Header    = openNSx(strcat(filename,'.',ext),'noread');
banks        = unique({NS_Header.ElectrodesInfo.ConnectorBank}); banks(ismember(banks,'E')) = []; % bank E is BNC cable inputs

for b = 1:length(banks)
    clear neural label 
    neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},banks{b}); 
    firstlabel   = cell2mat({NS_Header.ElectrodesInfo(find(neural,1,'first')).Label}); 
    if str2double(firstlabel(3:4)) < 2
        sortdirection = 'ascending'; 
    else
        sortdirection = 'descending'; 
    end
end

if any(strfind(firstlabel,'eD'))
    ebank = 'eD'; %electrode bank
elseif any(strfind(firstlabel, 'eC'))
    ebank = 'eC';
elseif any(strfind(firstlabel, 'eB')) 
    ebank = 'eB';
else 
    error('Could not find ebank')
end

% load LFP with NS2 file

lfp = getLFP(filename,ext,ebank,sortdirection);


%% LOAD LFP and analog MUA with NS6 file.

clear ext NS_header banks neural 
% Read in NS Header
ext          = 'ns6'; 
NS_Header    = openNSx(strcat(filename,'.',ext),'noread');

% get basic info about recorded data
neural       = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},ebank(2)); % logicals where contact bank name matches electrode of interest
N.neural     = sum(neural); % number of neural channels 
NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label}; %get labels
Fs           = NS_Header.MetaTags.SamplingFreq; % get sampling frequency
nyq          = Fs/2; 
r            = Fs/1000; 

% counters
clear nct
nct = 0;

% process data electrode by electrode
for e = 1:length(neural)
    
    if neural(e) == 1    % neural is a vector of logicals, 1 = contacts we want

        nct = nct+1;
        
        % open data for this channel. 
        clear NS DAT
        electrode = sprintf('c:%u',e);
        NS        = openNSx(strcat(filename,'.',ext),electrode,'read','uV');
        DAT       = NS.Data; NS.Data = [];  % this is the whole signal on one channel, 30 kHz!
        
        
        % preallocate data matrices 
        if nct == 1
            N.samples = length(DAT); 
            LFP       = nan(ceil(N.samples/r),N.neural); % preallocating for downsampled data
            MUA       = nan(ceil(N.samples/r),N.neural);
        end
        
        % extract the LFP. 
        clear lpc lWn bwb bwa lpLFP
        lpc       = 200; %low pass cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        lpLFP      = filtfilt(bwb,bwa,DAT);  %low pass filter 
        
        % extract the MUA:
        clear hpc hWn bwb bwa hpMUA
        hpc       = 750;  %high pass cutoff
        hWn       = hpc/nyq;
        [bwb,bwa] = butter(4,hWn,'high');
        hpMUA     = filtfilt(bwb,bwa,DAT); %high pass filter
        
        % low pass at 5000 Hz and rectify (take the absolute value) 
        clear lpc lWn bwb bwa 
        lpc       = 5000;  % cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        hpMUA     = abs(filtfilt(bwb,bwa,hpMUA)); %low pass filter &rectify
        
        % low pass filter at x Hz. 
        clear lpc lWn bwb bwa lpMUA
        lpc       = 200; %low pass cutoff
        lWn       = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low'); 
        lpMUA     = filtfilt(bwb,bwa,hpMUA);  %low pass filter to smooth
        
      
        % decimate both LFP and analog MUA (aMUA) to get 1kHz samp freq
        MUA(:,nct) = decimate(lpMUA,r); 
        LFP(:,nct) = decimate(lpLFP,r); 
        
        clear DAT 
        
    end
    
end 

% sort electrode contacts in ascending order:
idx    = nan(length(NeuralLabels),1); 
for ch = 1:length(NeuralLabels)
    chname  = strcat(sprintf('%s',ebank),sprintf('%02d',ch));
    idx(ch) = find(contains(NeuralLabels,chname));
    
end

% not 100% positive that this is working as intended. Though the datalog
% for layers does appear to match the resulting CSD 
switch sortdirection
    case 'ascending'
        MUA = MUA(:,idx);
        LFP = LFP(:,idx);
        sortedLabels = NeuralLabels(idx); 
    case 'descending'
        MUA = MUA(:,flipud(idx)); %kacie used fliplr here but it wasn't flipping properly 
        LFP = LFP(:,flipud(idx)); % same as above
        sortedLabels = NeuralLabels(flipud(idx)); % same as above
end
%% calculate CSD 
% calculate CSD before triggering to trials OR on the trial data BUT not on
% the mean LFP.

CSD = mod_iCSD(LFP')';  % this function takes LFP in channels x samples so let's transpose LFP and then flip it right back 
                        % feed in units of microV and get back units of
                        % nA/mm^3


CSD = padarray(CSD,[0 1],NaN,'replicate'); % pad array if you want to keep the matrix the same size on the channel
                                           % dimension as the other matrices

%% trigger the neural data to the event codes of interest                                         

offset = round(((STIM.offsets(1)-STIM.onsets(1))/(30)),0); % stimulus offset as calculated from grating text file.
pre   = -50; % 50ms before stim onset
post  = (round(10^-2*offset)/10^-2)+100; % ~100 ms after stim offset

STIM.LFP  = trigData(LFP,STIM.onsetsdown,-pre,post); %pre variable is in absolute units 
STIM.CSD  = trigData(CSD,STIM.onsetsdown,-pre,post); 
STIM.aMUA = trigData(MUA,STIM.onsetsdown,-pre,post); 

%% Averaging across trials & baseline correct
clear avg
STIM.avg.LFP = mean(STIM.LFP,3);
STIM.avg.aMUA = mean(STIM.aMUA,3);
STIM.avg.CSD = mean(STIM.CSD,3);

% These aren't currently used because I'm about to convert to percent
% change from baseline (next section). 

[STIM.bsl.LFP] = BMbasecorrect(STIM.avg.LFP); 
[STIM.bsl.aMUA] = BMbasecorrect(STIM.avg.aMUA); 
[STIM.bsl.CSD] = BMbasecorrect(STIM.avg.CSD);

%% Data conversion 
% Converting the aMUA raw neural response to either percent change from
% baseline or z score (baseline as Xbar). 

% pre-allocate
clear cMUA
STIM.cMUA = nan(size(STIM.aMUA,1),size(STIM.aMUA,2),size(STIM.aMUA,3));
STIM.zMUA = nan(size(STIM.aMUA,1),size(STIM.aMUA,2),size(STIM.aMUA,3));
% aMUA conversion to either z-score or percent change
clear t c
for t = 1:size(STIM.aMUA,3)
    for c = 1:size(STIM.aMUA,2)
        STIM.zMUA(:,c,t) = (STIM.aMUA(:,c,t)-mean(STIM.aMUA(25:75,c,t)))./(std(STIM.aMUA(25:75,c,t))); %z score
%         STIM.cMUA(:,c,t) = (STIM.aMUA(:,c,t)-mean(STIM.aMUA(25:75,c,t)))./(mean(STIM.aMUA(25:75,c,t)))*100; %percent change
        STIM.cMUA(:,c,t) = (STIM.aMUA(:,c,t)-mean(STIM.aMUA(25:75,c,t))); % baseline corrected raw
    end
end

%% Defining conditions of interest
% This is to create structures for the different conditions (including
% contrast levels and eye to which stimulus was shown. 

STIM.levels = unique(STIM.contrast);  % variable that contains all contrast levels
clear STIM.BINconditions STIM.DEconditions STIM.NDEconditions
clear i STIM.conditions.DE
for i = 1:length(STIM.levels)  % monocular (DE) contrast conditions
STIM.conditions.DE(i,:) = STIM.contrast == STIM.levels(i) & STIM.fixedc == 0; 
end

clear i STIM.conditions.NDE

for i = 1:length(STIM.levels)  % monocular (NDE) contrast conditions
STIM.conditions.NDE(i,:) = STIM.contrast == 0 & STIM.fixedc == STIM.levels(i); 
end

clear i STIM.conditions.BIN
for i = 1:length(STIM.levels)  % binocular contrast (same contrast in two eyes)
STIM.conditions.BIN(i,:) = STIM.contrast == STIM.levels(i) & STIM.fixedc == STIM.levels(i); 
end

for i = 1:length(STIM.levels)  % binocular contrast (same contrast in two eyes)
STIM.conditions.DI_DElow(i,:) = STIM.contrast == STIM.levels(2) & STIM.fixedc == STIM.levels(i); 
end

for i = 1:length(STIM.levels)  % binocular contrast (same contrast in two eyes)
STIM.conditions.DI_DEmed(i,:) = STIM.contrast == STIM.levels(3) & STIM.fixedc == STIM.levels(i); 
end

for i = 1:length(STIM.levels)  % binocular contrast (same contrast in two eyes)
STIM.conditions.DI_DEhigh(i,:) = STIM.contrast == STIM.levels(4) & STIM.fixedc == STIM.levels(i); 
end

clear STIM.conditions.DI 

STIM.conditions.DI_exclusive(1,:) = STIM.contrast == STIM.levels(2) & STIM.fixedc == STIM.levels(3); % DE22NDE45
STIM.conditions.DI_exclusive(2,:) = STIM.contrast == STIM.levels(2) & STIM.fixedc == STIM.levels(4); % DE22NDE90
STIM.conditions.DI_exclusive(3,:) = STIM.contrast == STIM.levels(3) & STIM.fixedc == STIM.levels(2); % DE45NDE22
STIM.conditions.DI_exclusive(4,:) = STIM.contrast == STIM.levels(3) & STIM.fixedc == STIM.levels(4); % DE45NDE90
STIM.conditions.DI_exclusive(5,:) = STIM.contrast == STIM.levels(4) & STIM.fixedc == STIM.levels(2); % DE90NDE22
STIM.conditions.DI_exclusive(6,:) = STIM.contrast == STIM.levels(4) & STIM.fixedc == STIM.levels(3); % DE90NDE45

%% Averaged trials by condition
% The above matrices were just logicals. Now time to take my full cMUA
% trials and use those logicals to trigger to and average by conditions of
% interest.

% aMUA by condition
clear m STIM.DE.cMUA STIM.NDE.cMUA STIM.BIN.cMUA 
for m = 1:size(STIM.conditions.DE,1)
    STIM.DE.cMUA(m).contrast = mean(STIM.cMUA(:,:,STIM.conditions.DE(m,:)),3); 
    STIM.NDE.cMUA(m).contrast = mean(STIM.cMUA(:,:,STIM.conditions.NDE(m,:)),3); 
    STIM.BIN.cMUA(m).contrast = mean(STIM.cMUA(:,:,STIM.conditions.BIN(m,:)),3);
    STIM.DI.DElow.cMUA(m).contrast = mean(STIM.cMUA(:,:,STIM.conditions.DI_DElow(m,:)),3);
    STIM.DI.DEmed.cMUA(m).contrast = mean(STIM.cMUA(:,:,STIM.conditions.DI_DEmed(m,:)),3);
    STIM.DI.DEhigh.cMUA(m).contrast = mean(STIM.cMUA(:,:,STIM.conditions.DI_DEhigh(m,:)),3);
end

clear m
for m = 1:size(STIM.conditions.DI_exclusive,1)
    STIM.DI_exclusive.cMUA(m).contrast = mean(STIM.cMUA(:,:,STIM.conditions.DI_exclusive(m,:)),3);
end

% CSD by condition
clear m Mon_CSD BIN_CSD
for m = 1:size(STIM.conditions.DE,1)
    STIM.DE.CSD(m).contrast = mean(STIM.CSD(:,:,STIM.conditions.DE(m,:)),3); 
    STIM.NDE.CSD(m).contrast = mean(STIM.CSD(:,:,STIM.conditions.NDE(m,:)),3);
    STIM.BIN.CSD(m).contrast = mean(STIM.CSD(:,:,STIM.conditions.BIN(m,:)),3);
    STIM.DI.DElow.CSD(m).contrast = mean(STIM.CSD(:,:,STIM.conditions.DI_DElow(m,:)),3);
    STIM.DI.DEmed.CSD(m).contrast = mean(STIM.CSD(:,:,STIM.conditions.DI_DEmed(m,:)),3);
    STIM.DI.DEhigh.CSD(m).contrast = mean(STIM.CSD(:,:,STIM.conditions.DI_DEhigh(m,:)),3);
end

for m = 1:size(STIM.conditions.DI_exclusive,1)
    STIM.DI.exclusive.CSD(m).contrast = mean(STIM.CSD(:,:,STIM.conditions.DI_exclusive(m,:)),3);
end

%% collapsing across time for each condition
% In order to get a single number to represent the contrast response for each contact, I can collapse
% across time. The following collapses across different lengths of time:
% full stim duration (80ms to offset), transient (80 to 150ms), and
% sustained (151ms to offset). 

clear i 
clear STIM.DE.coll STIM.NDE.coll STIM.BIN.coll
for i=1:size(STIM.DE.cMUA,2)
    STIM.DE.coll.full(i,:)  = mean(STIM.DE.cMUA(i).contrast(80:offset,:),1);
    STIM.NDE.coll.full(i,:)  = mean(STIM.NDE.cMUA(i).contrast(80:offset,:),1);
    STIM.BIN.coll.full(i,:)  = mean(STIM.BIN.cMUA(i).contrast(80:offset,:),1);
    STIM.DI.DElow.coll.full(i,:)  = mean(STIM.DI.DElow.cMUA(i).contrast(80:offset,:),1);
    STIM.DI.DEmed.coll.full(i,:)  = mean(STIM.DI.DEmed.cMUA(i).contrast(80:offset,:),1);
    STIM.DI.DEhigh.coll.full(i,:)  = mean(STIM.DI.DEhigh.cMUA(i).contrast(80:offset,:),1);
    STIM.DE.coll.transient(i,:)  = mean(STIM.DE.cMUA(i).contrast(80:150,:),1);
    STIM.NDE.coll.transient(i,:)  = mean(STIM.NDE.cMUA(i).contrast(80:150,:),1);
    STIM.BIN.coll.transient(i,:)  = mean(STIM.BIN.cMUA(i).contrast(80:150,:),1);
    STIM.DI.DElow.coll.transient(i,:)  = mean(STIM.DI.DElow.cMUA(i).contrast(80:150,:),1);
    STIM.DI.DEmed.coll.transient(i,:)  = mean(STIM.DI.DEmed.cMUA(i).contrast(80:150,:),1);
    STIM.DI.DEhigh.coll.transient(i,:)  = mean(STIM.DI.DEhigh.cMUA(i).contrast(80:150,:),1);
    STIM.DE.coll.sustained(i,:)  = mean(STIM.DE.cMUA(i).contrast(151:offset,:),1);
    STIM.NDE.coll.sustained(i,:)  = mean(STIM.NDE.cMUA(i).contrast(151:offset,:),1);
    STIM.BIN.coll.sustained(i,:)  = mean(STIM.BIN.cMUA(i).contrast(151:offset,:),1);
    STIM.DI.DElow.coll.sustained(i,:)  = mean(STIM.DI.DElow.cMUA(i).contrast(151:offset,:),1);
    STIM.DI.DEmed.coll.sustained(i,:)  = mean(STIM.DI.DEmed.cMUA(i).contrast(151:offset,:),1);
    STIM.DI.DEhigh.coll.sustained(i,:)  = mean(STIM.DI.DEhigh.cMUA(i).contrast(151:offset,:),1);
end

clear i
for i=1:size(STIM.DI.exclusive.cMUA,2)
    STIM.DI.exclusive.coll.full(i,:)  = mean(STIM.DI.exclusive.cMUA(i).contrast(80:offset,:),1);
    STIM.DI.exclusive.coll.transient(i,:)  = mean(STIM.DI.exclusive.cMUA(i).contrast(80:150,:),1);
    STIM.DI.exclusive.coll.sustained(i,:)  = mean(STIM.DI.exclusive.cMUA(i).contrast(151:offset,:),1);
end

%% Binning contacts into V1 layers (informed by datalogs)
% Seperated into full stim duration, transient response, and sustained
% response

STIM.DE.layers.full = [mean(STIM.DE.coll.full(:,STIM.laminae.supra),2),mean(STIM.DE.coll.full(:,STIM.laminae.gran),2),mean(STIM.DE.coll.full(:,STIM.laminae.infra),2)];
STIM.DE.layers.transient = [mean(STIM.DE.coll.transient(:,STIM.laminae.supra),2),mean(STIM.DE.coll.transient(:,STIM.laminae.gran),2),mean(STIM.DE.coll.transient(:,STIM.laminae.infra),2)];
STIM.DE.layers.sustained = [mean(STIM.DE.coll.sustained(:,STIM.laminae.supra),2),mean(STIM.DE.coll.sustained(:,STIM.laminae.gran),2),mean(STIM.DE.coll.sustained(:,STIM.laminae.infra),2)];

STIM.NDE.layers.full = [mean(STIM.NDE.coll.full(:,STIM.laminae.supra),2),mean(STIM.NDE.coll.full(:,STIM.laminae.gran),2),mean(STIM.NDE.coll.full(:,STIM.laminae.infra),2)];
STIM.NDE.layers.transient = [mean(STIM.NDE.coll.transient(:,STIM.laminae.supra),2),mean(STIM.NDE.coll.transient(:,STIM.laminae.gran),2),mean(STIM.NDE.coll.transient(:,STIM.laminae.infra),2)];
STIM.NDE.layers.sustained = [mean(STIM.NDE.coll.sustained(:,STIM.laminae.supra),2),mean(STIM.NDE.coll.sustained(:,STIM.laminae.gran),2),mean(STIM.NDE.coll.sustained(:,STIM.laminae.infra),2)];

STIM.BIN.layers.full = [mean(STIM.BIN.coll.full(:,STIM.laminae.supra),2),mean(STIM.BIN.coll.full(:,STIM.laminae.gran),2), mean(STIM.BIN.coll.full(:,STIM.laminae.infra),2)];
STIM.BIN.layers.transient = [mean(STIM.BIN.coll.transient(:,STIM.laminae.supra),2),mean(STIM.BIN.coll.transient(:,STIM.laminae.gran),2),mean(STIM.BIN.coll.transient(:,STIM.laminae.infra),2)];
STIM.BIN.layers.sustained = [mean(STIM.BIN.coll.sustained(:,STIM.laminae.supra),2),mean(STIM.BIN.coll.sustained(:,STIM.laminae.gran),2),mean(STIM.BIN.coll.sustained(:,STIM.laminae.infra),2)];

STIM.DI.DElow.layers.full = [mean(STIM.DI.DElow.coll.full(:,STIM.laminae.supra),2),mean(STIM.DI.DElow.coll.full(:,STIM.laminae.gran),2), mean(STIM.DI.DElow.coll.full(:,STIM.laminae.infra),2)];
STIM.DI.DElow.layers.transient = [mean(STIM.DI.DElow.coll.transient(:,STIM.laminae.supra),2),mean(STIM.DI.DElow.coll.transient(:,STIM.laminae.gran),2),mean(STIM.DI.DElow.coll.transient(:,STIM.laminae.infra),2)];
STIM.DI.DElow.layers.sustained = [mean(STIM.DI.DElow.coll.sustained(:,STIM.laminae.supra),2),mean(STIM.DI.DElow.coll.sustained(:,STIM.laminae.gran),2),mean(STIM.DI.DElow.coll.sustained(:,STIM.laminae.infra),2)];

STIM.DI.DEmed.layers.full = [mean(STIM.DI.DEmed.coll.full(:,STIM.laminae.supra),2),mean(STIM.DI.DEmed.coll.full(:,STIM.laminae.gran),2), mean(STIM.DI.DEmed.coll.full(:,STIM.laminae.infra),2)];
STIM.DI.DEmed.layers.transient = [mean(STIM.DI.DEmed.coll.transient(:,STIM.laminae.supra),2),mean(STIM.DI.DEmed.coll.transient(:,STIM.laminae.gran),2),mean(STIM.DI.DEmed.coll.transient(:,STIM.laminae.infra),2)];
STIM.DI.DEmed.layers.sustained = [mean(STIM.DI.DEmed.coll.sustained(:,STIM.laminae.supra),2),mean(STIM.DI.DEmed.coll.sustained(:,STIM.laminae.gran),2),mean(STIM.DI.DEmed.coll.sustained(:,STIM.laminae.infra),2)];

STIM.DI.DEhigh.layers.full = [mean(STIM.DI.DEhigh.coll.full(:,STIM.laminae.supra),2),mean(STIM.DI.DEhigh.coll.full(:,STIM.laminae.gran),2), mean(STIM.DI.DEhigh.coll.full(:,STIM.laminae.infra),2)];
STIM.DI.DEhigh.layers.transient = [mean(STIM.DI.DEhigh.coll.transient(:,STIM.laminae.supra),2),mean(STIM.DI.DEhigh.coll.transient(:,STIM.laminae.gran),2),mean(STIM.DI.DEhigh.coll.transient(:,STIM.laminae.infra),2)];
STIM.DI.DEhigh.layers.sustained = [mean(STIM.DI.DEhigh.coll.sustained(:,STIM.laminae.supra),2),mean(STIM.DI.DEhigh.coll.sustained(:,STIM.laminae.gran),2),mean(STIM.DI.DEhigh.coll.sustained(:,STIM.laminae.infra),2)];

STIM.DI.exclusive.layers.full = [mean(STIM.DI.exclusive.coll.full(:,STIM.laminae.supra),2),mean(STIM.DI.exclusive.coll.full(:,STIM.laminae.gran),2), mean(STIM.DI.exclusive.coll.full(:,STIM.laminae.infra),2)];
STIM.DI.exclusive.layers.transient = [mean(STIM.DI.exclusive.coll.transient(:,STIM.laminae.supra),2),mean(STIM.DI.exclusive.coll.transient(:,STIM.laminae.gran),2),mean(STIM.DI.exclusive.coll.transient(:,STIM.laminae.infra),2)];
STIM.DI.exclusive.layers.sustained = [mean(STIM.DI.exclusive.coll.sustained(:,STIM.laminae.supra),2),mean(STIM.DI.exclusive.coll.sustained(:,STIM.laminae.gran),2),mean(STIM.DI.exclusive.coll.sustained(:,STIM.laminae.infra),2)];

%% Plotting: (SNAPSHOT)
% All averaged, baseline corrected trials (LFP, aMUA, CSD, iCSD)

STIM.refwin = pre:post; % reference window for line plotting
STIM.channels = 1:nct;  % how many channels (nct is a predefined variable with the exact number of channels

h1 = figure('position',[15,135,1200,500]);
clear i
avg_fields = fieldnames(STIM.avg);
for i = 1:length(avg_fields)
subplot(1,4,i)
f_ShadedLinePlotbyDepthMod((STIM.avg.(avg_fields{i})),0:(1/(numel(STIM.channels))):1,STIM.refwin, STIM.channels, 1); % this function baseline corrects and scales
hold on
plot([0 0], ylim,'k')
plot([offset offset], ylim,'k','linestyle','-.','linewidth',0.5)
title(avg_fields{i})
xlabel('time (ms)')
ylabel('contacts indexed down from surface')
hold off
end

bAVG_iCSD = filterCSD(STIM.bsl.CSD')';

h = subplot(1,4,4);
imagesc(STIM.refwin,STIM.channels,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(0); set(v,'color','k','linestyle','-','linewidth',1);
set(gca,'tickdir','out');  
climit = max(abs(get(gca,'CLim'))*1);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
plot([offset offset], ylim,'k','linestyle','-.','linewidth',0.5)
title('Interpolated CSD')
xlabel('time (ms)')
clrbar = colorbar; clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts indexed down from surface');
hold off

sgtitle({'All trials triggered to stim onset',BRdatafile}, 'Interpreter', 'none');
set(h,'position',[0.7483,0.1253,0.1055,0.6826]);

% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_snapshot',BRdatafile), '-jpg', '-transparent');

%% Dominant Eye Contrast responses (Contrasts)
% Dominant eye response to stimulus with different contrast levels

figure('position',[15,135,1200,500]);
subplot(1,length(STIM.levels),numel(STIM.levels))
hold on
contrastValue = max(STIM.levels);
global scalingfactor
f_ShadedLinePlotbyDepth(mean(STIM.aMUA(:,:,STIM.conditions.DE(numel(STIM.levels),:)),3),0:(1/(numel(STIM.channels))):1,STIM.refwin,STIM.channels,1,1);
title('1 contrast in DE');
xlabel('time (ms)');
hold off

clear i 
for i = 1:length(STIM.levels)-1
    subplot(1,length(STIM.levels),i);
    f_ShadedLinePlotbyDepth_BAM(mean(STIM.aMUA(:,:,STIM.conditions.DE(i,:)),3),0:(1/(numel(STIM.channels))):1,STIM.refwin,STIM.channels,1,1,false,scalingfactor);
   
plot([0 0], ylim,'k')
plot([offset offset], ylim,'k')
if i == 1
    title({STIM.levels(i),' STIM.levels in both eyes'});
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
    bar(STIM.BIN.coll.transient(:,selectchannels(c)),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
    hold on
    bar(STIM.DE.coll.transient(:,selectchannels(c)),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
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

sgtitle({'Monocular vs Binocular contrast response function (MUA)'...
    'Responses collapsed across initial 100ms',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_bar-contrasts-transient',BRdatafile), '-jpg', '-transparent');

%% All contacts, collapsed across time, sectioned by contrast level (Tightplot)
figure('position',[185,41.666666666666664,645.3333333333333,599.3333333333333]);

[ha, pos] = tight_subplot(2,3,[0.005 .03],[.10 .2],[.05 .05]); %channels, columns, [spacing], [bottom and top margin], [left and right margin]
clear c
for c = 1:3
    
    axes(ha(c)); % ha is a variable that gets the axis of each subplot
    plot(fliplr(STIM.DE.coll.transient(c+1,:)),STIM.channels,'b','linewidth',.5);
    hold on
    plot(fliplr(STIM.BIN.coll.transient(c+1,:)),STIM.channels,'.-r','linewidth',0.5);
    hline(STIM.channels(end)-STIM.laminae.gran(end),'-.')
    xlim([-2 10])
    yticks('')
    yticklabels(fliplr(1:nct))
    ylim([1 nct])
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
    plot(fliplr(STIM.DE.coll.sustained(c-2,:)),STIM.channels,'b','linewidth',.5);
    hold on
    plot(fliplr(STIM.BIN.coll.sustained(c-2,:)),STIM.channels,'.-r','linewidth',0.5);
    hline(STIM.channels(end)-STIM.laminae.gran(end),'-.')
    xlim([-10 80])
    yticks('')
    yticklabels(fliplr(1:nct))
    ylim([1 nct])
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
cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_tightplot',BRdatafile), '-jpg', '-transparent');

%% Collapsed contrast responses for monocular and binocular stim (contrast lines)
% Designed to show Monocular vs binocular contrast responses alongside CSD
% for layer reference. 
figure('position',[185 150 887 450]);

h = subplot(1,4,1);
bAVG_iCSD = filterCSD(STIM.bsl.CSD')';
imagesc(STIM.refwin,STIM.channels,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(0); set(v,'color','k','linestyle','-','linewidth',1);
set(gca,'tickdir','out');  
climit = max(abs(get(gca,'CLim'))*1);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
plot([offset offset], ylim,'k','linestyle','-.','linewidth',0.5)
title('iCSD: All Trials')
xlabel('time (ms)')
clrbar = colorbar; %clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts indexed down from surface');
set(h,'position',[0.065388951521984,0.097526988745119,0.145749605022835,0.722586325702473]);
hold off

subplot(1,4,3)
plot(fliplr(STIM.BIN.coll.transient(:,:)),STIM.channels);
hold on
%ylabel('contacts indexed down from surface');
set(gca,'box','off');
%yticklabels({'','20','15','10','5'});
xlabel('Percent change');
grid on
climit = max(abs(get(gca,'xlim'))*1);
xlim([-5 climit*1.2]);
yticks(1:nct)
yticklabels(fliplr(1:nct))
ylim([1 nct])
title('Binocular');
hold off

subplot(1,4,2)
hold on
plot(fliplr(STIM.DE.coll.transient(:,:)),STIM.channels);
set(gca,'box','off');
grid on
xlim([-5 climit*1.2]);
yticks(1:nct)
yticklabels(fliplr(1:nct))
ylim([1 nct])
xlabel('Percent change');
grid on
title('Monocular');
hold off

STIM.calc.contacts.subtractionDE.full = STIM.BIN.coll.full(:,:)-STIM.DE.coll.full(:,:);
STIM.calc.contacts.subtractionDE.transient = STIM.BIN.coll.transient(:,:)-STIM.DE.coll.transient(:,:);
STIM.calc.contacts.subtractionDE.sustained = STIM.BIN.coll.sustained(:,:)-STIM.DE.coll.sustained(:,:);

STIM.calc.contacts.subtractionNDE.full = STIM.BIN.coll.full(:,:)-STIM.NDE.coll.full(:,:);
STIM.calc.contacts.subtractionNDE.transient = STIM.BIN.coll.transient(:,:)-STIM.NDE.coll.transient(:,:);
STIM.calc.contacts.subtractionNDE.sustained = STIM.BIN.coll.sustained(:,:)-STIM.NDE.coll.sustained(:,:);

subplot(1,4,4)
plot(fliplr(STIM.calc.contacts.subtractionDE.transient),STIM.channels);
hold on 
grid on
xlim([-5 25]);
set(gca,'box','off');
yticks(1:nct)
yticklabels(fliplr(1:nct))
ylim([1 nct])
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
subplot(2,3,1)
clear c
for c = 1:length(STIM.levels)
plot(STIM.refwin,mean(STIM.DE.cMUA(c).contrast(:,STIM.laminae.supra),2),'color','b')
hold on
plot(STIM.refwin,mean(STIM.BIN.cMUA(c).contrast(:,STIM.laminae.supra),2),'color','r')
%ylimit = max(abs(get(gcf,'ylim')));
ylimit = 100;
set(gca,'ylim',[-10 ylimit],'Box','off','TickDir','out')
end
ylabel({'Percent change'...
    'from baseline'});
xlabel('time (ms)');
title('Supragranular');

subplot(2,3,2)
clear c
for c = 1:length(STIM.levels)
plot(STIM.refwin,mean(STIM.DE.cMUA(c).contrast(:,STIM.laminae.gran),2),'color','b')
hold on
plot(STIM.refwin,mean(STIM.BIN.cMUA(c).contrast(:,STIM.laminae.gran),2),'color','r')
set(gca,'ylim',[-10 ylimit],'Box','off','TickDir','out')
end

xlabel('time (ms)');
title('Granular');

subplot(2,3,3)
clear c
for c = 1:length(STIM.levels)
plot(STIM.refwin,mean(STIM.DE.cMUA(c).contrast(:,STIM.laminae.infra),2),'color','b')
hold on
plot(STIM.refwin,mean(STIM.BIN.cMUA(c).contrast(:,STIM.laminae.infra),2),'color','r')
set(gca,'ylim',[-10 ylimit],'Box','off','TickDir','out')
end

xlabel('time (ms)');
title('Infragranular');

subplot(2,3,4)
clear c
for c = 1:length(STIM.levels)
plot(STIM.refwin,smooth(mean(STIM.BIN.cMUA(c).contrast(:,STIM.laminae.supra),2)-(mean(STIM.DE.cMUA(c).contrast(:,STIM.laminae.supra),2)),.1))
hold on
ylimit = max(abs(get(gcf,'ylim')));
set(gca,'ylim',[-10 ylimit/5],'Box','off','TickDir','out')
end
ylabel({'Percent difference'...
    '(bin - mon)'});
xlabel('time (ms)');

subplot(2,3,5)
clear c
for c = 1:length(STIM.levels)
plot(STIM.refwin,smooth(mean(STIM.BIN.cMUA(c).contrast(:,STIM.laminae.gran),2)-(mean(STIM.DE.cMUA(c).contrast(:,STIM.laminae.gran),2)),.1))
hold on
ylimit = max(abs(get(gcf,'ylim')));
set(gca,'ylim',[-10 ylimit/5],'Box','off','TickDir','out')
end

xlabel('time (ms)');

subplot(2,3,6)
clear c
for c = 1:length(STIM.levels)
plot(STIM.refwin,smooth(mean(STIM.BIN.cMUA(c).contrast(:,STIM.laminae.infra),2)-(mean(STIM.DE.cMUA(c).contrast(:,STIM.laminae.infra),2)),.1))
hold on
ylimit = max(abs(get(gcf,'ylim')));
set(gca,'ylim',[-10 ylimit/5],'Box','off','TickDir','out')
end

xlabel('time (ms)');

sgtitle({'V1 laminar contrast response profiles: monocular (blue) vs binocular (red)'...
   ,BRdatafile},'Interpreter','none');

% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_timeplots',BRdatafile), '-jpg', '-transparent');

%% Bar Graphs of Binned Layers (Binned Layers)
% Using bar graphs to represent monocular vs binocular response layer
% differences and fold change from one to the other

rcontrast = round(STIM.levels,2,'significant');

% subtraction calculations
STIM.calc.layers.subtractionDE.full = ((STIM.BIN.layers.full(:,:)-STIM.DE.layers.full(:,:)));
STIM.calc.layers.subtractionDE.transient = ((STIM.BIN.layers.transient(:,:)-STIM.DE.layers.transient(:,:)));
STIM.calc.layers.subtractionDE.sustained = ((STIM.BIN.layers.sustained(:,:)-STIM.DE.layers.sustained(:,:)));

% fold change calculations
STIM.calc.layers.fchange.full = ((STIM.BIN.layers.full(:,:)-STIM.DE.layers.full(:,:))./(STIM.DE.layers.full(:,:)));
STIM.calc.layers.fchange.transient = ((STIM.BIN.layers.transient(:,:)-STIM.DE.layers.transient(:,:))./(STIM.DE.layers.transient(:,:)));
STIM.calc.layers.fchange.sustained = ((STIM.BIN.layers.sustained(:,:)-STIM.DE.layers.sustained(:,:))./(STIM.DE.layers.sustained(:,:)));


figure('Position', [148,73,633,487]);
subplot(3,3,1)
bar(STIM.BIN.layers.full(:,1),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(STIM.DE.layers.full(:,1),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-5 45]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('percent change');
title('Supragranular');
hold off

subplot(3,3,2)
bar(STIM.BIN.layers.full(:,2),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(STIM.DE.layers.full(:,2),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-5 45]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('percent change');
title('Granular');
hold off

subplot(3,3,3)
bar(STIM.BIN.layers.full(:,3),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(STIM.DE.layers.full(:,3),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-5 45]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('percent change');
title('Infragranular');
hold off

subplot(3,3,4)
bar(STIM.calc.layers.subtractionDE.full(:,1),0.8,'FaceColor',[0.20, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-5 15]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('percent difference');
%title('Supragranular');
hold off

subplot(3,3,5)
bar(STIM.calc.layers.subtractionDE.full(:,2),0.8,'FaceColor',[0.20, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-5 15]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('percent difference');
%title('Granular');
hold off

subplot(3,3,6)
bar(STIM.calc.layers.subtractionDE.full(:,3),0.8,'FaceColor',[0.20, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-5 15]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('percent difference');
hold off

subplot(3,3,7)
bar(STIM.calc.layers.fchange.full(:,1),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-.2 1]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('fold change');
%title('Supragranular');
hold off

subplot(3,3,8)
bar(STIM.calc.layers.fchange.full(:,2),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-.2 1]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('fold change');
%title('Granular');
hold off

subplot(3,3,9)
bar(STIM.calc.layers.fchange.full(:,3),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-.2 1]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('fold change');
hold off

sgtitle({'Binned contacts by layer | aMUA responses',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_binned-layers-new',BRdatafile), '-jpg', '-transparent');

%% Transient vs sustained, Fold change Semilogx

figure('position',[360,450.3,560,167.6]);
subplot(1, 3, 1)
format bank;
semilogx(STIM.levels,STIM.calc.layers.fchange.transient(:,1),'-.k');
hold on
semilogx(STIM.levels,STIM.calc.layers.fchange.sustained(:,1),'k');
%semilogx(contrast,full(:,1),'-o');
xticklabels('');
ylim([-.2 1]);
xlabel('contrast level')
ylabel('Fold change');
title('Supragranular');
hold off

subplot(1, 3, 2)
format bank;
semilogx(STIM.levels,STIM.calc.layers.fchange.transient(:,2),'-.k');
hold on
semilogx(STIM.levels,STIM.calc.layers.fchange.sustained(:,2),'k');
xticklabels('');
xlabel('contrast level')
ylim([-.2 1]);
title('Granular');
hold off

subplot(1, 3, 3)
format bank;
semilogx(STIM.levels,STIM.calc.layers.fchange.transient(:,3), '-.k');
hold on
semilogx(STIM.levels,STIM.calc.layers.fchange.sustained(:,3),'k');
xticklabels('');
ylim([-.2 1]);
xlabel('contrast level')
title('Infragranular');
hold off
%legend('transient','sustained','full','Location','southoutside','orientation','vertical');

sgtitle({'Fold change from monocular to binocular: transient vs sustained',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_layers-transvsust',BRdatafile), '-jpg', '-transparent');

clear h1

%% Saving Workspace to D drive

Prompt = ('Would you like to save the workspace? (y/n)');
str = input(Prompt,'s');
if str == 'n' || str == 'N'
    error('Did not save workspace');
else 
    fprintf('\nSaving workspace...\n');
end

cd('D:\mcosinteroc\')
save(sprintf('%s',BRdatafile),'STIM','offset','BRdatafile','nct');
%savex(sprintf('%s',BRdatafile),'lfp','hpMUA','LFP','lpLFP','lpMUA','MUA');

fprintf('Workspace saved');