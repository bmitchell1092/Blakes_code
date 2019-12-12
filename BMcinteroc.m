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

BRdatafile = '160510_E_mcosinteroc001'; 
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

% sort the electrode contacts
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

STIM.LFP.raw  = trigData(LFP,STIM.onsetsdown,-pre,post); %pre variable is in absolute units 
STIM.CSD.raw  = trigData(CSD,STIM.onsetsdown,-pre,post); 
STIM.aMUA.raw = trigData(MUA,STIM.onsetsdown,-pre,post); 

%% Averaging across trials & baseline correct
clear avg
STIM.avg.LFP = mean(STIM.LFP.raw,3);
STIM.avg.aMUA = mean(STIM.aMUA.raw,3);
STIM.avg.CSD = mean(STIM.CSD.raw,3);

% These aren't currently used because I'm about to convert to percent
% change from baseline (next section). 

[STIM.LFP.bsl] = BMbasecorrect(STIM.avg.LFP); 
%[STIM.aMUA.bsl] = BMbasecorrect(STIM.avg.aMUA); 
[STIM.CSD.bsl] = BMbasecorrect(STIM.avg.CSD);

%% Data conversion 
% Converting the aMUA raw neural response to either percent change from
% baseline or z score (baseline as Xbar). 

% pre-allocate
STIM.aMUA.pc = nan(size(STIM.aMUA.raw,1),size(STIM.aMUA.raw,2),size(STIM.aMUA.raw,3));
STIM.aMUA.z = nan(size(STIM.aMUA.raw,1),size(STIM.aMUA.raw,2),size(STIM.aMUA.raw,3));
STIM.aMUA.bsl = nan(size(STIM.aMUA.raw,1),size(STIM.aMUA.raw,2),size(STIM.aMUA.raw,3));
% aMUA conversion to either z-score, percent change, or just baseline
% corrected raw.
clear t c
for t = 1:size(STIM.aMUA.raw,3)
    for c = 1:size(STIM.aMUA.raw,2)
        STIM.aMUA.z(:,c,t) = (STIM.aMUA.raw(:,c,t)-mean(STIM.aMUA.raw(25:75,c,t)))./(std(STIM.aMUA.raw(25:75,c,t))); %z score
        STIM.aMUA.pc(:,c,t) = (STIM.aMUA.raw(:,c,t)-mean(STIM.aMUA.raw(25:75,c,t)))./(mean(STIM.aMUA.raw(25:75,c,t)))*100; %percent change
        STIM.aMUA.bsl(:,c,t) = (STIM.aMUA.raw(:,c,t)-mean(STIM.aMUA.raw(25:75,c,t))); % baseline corrected raw
    end
end

%% Defining conditions of interest
% This is to create structures for the different conditions (including
% contrast levels and eye to which stimulus was shown. 

STIM.levels = unique(STIM.contrast);  % retrieves unique contrast levels
STIM.ori = unique(STIM.tilt); % retrieves unique orientations

clear i
for i = 1:length(STIM.levels)  % monocular (DE) contrast conditions
STIM.conditions.DE(i,:) = STIM.contrast == STIM.levels(i) & STIM.fixedc == 0; 
end

clear i 
for i = 1:length(STIM.levels)  % monocular (NDE) contrast conditions
STIM.conditions.NDE(i,:) = STIM.contrast == 0 & STIM.fixedc == STIM.levels(i); 
end

% dioptic contrast (same contrast in two eyes)
for i = 1:length(STIM.levels)  
STIM.conditions.BIN(i,:) = STIM.contrast == STIM.levels(i) & STIM.fixedc == STIM.levels(i); 
end

% dichoptic contrast (different contrast in two eyes)
STIM.conditions.DI(1,:) = STIM.contrast == STIM.levels(2) & STIM.fixedc == STIM.levels(3); % DE22NDE45
STIM.conditions.DI(2,:) = STIM.contrast == STIM.levels(2) & STIM.fixedc == STIM.levels(4); % DE22NDE90
STIM.conditions.DI(3,:) = STIM.contrast == STIM.levels(3) & STIM.fixedc == STIM.levels(2); % DE45NDE22
STIM.conditions.DI(4,:) = STIM.contrast == STIM.levels(3) & STIM.fixedc == STIM.levels(4); % DE45NDE90
STIM.conditions.DI(5,:) = STIM.contrast == STIM.levels(4) & STIM.fixedc == STIM.levels(2); % DE90NDE22
STIM.conditions.DI(6,:) = STIM.contrast == STIM.levels(4) & STIM.fixedc == STIM.levels(3); % DE90NDE45

%% Averaged trials by condition
% The above matrices were just logicals. Now time to take my full cMUA
% trials and use those logicals to trigger to and average by conditions of
% interest.

% aMUA by condition
clear m 
for m = 1:size(STIM.conditions.DE,1)
    STIM.DE.aMUA.pc(1).all(m,:,:) = mean(STIM.aMUA.pc(:,:,STIM.conditions.DE(m,:)),3); 
    STIM.NDE.aMUA.pc(1).all(m,:,:) = mean(STIM.aMUA.pc(:,:,STIM.conditions.NDE(m,:)),3); 
    STIM.BIN.aMUA.pc(1).all(m,:,:) = mean(STIM.aMUA.pc(:,:,STIM.conditions.BIN(m,:)),3);
end

clear m
for m = 1:size(STIM.conditions.DI,1)
    STIM.DI.aMUA.pc(1).all(m,:,:) = mean(STIM.aMUA.pc(:,:,STIM.conditions.DI(m,:)),3);
end

% CSD by condition
clear m 
for m = 1:size(STIM.conditions.DE,1)
    STIM.DE.CSD.raw(1).all(m,:,:) = mean(STIM.CSD.raw(:,:,STIM.conditions.DE(m,:)),3); 
    STIM.NDE.CSD.raw(1).all(m,:,:) = mean(STIM.CSD.raw(:,:,STIM.conditions.NDE(m,:)),3);
    STIM.BIN.CSD.raw(1).all(m,:,:) = mean(STIM.CSD.raw(:,:,STIM.conditions.BIN(m,:)),3);
end

for m = 1:size(STIM.conditions.DI,1)
    STIM.DI.CSD.raw(1).all(m,:,:) = mean(STIM.CSD.raw(:,:,STIM.conditions.DI(m,:)),3);
end

%% Binning into layers

laminae_fields = fieldnames(STIM.laminae);
clear l
for L = 1:length(laminae_fields)
STIM.DE.aMUA.pc(1).layers(:,:,L) = mean(STIM.DE.aMUA.pc(1).all(:,:,STIM.laminae.(laminae_fields{L})),3);
STIM.NDE.aMUA.pc(1).layers(:,:,L) = mean(STIM.NDE.aMUA.pc(1).all(:,:,STIM.laminae.(laminae_fields{L})),3);
STIM.BIN.aMUA.pc(1).layers(:,:,L) = mean(STIM.BIN.aMUA.pc(1).all(:,:,STIM.laminae.(laminae_fields{L})),3);
STIM.DI.aMUA.pc(1).layers(:,:,L) = mean(STIM.DI.aMUA.pc(1).all(:,:,STIM.laminae.(laminae_fields{L})),3);
end

%% collapsing across time for each condition
% In order to get a single number to represent the contrast response for each contact, I can collapse
% across time. The following code collapses across different lengths of time:
% full stim duration (80ms to offset), transient (80 to 150ms), and
% sustained (151ms to offset). 
full = 80:offset;
transient = 80:150;
sustained = 151:offset;

% Collapsed aMUA percent change
clear i 
for i=1:size(STIM.DE.aMUA.pc(1).all,1)
    STIM.DE.aMUA.pc(1).coll(i,1,:)  = mean(STIM.DE.aMUA.pc(1).all(i,full,:),2);
    STIM.NDE.aMUA.pc(1).coll(i,1,:)  = mean(STIM.NDE.aMUA.pc(1).all(i,full,:),2);
    STIM.BIN.aMUA.pc(1).coll(i,1,:)  = mean(STIM.BIN.aMUA.pc(1).all(i,full,:),2);
    STIM.DE.aMUA.pc(1).coll(i,2,:)  = mean(STIM.DE.aMUA.pc(1).all(i,transient,:),2);
    STIM.NDE.aMUA.pc(1).coll(i,2,:)  = mean(STIM.NDE.aMUA.pc(1).all(i,transient,:),2);
    STIM.BIN.aMUA.pc(1).coll(i,2,:)  = mean(STIM.BIN.aMUA.pc(1).all(i,transient,:),2);
    STIM.DE.aMUA.pc(1).coll(i,3,:)  = mean(STIM.DE.aMUA.pc(1).all(i,sustained,:),2);
    STIM.NDE.aMUA.pc(1).coll(i,3,:)  = mean(STIM.NDE.aMUA.pc(1).all(i,sustained,:),2);
    STIM.BIN.aMUA.pc(1).coll(i,3,:)  = mean(STIM.BIN.aMUA.pc(1).all(i,sustained,:),2);
end

clear i
for i=1:size(STIM.DI.aMUA.pc(1).all,1)
    STIM.DI.aMUA.pc(1).coll(i,1,:)  = mean(STIM.DI.aMUA.pc(1).all(i,full,:),2);
    STIM.DI.aMUA.pc(1).coll(i,2,:)  = mean(STIM.DI.aMUA.pc(1).all(i,transient,:),2);
    STIM.DI.aMUA.pc(1).coll(i,3,:)  = mean(STIM.DI.aMUA.pc(1).all(i,sustained,:),2);
end

%% Binning contacts into V1 layers (informed by datalogs)
% Seperated into full stim duration, transient response, and sustained
% response

laminae_fields = fieldnames(STIM.laminae);
clear l
for L = 1:length(laminae_fields)
STIM.DE.aMUA.pc(1).coll_layers(:,:,L) = mean(STIM.DE.aMUA.pc(1).coll(:,:,STIM.laminae.(laminae_fields{L})),3);
STIM.NDE.aMUA.pc(1).coll_layers(:,:,L) = mean(STIM.NDE.aMUA.pc(1).coll(:,:,STIM.laminae.(laminae_fields{L})),3);
STIM.BIN.aMUA.pc(1).coll_layers(:,:,L) = mean(STIM.BIN.aMUA.pc(1).coll(:,:,STIM.laminae.(laminae_fields{L})),3);
STIM.DI.aMUA.pc(1).coll_layers(:,:,L) = mean(STIM.DI.aMUA.pc(1).coll(:,:,STIM.laminae.(laminae_fields{L})),3);
end

%% Calculations (STIM.CALC)

% Subtractions: coll
STIM.calc.aMUA.pc.subtract.BIN_DE.coll = STIM.BIN.aMUA.pc.coll(:,:,:)-STIM.DE.aMUA.pc.coll(:,:,:);
STIM.calc.aMUA.pc.subtract.BIN_NDE.coll = STIM.BIN.aMUA.pc.coll(:,:,:)-STIM.NDE.aMUA.pc.coll(:,:,:);
STIM.calc.aMUA.pc.subtract.DE_NDE.coll = STIM.DE.aMUA.pc.coll(:,:,:)-STIM.NDE.aMUA.pc.coll(:,:,:);

% Subtractions: coll_layer
STIM.calc.aMUA.pc.subtract.BIN_DE.coll_layers = (STIM.BIN.aMUA.pc.coll_layers(:,:,:)-STIM.DE.aMUA.pc.coll_layers(:,:,:));
STIM.calc.aMUA.pc.subtract.BIN_NDE.coll_layers = (STIM.BIN.aMUA.pc.coll_layers(:,:,:)-STIM.NDE.aMUA.pc.coll_layers(:,:,:));
STIM.calc.aMUA.pc.subtract.DE_NDE.coll_layers = (STIM.DE.aMUA.pc.coll_layers(:,:,:)-STIM.NDE.aMUA.pc.coll_layers(:,:,:));

% Fold change: coll_layer
STIM.calc.aMUA.pc.fchange.BIN_DE = ((STIM.BIN.aMUA.pc.coll_layers(:,:,:)-STIM.DE.aMUA.pc.coll_layers(:,:,:))./(STIM.DE.aMUA.pc.coll_layers(:,:,:)));
STIM.calc.aMUA.pc.fchange.BIN_NDE = ((STIM.BIN.aMUA.pc.coll_layers(:,:,:)-STIM.DE.aMUA.pc.coll_layers(:,:,:))./(STIM.NDE.aMUA.pc.coll_layers(:,:,:)));
STIM.calc.aMUA.pc.fchange.DE_NDE = ((STIM.DE.aMUA.pc.coll_layers(:,:,:)-STIM.NDE.aMUA.pc.coll_layers(:,:,:))./(STIM.NDE.aMUA.pc.coll_layers(:,:,:)));

%% Models (LSM | QSM)
% Model Prediction: Linear Summation (LSM)

% LSMvsBIN
STIM.BIN.aMUA.pc_LSM.all = STIM.DE.aMUA.pc.all(:,:,:)+STIM.NDE.aMUA.pc.all(:,:,:);
STIM.BIN.aMUA.pc_LSM.layers = STIM.DE.aMUA.pc.layers(:,:,:)+STIM.NDE.aMUA.pc.layers(:,:,:);
STIM.BIN.aMUA.pc_LSM.coll = STIM.DE.aMUA.pc.coll(:,:,:)+STIM.NDE.aMUA.pc.coll(:,:,:);
STIM.BIN.aMUA.pc_LSM.coll_layers = STIM.DE.aMUA.pc.coll_layers(:,:,:)+STIM.NDE.aMUA.pc.coll_layers(:,:,:);

% LSMvsDI

sDE = [2 2 3 3 4 4]; % code of numbers to loop thru for dichoptic conditions
sNDE = [3 4 2 4 2 3];

clear i
for i = 1:length(sDE)
STIM.DI.aMUA.pc_LSM.all(i,:,:) = STIM.DE.aMUA.pc.all(sDE(i),:,:) + STIM.NDE.aMUA.pc.all(sNDE(i),:,:);
STIM.DI.aMUA.pc_LSM.layers(i,:,:) = STIM.DE.aMUA.pc.layers(sDE(i),:,:) + STIM.NDE.aMUA.pc.layers(sNDE(i),:,:);
STIM.DI.aMUA.pc_LSM.coll(i,:,:) = STIM.DE.aMUA.pc.coll(sDE(i),:,:) + STIM.NDE.aMUA.pc.coll(sNDE(i),:,:);
STIM.DI.aMUA.pc_LSM.coll_layers(i,:,:) = STIM.DE.aMUA.pc.coll_layers(sDE(i),:,:) + STIM.NDE.aMUA.pc.coll_layers(sNDE(i),:,:);
end

% Model Prediction: QSM (Quadratic summation)
% QSMvsBIN
STIM.BIN.aMUA.pc_QSM.all = sqrt((STIM.DE.aMUA.pc.all(:,:,:).^2) + (STIM.NDE.aMUA.pc.all(:,:,:).^2));
STIM.BIN.aMUA.pc_QSM.layers = sqrt((STIM.DE.aMUA.pc.layers(:,:,:).^2) + (STIM.NDE.aMUA.pc.layers(:,:,:).^2));
STIM.BIN.aMUA.pc_QSM.coll = sqrt((STIM.DE.aMUA.pc.coll(:,:,:).^2) + (STIM.NDE.aMUA.pc.coll(:,:,:).^2));
STIM.BIN.aMUA.pc_QSM.coll_layers = sqrt((STIM.DE.aMUA.pc.coll_layers(:,:,:).^2) + (STIM.NDE.aMUA.pc.coll_layers(:,:,:).^2));

% QSMvsDI
for i = 1:length(sDE)
STIM.DI.aMUA.pc_QSM.all(i,:,:) = sqrt((STIM.DE.aMUA.pc.all(sDE(i),:,:).^2) + (STIM.NDE.aMUA.pc.all(sNDE(i),:,:).^2));
STIM.DI.aMUA.pc_QSM.layers(i,:,:) = sqrt((STIM.DE.aMUA.pc.layers(sDE(i),:,:).^2) + (STIM.NDE.aMUA.pc.layers(sNDE(i),:,:).^2));
STIM.DI.aMUA.pc_QSM.coll(i,:,:) = sqrt((STIM.DE.aMUA.pc.coll(sDE(i),:,:).^2) + (STIM.NDE.aMUA.pc.coll(sNDE(i),:,:).^2));
STIM.DI.aMUA.pc_QSM.coll_layers(i,:,:) = sqrt((STIM.DE.aMUA.pc.coll_layers(sDE(i),:,:).^2) + (STIM.NDE.aMUA.pc.coll_layers(sNDE(i),:,:).^2));
end
%%

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

bAVG_iCSD = filterCSD(STIM.CSD.bsl')';

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
f_ShadedLinePlotbyDepth(mean(STIM.aMUA.raw(:,:,STIM.conditions.DE(numel(STIM.levels),:)),3),0:(1/(numel(STIM.channels))):1,STIM.refwin,STIM.channels,1,1);
title('.90 contrast in DE');
xlabel('time (ms)');
hold off

clear i 
for i = 1:length(STIM.levels)-1
    subplot(1,length(STIM.levels),i);
    f_ShadedLinePlotbyDepth_BAM(mean(STIM.aMUA.raw(:,:,STIM.conditions.DE(i,:)),3),0:(1/(numel(STIM.channels))):1,STIM.refwin,STIM.channels,1,1,false,scalingfactor);
   
plot([0 0], ylim,'k')
plot([offset offset], ylim,'k')
if i == 1
    title({STIM.levels(i),' contrast in both eyes'});
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
    bar(STIM.BIN.aMUA.pc.coll(:,2,selectchannels(c)),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
    hold on
    bar(STIM.DE.aMUA.pc.coll(:,2,selectchannels(c)),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
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
    plot(flipud(squeeze(STIM.DE.aMUA.pc.coll(c+1,2,:))),STIM.channels,'b','linewidth',.5);
    hold on
    plot(flipud(squeeze(STIM.BIN.aMUA.pc.coll(c+1,2,:))),STIM.channels,'.-r','linewidth',0.5);
    hline(STIM.channels(end)-STIM.laminae.gran(end),'-.')
    xlim([-10 100])
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
    plot(flipud(squeeze(STIM.DE.aMUA.pc.coll(c-2,3,:))),STIM.channels,'b','linewidth',.5);
    hold on
    plot(flipud(squeeze(STIM.BIN.aMUA.pc.coll(c-2,3,:))),STIM.channels,'.-r','linewidth',0.5);
    hline(STIM.channels(end)-STIM.laminae.gran(end),'-.')
    xlim([-10 100])
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
%cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
%export_fig(sprintf('%s_tightplot',BRdatafile), '-jpg', '-transparent');

%% Collapsed contrast responses for monocular and binocular stim (contrast lines)
% Designed to show Monocular vs binocular contrast responses alongside CSD
% for layer reference. 
figure('position',[185 150 887 450]);

h = subplot(1,4,1);
bAVG_iCSD = filterCSD(STIM.CSD.bsl')';
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
plot(fliplr(squeeze(STIM.BIN.aMUA.pc.coll(:,2,:))),STIM.channels);
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
plot(fliplr(squeeze(STIM.DE.aMUA.pc.coll(:,2,:))),STIM.channels);
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

subplot(1,4,4)
plot(fliplr(squeeze(STIM.calc.aMUA.pc.subtract.BIN_DE.coll(:,2,:))),STIM.channels);
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

clear c L
for L = 1:3
    subplot(2,3,L)
    plot(STIM.refwin,STIM.DE.aMUA.pc.layers(:,:,L),'color','b')
    hold on
    plot(STIM.refwin,STIM.BIN.aMUA.pc.layers(:,:,L),'color','r')
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
    plot(STIM.refwin,((smooth(STIM.BIN.aMUA.pc.layers(c,:,L)))-(smooth(STIM.DE.aMUA.pc.layers(c,:,L)))),'linewidth',.1);
    hold on
    %ylimit = max(abs(get(gcf,'ylim')));
    ylimit = 50;
    set(gca,'ylim',[-10 ylimit],'Box','off','TickDir','out')
    ylabel({'Percent change'...
    'from baseline'});
    xlabel('time (ms)');
    end
end

sgtitle({'V1 laminae contrast responses: monoptic (blue) vs dioptic (red) stimulation'...
   ,BRdatafile},'Interpreter','none');

% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_timeplots',BRdatafile), '-jpg', '-transparent');

%% Bar Graph: Coll_layer; DE vs NDE vs BIN
% Using bar graphs to represent monocular vs binocular response layer
% differences and fold change from one to the other

rcontrast = round(STIM.levels,2,'significant');

figure('Position', [155,98,965,487]);

labels = {0, .22, .45, .9,[],0,.22,.45,.9}; format bank;

clear L c
for L = 1:3
subplot(2,3,L)
bar([STIM.BIN.aMUA.pc.coll_layers(:,2,L);NaN;STIM.BIN.aMUA.pc.coll_layers(:,3,L)],0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar([STIM.DE.aMUA.pc.coll_layers(:,2,L);NaN;STIM.DE.aMUA.pc.coll_layers(:,3,L)],0.6,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
bar([STIM.NDE.aMUA.pc.coll_layers(:,2,L);NaN;STIM.NDE.aMUA.pc.coll_layers(:,3,L)],0.4,'FaceColor',[.2, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-5 55]);
xticklabels(labels)
xlabel('contrast')
ylabel('Relative aMUA response');
    if L == 1
        title('Supragranular');
        legend('BIN','DE','NDE','location','northeast');
    elseif L == 2
            title('Granular');
        else 
            title('Infragranular');
    end
hold off
end

clear L c
for L = 1:3
subplot(2,3,L+3)
bar([STIM.calc.aMUA.pc.subtract.BIN_NDE.coll_layers(:,2,L);NaN;STIM.calc.aMUA.pc.subtract.BIN_NDE.coll_layers(:,3,L)],0.4,'FaceColor',[0.2500, 0.250, 0.2980],'EdgeColor','k','LineWidth',0.8);
hold on
bar([STIM.calc.aMUA.pc.subtract.BIN_DE.coll_layers(:,2,L);NaN;STIM.calc.aMUA.pc.subtract.BIN_DE.coll_layers(:,3,L)],0.8,'FaceColor',[0.8500, 0.5250, 0.7980],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-5 30]);
xticklabels(labels)
xlabel('contrast')
ylabel('difference');
hold off
end

sgtitle({'Binned contacts by layer | aMUA responses',BRdatafile},'Interpreter','none');

% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_binned-layers-new',BRdatafile), '-jpg', '-transparent');

%% Bar Graphs: Model Predictions vs BINOCULAR response

figure('Position', [155,98,965,487]);

clear L c
for L = 1:3
subplot(3,3,L)
bar([STIM.DE.aMUA.pc.coll_layers(:,2,L);NaN;STIM.DE.aMUA.pc.coll_layers(:,3,L)],0.8,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
hold on
bar([STIM.NDE.aMUA.pc.coll_layers(:,2,L);NaN;STIM.NDE.aMUA.pc.coll_layers(:,3,L)],0.6,'FaceColor',[.2, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-5 80]);
xticklabels(labels)
xlabel('contrast')
ylabel('aMUA response');
    if L == 1
        title('Supragranular');
        lgd = legend('DE','NDE','location','northwest');
        lgd.FontSize = 4;
    elseif L == 2
            title('Granular');
        else 
            title('Infragranular');
    end
hold off
end

clear L c
for L = 1:3
subplot(3,3,L+3)
bar([STIM.BIN.aMUA.pc_LSM.coll_layers(:,2,L);NaN;STIM.BIN.aMUA.pc_LSM.coll_layers(:,3,L)],0.8,'FaceColor',[1, 1, 1],'linestyle','--','EdgeColor','k','LineWidth',0.8);
hold on
bar([STIM.BIN.aMUA.pc.coll_layers(:,2,L);NaN;STIM.BIN.aMUA.pc.coll_layers(:,3,L)],0.6,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
bar([STIM.BIN.aMUA.pc_QSM.coll_layers(:,2,L);NaN;STIM.BIN.aMUA.pc_QSM.coll_layers(:,3,L)],0.4,'FaceColor',[.85, .325, .098],'linestyle',':','EdgeColor','k','LineWidth',1);
set(gca,'box','off');
ylim([-5 80]);
xticklabels(labels)
xlabel('contrast')
ylabel('aMUA response');
hold off
    if L == 1
        lgd = legend('LSM','BIN','QSM','location','northwest');
        lgd.FontSize = 4;
    end
end

clear i
for i = 1:3
subplot(3,3,i+6)
bar([STIM.BIN.aMUA.pc_LSM.coll_layers(:,2,L)-STIM.BIN.aMUA.pc.coll_layers(:,2,L);NaN;STIM.BIN.aMUA.pc_LSM.coll_layers(:,3,L)-STIM.BIN.aMUA.pc.coll_layers(:,3,L)],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',0.8,'linestyle','--');
hold on
bar([STIM.BIN.aMUA.pc_QSM.coll_layers(:,2,L)-STIM.BIN.aMUA.pc.coll_layers(:,2,L);NaN;STIM.BIN.aMUA.pc_QSM.coll_layers(:,3,L)-STIM.BIN.aMUA.pc.coll_layers(:,3,L)],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',1,'linestyle',':');
set(gca,'box','off');
ylim([-15 30]);
xticklabels(labels);
xlabel('contrast')
ylabel('difference from model');
hold off
end

sgtitle({'Binned contacts by layer | aMUA responses',BRdatafile},'Interpreter','none');

%%

figure('position',[145,88.33333333333333,786.6666666666666,541.3333333333333]);
clear i L
for L = 1:3
subplot(3,3,L)
plot(smooth(STIM.refwin,STIM.BIN.aMUA.pc_QSM.layers(2,:,L)-STIM.BIN.aMUA.pc.layers(2,:,L),0.1,'rloess'),'-b','linewidth',2);
hold on
plot(smooth(STIM.refwin,STIM.BIN.aMUA.pc_LSM.layers(2,:,L)-STIM.BIN.aMUA.pc.layers(2,:,L),0.1,'rloess'),'-r','linewidth',2);
ylim([-20 40])
xlim([0 600])
xlabel('time (ms)');
ylabel('difference');
hline(0,'-.k')
vline(offset,'-.k')
    if L == 1
        title('Supragranular');
        lgd = legend('QSM','LSM','location','northwest');
        lgd.FontSize = 4;
    elseif L == 2
            title('Granular');
        else 
            title('Infragranular');
    end
end

clear i L
for L = 1:3
subplot(3,3,L+3)
plot(smooth(STIM.refwin,STIM.BIN.aMUA.pc_QSM.layers(3,:,L)-STIM.BIN.aMUA.pc.layers(3,:,L),0.1,'rloess'),'-b','linewidth',2);
hold on
plot(smooth(STIM.refwin,STIM.BIN.aMUA.pc_LSM.layers(3,:,L)-STIM.BIN.aMUA.pc.layers(3,:,L),0.1,'rloess'),'-r','linewidth',2);
ylim([-20 40])
xlim([0 600])
xlabel('time (ms)');
ylabel('difference');
hline(0,'-.k')
vline(offset,'k')
end

clear i L
for L = 1:3
subplot(3,3,L+6)
plot(smooth(STIM.refwin,STIM.BIN.aMUA.pc_QSM.layers(4,:,L)-STIM.BIN.aMUA.pc.layers(4,:,L),0.1,'rloess'),'-b','linewidth',2);
hold on
plot(smooth(STIM.refwin,STIM.BIN.aMUA.pc_LSM.layers(4,:,L)-STIM.BIN.aMUA.pc.layers(4,:,L),0.1,'rloess'),'-r','linewidth',2);
ylim([-20 40])
xlim([0 600])
xlabel('time (ms)');
ylabel('difference');
hline(0,'-.k')
vline(offset,'k')
hold off
end

sgtitle("idk");

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