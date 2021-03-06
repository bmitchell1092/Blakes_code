%% BMcinteroc
% script to: select an input file, load in .nev,
% pull out stim onset, load in .ns2 and .ns6 and sync event codes to the
% raw neural data. Each stim onset with animal fixation for the duration of presentation is a trial. 
% Generate LFP, aMUA, and CSD --triggered to a reference
% window. Average across trials and baseline corrected when necessary. 
clear

loop = 0; % 1 = loop thru all files; 0 requires file to be specified
%% Start

if loop == true
    myFolder = 'D:\DATA (unused)\';  % Specify the folder where the files live.

    % Check to make sure that folder actually exists.  Warn user if it doesn't.
    if ~isfolder(myFolder)
      errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
      uiwait(warndlg(errorMessage));
      return;
    end

    % Get a list of all files in the folder with the desired file name pattern.
    filePattern = fullfile(myFolder, '*.ns6'); 
    matFiles = dir(filePattern);
    for k = 1:length(matFiles)
    baseFileName{k} = matFiles(k).name;
    fullFileName{k} = fullfile(myFolder, baseFileName{k});
    [~, name{k}, ~] = fileparts(baseFileName{k});
    end
end

%% Optional (Loop thru all sessions on D: Drive)
clear z
% for z = 2:6



%% Establish directories and set path

cd('C:\Users\bmitc\')

if strcmp(getenv('USER'),'maierav')                                      %retrieves environment variable 'USER' 
    npmkdir  = '/Users/alex 1/Desktop/LAB/Brock/OLD/NPMK-4.5.3.0/NPMK/'; %directory for Alex's machine
    nbanalysisdir   = '/Users/alex 1/Desktop/LAB/bootcamp/nbanalysis/';  %directory for Alex's machine
    datadir  = '/Users/alex 1/Desktop/LAB/';                             %directory for the stored data
else
    npmkdir  = '/users/bmitc/Documents/MATLAB/NPMK/';                    %neural processing matlab kit (NPMK)
    nbanalysisdir   = '/users/bmitc/Documents/MATLAB/nbanalysis/';       %directory with various tools for opening, loading, and processing 
    %datadir  = '/users/bmitc/Box Sync/DATA/';                            %this is my Vanderbilt Box sync 
    %datadir = 'users/bmitc/Documents/MATLAB/data/';
    datadir = 'D:\DATA (unused)\';
end

addpath(genpath(npmkdir))
addpath(genpath(nbanalysisdir))
addpath(genpath(datadir))

if loop == false
    BRdatafile = '160427_E_mcosinteroc001'; 
else 

    BRdatafile = char(name(z));
end

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
        MUA = MUA(:,flipud(idx)); % kacie used fliplr here but it wasn't flipping properly for me 
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
STIM.off = round(((STIM.offsets(1)-STIM.onsets(1))/(30)),0); % stimulus offset as calculated from grating text file.

pre   = -50; % 50ms before stim onset 
post = 600; % not all stim offsets are the same, but this is a consistent window

STIM.LFP.raw = trigData(LFP,STIM.onsetsdown,-pre,post); %pre variable is in absolute units 
STIM.CSD.raw  = trigData(CSD,STIM.onsetsdown,-pre,post); 
STIM.aMUA.raw = trigData(MUA,STIM.onsetsdown,-pre,post); 

%% Reference window (for plotting)

STIM.refwin = pre:post;

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
% Converting the aMUA raw neural response to either percent change or z score (baseline as Xbar). 

% pre-allocate
STIM.aMUA.pc = nan(size(STIM.aMUA.raw,1),size(STIM.aMUA.raw,2),size(STIM.aMUA.raw,3));
STIM.aMUA.z = nan(size(STIM.aMUA.raw,1),size(STIM.aMUA.raw,2),size(STIM.aMUA.raw,3));
STIM.aMUA.bsl = nan(size(STIM.aMUA.raw,1),size(STIM.aMUA.raw,2),size(STIM.aMUA.raw,3));

% aMUA conversion
clear t c
for t = 1:size(STIM.aMUA.raw,3)
    for c = 1:size(STIM.aMUA.raw,2)
        STIM.aMUA.z(:,c,t) = (STIM.aMUA.raw(:,c,t)-mean(STIM.aMUA.raw(25:75,c,t)))./(std(STIM.aMUA.raw(25:75,c,t))); %z score
        STIM.aMUA.pc(:,c,t) = (STIM.aMUA.raw(:,c,t)-mean(STIM.aMUA.raw(25:75,c,t)))./(mean(STIM.aMUA.raw(25:75,c,t)))*100; %percent change
        STIM.aMUA.bsl(:,c,t) = (STIM.aMUA.raw(:,c,t)-mean(STIM.aMUA.raw(25:75,c,t))); % baseline corrected raw
    end
end

%% Defining conditions of interest
% This is to create logical matrices for the different conditions
% Independent variables of interest: contrast, ocularity, and orientation 
% michele : here's how i'm currently defining stimulus conditions: 2
% monocular conditions (DE and NDE), two binocular conditions (BIN and DI).
% Currently, orientation and phase are being disregarded  

STIM.levels = unique(STIM.contrast);  % retrieves unique contrast levels
STIM.ori = unique(STIM.tilt); % retrieves unique orientations

clear i
for i = 1:length(STIM.levels)  % monocular (DE) contrast conditions
STIM.conditions.DE(i,:) = STIM.contrast == STIM.levels(i) ...
    & STIM.fixedc == 0; % STIM.contrast is dominant eye; STIM.fixedc is non-dominant eye; 
    %& STIM.tilt == STIM.ori(2);
end

clear i 
for i = 1:length(STIM.levels)  % monocular (NDE) contrast conditions
STIM.conditions.NDE(i,:) = STIM.contrast == 0 ...
    & STIM.fixedc == STIM.levels(i);
    %& STIM.tilt == STIM.ori(2);
 
end

% dioptic contrast (same contrast in two eyes)
for i = 1:length(STIM.levels)  
STIM.conditions.BIN(i,:) = STIM.contrast == STIM.levels(i) ...
    & STIM.fixedc == STIM.levels(i);
    %& STIM.tilt == STIM.ori(2);

end

% dichoptic contrast (different contrast in two eyes); in its current state, does not control for
% different orientations. 
STIM.conditions.DI(1,:) = STIM.contrast == STIM.levels(2) & STIM.fixedc == STIM.levels(3); % DE22 NDE45
STIM.conditions.DI(2,:) = STIM.contrast == STIM.levels(3) & STIM.fixedc == STIM.levels(2); % DE45 NDE22
STIM.conditions.DI(3,:) = STIM.contrast == STIM.levels(2) & STIM.fixedc == STIM.levels(4); % DE22 NDE90
STIM.conditions.DI(4,:) = STIM.contrast == STIM.levels(4) & STIM.fixedc == STIM.levels(2); % DE90 NDE22
STIM.conditions.DI(5,:) = STIM.contrast == STIM.levels(3) & STIM.fixedc == STIM.levels(4); % DE45 NDE90
STIM.conditions.DI(6,:) = STIM.contrast == STIM.levels(4) & STIM.fixedc == STIM.levels(3); % DE90 NDE45

%% Averaged trials by condition
% The above code defined logicals for each condition. The following will 
% use those logicals to create corresponding data matrices
% michele2 : Below, I create a new field in the STIM structure to average
% across all the trials with the stimulus features for which I made logical vectors in the
% previous section. 

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
    STIM.DE.CSD.raw(m,:,:) = mean(STIM.CSD.raw(:,:,STIM.conditions.DE(m,:)),3); 
    STIM.NDE.CSD.raw(m,:,:) = mean(STIM.CSD.raw(:,:,STIM.conditions.NDE(m,:)),3);
    STIM.BIN.CSD.raw(m,:,:) = mean(STIM.CSD.raw(:,:,STIM.conditions.BIN(m,:)),3);
end

for m = 1:size(STIM.conditions.DI,1)
    STIM.DI.CSD.raw(m,:,:) = mean(STIM.CSD.raw(:,:,STIM.conditions.DI(m,:)),3);
end

% BSL CSD
        STIM.DE.CSD.bsl(:,:,:) = (STIM.DE.CSD.raw(:,:,:)-mean(STIM.DE.CSD.raw(:,25:75,:),2)); % baseline corrected raw
        STIM.NDE.CSD.bsl(:,:,:) = (STIM.NDE.CSD.raw(:,:,:)-mean(STIM.NDE.CSD.raw(:,25:75,:),2)); % baseline corrected raw
        STIM.BIN.CSD.bsl(:,:,:) = (STIM.BIN.CSD.raw(:,:,:)-mean(STIM.BIN.CSD.raw(:,25:75,:),2)); % baseline corrected raw
        STIM.DI.CSD.bsl(:,:,:) = (STIM.DI.CSD.raw(:,:,:)-mean(STIM.DI.CSD.raw(:,25:75,:),2)); % baseline corrected raw

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
full = 80:STIM.off;
transient = 80:150;
sustained = 151:STIM.off;

% Collapsed aMUA percent change
clear i 
for i=1:size(STIM.DE.aMUA.pc(1).all,1) % probably don't need this for loop
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

laminae_fields = fieldnames(STIM.laminae);
clear l
for L = 1:length(laminae_fields)
STIM.DE.aMUA.pc(1).coll_layers(:,:,L) = mean(STIM.DE.aMUA.pc(1).coll(:,:,STIM.laminae.(laminae_fields{L})),3);
STIM.NDE.aMUA.pc(1).coll_layers(:,:,L) = mean(STIM.NDE.aMUA.pc(1).coll(:,:,STIM.laminae.(laminae_fields{L})),3);
STIM.BIN.aMUA.pc(1).coll_layers(:,:,L) = mean(STIM.BIN.aMUA.pc(1).coll(:,:,STIM.laminae.(laminae_fields{L})),3);
STIM.DI.aMUA.pc(1).coll_layers(:,:,L) = mean(STIM.DI.aMUA.pc(1).coll(:,:,STIM.laminae.(laminae_fields{L})),3);
end

%% Calculations (STIM.CALC)
% These calculations are less important, considering I don't necessarily
% need to store this information (I can just do it for plots).

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

%% Models (LSM | QSM | NRM)
% Model Prediction: Linear Summation (LSM)

% LSMvsBIN
STIM.BIN.aMUA.pc_LSM.all = STIM.DE.aMUA.pc.all(:,:,:)+STIM.NDE.aMUA.pc.all(:,:,:);
STIM.BIN.aMUA.pc_LSM.layers = STIM.DE.aMUA.pc.layers(:,:,:)+STIM.NDE.aMUA.pc.layers(:,:,:);
STIM.BIN.aMUA.pc_LSM.coll = STIM.DE.aMUA.pc.coll(:,:,:)+STIM.NDE.aMUA.pc.coll(:,:,:);
STIM.BIN.aMUA.pc_LSM.coll_layers = STIM.DE.aMUA.pc.coll_layers(:,:,:) + STIM.NDE.aMUA.pc.coll_layers(:,:,:);
STIM.BIN.CSD.LSM = STIM.DE.CSD.bsl(:,:,:) + STIM.NDE.CSD.bsl(:,:,:);

% LSMvsDI

sDE = [2 3 2 4 3 4]; % code of numbers to loop thru for dichoptic conditions
sNDE = [3 2 4 2 4 3];

clear i
for i = 1:length(sDE)
STIM.DI.aMUA.pc_LSM.all(i,:,:) = STIM.DE.aMUA.pc.all(sDE(i),:,:) + STIM.NDE.aMUA.pc.all(sNDE(i),:,:);
STIM.DI.aMUA.pc_LSM.layers(i,:,:) = STIM.DE.aMUA.pc.layers(sDE(i),:,:) + STIM.NDE.aMUA.pc.layers(sNDE(i),:,:);
STIM.DI.aMUA.pc_LSM.coll(i,:,:) = STIM.DE.aMUA.pc.coll(sDE(i),:,:) + STIM.NDE.aMUA.pc.coll(sNDE(i),:,:);
STIM.DI.aMUA.pc_LSM.coll_layers(i,:,:) = STIM.DE.aMUA.pc.coll_layers(sDE(i),:,:) + STIM.NDE.aMUA.pc.coll_layers(sNDE(i),:,:);
STIM.DI.CSD.LSM(i,:,:) = (STIM.DE.CSD.bsl(sDE(i),:,:) + STIM.NDE.CSD.bsl(sNDE(i),:,:));
end

% Model Prediction: QSM (Quadratic summation)
% QSMvsBIN
STIM.BIN.aMUA.pc_QSM.all = sqrt((STIM.DE.aMUA.pc.all(:,:,:).^2) + (STIM.NDE.aMUA.pc.all(:,:,:).^2));
STIM.BIN.aMUA.pc_QSM.layers = sqrt((STIM.DE.aMUA.pc.layers(:,:,:).^2) + (STIM.NDE.aMUA.pc.layers(:,:,:).^2));
STIM.BIN.aMUA.pc_QSM.coll = sqrt((STIM.DE.aMUA.pc.coll(:,:,:).^2) + (STIM.NDE.aMUA.pc.coll(:,:,:).^2));
STIM.BIN.aMUA.pc_QSM.coll_layers = sqrt((STIM.DE.aMUA.pc.coll_layers(:,:,:).^2) + (STIM.NDE.aMUA.pc.coll_layers(:,:,:).^2));
STIM.BIN.CSD.QSM = -sqrt((STIM.DE.CSD.bsl(:,:,:).^2) + (STIM.NDE.CSD.bsl(:,:,:).^2));

% QSMvsDI
clear i
for i = 1:length(sDE)
STIM.DI.aMUA.pc_QSM.all(i,:,:) = sqrt((STIM.DE.aMUA.pc.all(sDE(i),:,:).^2) + (STIM.NDE.aMUA.pc.all(sNDE(i),:,:).^2));
STIM.DI.aMUA.pc_QSM.layers(i,:,:) = sqrt((STIM.DE.aMUA.pc.layers(sDE(i),:,:).^2) + (STIM.NDE.aMUA.pc.layers(sNDE(i),:,:).^2));
STIM.DI.aMUA.pc_QSM.coll(i,:,:) = sqrt((STIM.DE.aMUA.pc.coll(sDE(i),:,:).^2) + (STIM.NDE.aMUA.pc.coll(sNDE(i),:,:).^2));
STIM.DI.aMUA.pc_QSM.coll_layers(i,:,:) = sqrt((STIM.DE.aMUA.pc.coll_layers(sDE(i),:,:).^2) + (STIM.NDE.aMUA.pc.coll_layers(sNDE(i),:,:).^2));
STIM.DI.CSD.QSM(i,:,:) = sqrt((STIM.DE.CSD.bsl(sDE(i),:,:).^2) + (STIM.NDE.CSD.bsl(sNDE(i),:,:).^2));
end

% NRMvsBIN

STIM.BIN.aMUA.pc_NRM.all = (STIM.DE.aMUA.pc.all(:,:,:)+STIM.NDE.aMUA.pc.all(:,:,:))./2;
STIM.BIN.aMUA.pc_NRM.layers = (STIM.DE.aMUA.pc.layers(:,:,:)+STIM.NDE.aMUA.pc.layers(:,:,:))./2;
STIM.BIN.aMUA.pc_NRM.coll = (STIM.DE.aMUA.pc.coll(:,:,:)+STIM.NDE.aMUA.pc.coll(:,:,:))./2;
STIM.BIN.aMUA.pc_NRM.coll_layers = (STIM.DE.aMUA.pc.coll_layers(:,:,:)+STIM.NDE.aMUA.pc.coll_layers(:,:,:))./2;
STIM.BIN.CSD.NRM = (STIM.DE.CSD.bsl(:,:,:)+STIM.NDE.CSD.bsl(:,:,:))./2;

% NRMvsDI
for i = 1:length(sDE) % sDE and sNDE are defined earlier in this section 
STIM.DI.aMUA.pc_NRM.all(i,:,:) = (STIM.DE.aMUA.pc.all(sDE(i),:,:) + STIM.NDE.aMUA.pc.all(sNDE(i),:,:))./2;
STIM.DI.aMUA.pc_NRM.layers(i,:,:) = (STIM.DE.aMUA.pc.layers(sDE(i),:,:) + STIM.NDE.aMUA.pc.layers(sNDE(i),:,:))./2;
STIM.DI.aMUA.pc_NRM.coll(i,:,:) = (STIM.DE.aMUA.pc.coll(sDE(i),:,:) + STIM.NDE.aMUA.pc.coll(sNDE(i),:,:))./2;
STIM.DI.aMUA.pc_NRM.coll_layers(i,:,:) = (STIM.DE.aMUA.pc.coll_layers(sDE(i),:,:) + STIM.NDE.aMUA.pc.coll_layers(sNDE(i),:,:))./2;
STIM.DI.CSD.NRM(i,:,:) = (STIM.DE.CSD.bsl(sDE(i),:,:) + STIM.NDE.CSD.bsl(sNDE(i),:,:))./2;
end

%% End 
if loop == false
fprintf('\nWe did it, gang.\n');
fprintf('\nNo explosion today.\n');
end

%% Plot: (SNAPSHOT)
% All averaged, baseline corrected trials (LFP, aMUA, CSD, iCSD)
pre = -50;
post = 600;
STIM.refwin = pre:post; % reference window for line plotting
STIM.channels = 1:size(STIM.aMUA.raw,2);  % how many channels (nct is a predefined variable with the exact number of channels

h1 = figure('position',[15,135,1200,500]);
clear i
avg_fields = fieldnames(STIM.avg);
for i = 1:length(avg_fields)
subplot(1,4,i)
f_ShadedLinePlotbyDepthMod((STIM.avg.(avg_fields{i})),0:(1/(numel(STIM.channels))):1,STIM.refwin, STIM.channels, 1); % this function baseline corrects and scales
hold on
plot([0 0], ylim,'k')
plot([STIM.off STIM.off], ylim,'k','linestyle','-.','linewidth',0.5)
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
climit = max(abs(get(gca,'CLim'))*.5);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
plot([STIM.off STIM.off], ylim,'k','linestyle','-.','linewidth',0.5)
title('Interpolated CSD')
xlabel('time (ms)')
clrbar = colorbar; clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts indexed down from surface');
hold off

sgtitle({'All trials triggered to stim onset',BRdatafile}, 'Interpreter', 'none');
set(h,'position',[0.7483,0.1253,0.1055,0.6826]);

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_snapshot',BRdatafile), '-jpg', '-transparent');

%% Saving Workspace to D drive
if loop == false
Prompt = ('Would you like to save the workspace? (y/n)');
str = input(Prompt,'s');
    if str == 'n' || str == 'N'
        error('Did not save workspace');
    else 
        fprintf('\nSaving workspace...\n');
    end
else 
    fprintf('\nSaving workspace...\n');
end

cd('D:\mcosinteroc_all-E\')
save(sprintf('%s',BRdatafile),'STIM','BRdatafile');

fprintf('%s.mat saved\n',BRdatafile);

%clearvars -except name loop
%end