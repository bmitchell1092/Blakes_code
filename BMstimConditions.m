%% mainframe rewrite
%%script to select an input file amongst a decision tree, load in .nev,
%%.ns2, and ns6, pull out stim onset and tie those timepoints to the
%%raw neural data. Each stim onset with animal fixation is a trial. Generate LFP, aMUA, and CSD --triggered to a reference
%%window. Average across trials. Baseline correct the data. Plot each
%%contact channel from -50 to 300ms relative to stim onset. 
clear
%% Establish directories and set path

if strcmp(getenv('USER'),'maierav')                                      %retrieves environment variable 'USER' 
    npmkdir  = '/Users/alex 1/Desktop/LAB/Brock/OLD/NPMK-4.5.3.0/NPMK/'; %directory for Alex's machine
    nbanalysisdir   = '/Users/alex 1/Desktop/LAB/bootcamp/nbanalysis/';  %directory for Alex's machine
    datadir  = '/Users/alex 1/Desktop/LAB/';                             %directory for the stored data
else
    npmkdir  = '/users/bmitc/Documents/MATLAB/NPMK/';                    %neural processing matlab kit (NPMK)
    nbanalysisdir   = '/users/bmitc/Documents/MATLAB/nbanalysis/';       %directory with various tools for opening, loading, and processing 
    datadir  = '/users/bmitc/Documents/MATLAB/data/';
end

addpath(genpath(npmkdir))
addpath(genpath(nbanalysisdir))
addpath(genpath(datadir))

BRdatafile = '170421_I_cinteroc003'; 
filename = [datadir BRdatafile];

%% Define stimulus patterns and select from among them

patterns   = {'rforidrft','rfsfdrft','posdisparitydrft','disparitydrft','cinterocdrft','coneinterocdrft','conedrft', ...
                'colorflicker','bwflicker','rfori','rfsize','cinteroc','color','rfsf','mcosinteroc','dotmapping'}; 
for p = 1:length(patterns) 
    
    pattern = patterns{p};    %pattern is individual elements of the above array
    
    if any(strfind(BRdatafile,pattern))         %if BRdatafile contains any string the same as pattern
        startlog = strfind(BRdatafile,pattern); %create a variable to store the number of those patches
        if ~isequal(BRdatafile(startlog:end-3),pattern), continue
        else
        match = patterns{p};
        end
    end
    
end

if isequal(match,'dotmapping')
    ext = '.gDotsXY_di';
elseif isequal(match,'161003_E_cinteroc')
    ext = ['.g' upper(match) 'DRFTGrating_di'];
else
    ext = ['.g' upper(match) 'Grating_di'];
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
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;          %we don't know why we subtract 128
EventSamples    = NEV.Data.SerialDigitalIO.TimeStamp;                   %Events in samples 
EventTimes      = floor(NEV.Data.SerialDigitalIO.TimeStampSec.*1000);   %floor rounds to nearest integer and then convert event to ms 
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSamples);            %sorts codes, samps or times into trials

%The following is a structure that only contains trials where the animal
%did not break fixation, and includes grating info and stimulus onsets. 

STIM            = sortStimandTimeData(grating,pEvC,pEvT,'stim'); 
STIM.onsetsdown = floor(STIM.onsets./30);

%% Load LFP with NS2 file
clear ext
ext = 'ns2';

%first retrieve sort direction and bank info
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
else 
    ebank = 'eB';
end

% load LFP with NS2 file

lfp = getLFP(filename,ext,ebank,sortdirection);

%% Create a new STIM.field for separate stim orientations

% {clear i count1 count2
% count1 = 0; count2 = 0;
% for i = 1:size(STIM.tilt,1)
%     if STIM.tilt(i) == 157.5 && STIM.eye(i) == 3
%         count1 = count1+1;
%         STIM.ori1onsets(count1,:) = STIM.onsetsdown(i);
%     elseif STIM.tilt(i) == 157.5 && STIM.eye(i) == 2
%         count2 = count2+1;
%         STIM.ori2onsets(count2,:) = STIM.onsetsdown(i);
%     end    
% end


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

switch sortdirection
    case 'ascending'
        MUA = MUA(:,idx);
        LFP = LFP(:,idx);
        sortedLabels = NeuralLabels(idx); 
    case 'descending'
        MUA = MUA(:,fliplr(idx));
        LFP = LFP(:,fliplr(idx));
        sortedLabels = NeuralLabels(fliplr(idx)); 
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

pre   = -50;
post  = 300;
%post  = ((STIM.offsets(1)-STIM.onsets(1))/30)+50; 

STIM.LFP  = trigData(LFP,STIM.onsetsdown,-pre,post); %pre variable is in absolute units 
STIM.CSD  = trigData(CSD,STIM.onsetsdown,-pre,post); 
STIM.aMUA = trigData(MUA,STIM.onsetsdown,-pre,post); 

%% Averaging across trials

data.STIM = fieldnames(STIM);
for avtr=1:numel(data.STIM)
    AVG.(data.STIM{avtr})  = mean(STIM.(data.STIM{avtr}),3);
end

%% Baseline correct & Plot! 
% LFP, aMUA, CSD, & iCSD

[bAVG_LFP] = BMbasecorrect(AVG.LFP);
[bAVG_aMUA] = BMbasecorrect(AVG.aMUA);
[bAVG_CSD] = BMbasecorrect(AVG.CSD);

refwin = pre:post;
channels = 1:nct;
fauxcd = 0:(1/24):1; %used in the f_shadedlineplotbydepth function as 'unknown' cortical depth

h1 = figure('position',[15,135,1200,500]);
subplot(1,4,1)
f_ShadedLinePlotbyDepth(AVG.LFP,[],refwin,channels,1.2)
hold on
plot([0 0], ylim,'k')
plot([1200 1200], ylim,'k')
title('LFP')
xlabel('time (ms)')
ylabel('contacts indexed down from surface')
hold off

subplot(1,4,2)
f_ShadedLinePlotbyDepth(AVG.aMUA,[],refwin,channels,1.2)
hold on
plot([0 0], ylim,'k')
plot([1200 1200], ylim,'k')
title('aMUA')
xlabel('time (ms)')
ylabel('contacts indexed down from surface')
hold off

subplot(1,4,3)
f_ShadedLinePlotbyDepth(AVG.CSD,[],refwin,channels,1.2)
hold on
plot([0 0], ylim,'k')
plot([((STIM.offsets(1)-STIM.onsets(1))/30)], ylim,'k')
title('CSD')
xlabel('time (ms)')
ylabel('contacts indexed down from surface')
hold off

% filtered CSD

bAVG_iCSD = filterCSD(AVG.CSD')';

subplot(1,4,4)
imagesc(refwin,channels,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(0); set(v,'color','k','linestyle','-.','linewidth',1);
set(gca,'tickdir','out');  
climit = max(abs(get(gca,'CLim'))*1);
set(gca,'CLim',[-700 700],'Box','off','TickDir','out')
%plot([0 0], ylim,'k')
plot([1200 1200], ylim,'k')
title('Interpolated CSD')
xlabel('time (ms)')
clrbar = colorbar; clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts indexed down from surface');
hold off

sgtitle({'All trials triggered to stim onset',BRdatafile}, 'Interpreter', 'none');

%export_fig('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures', '-png');

%% Plotting different conditions

%% Vector of logicals for conditions of interest
%f = fields(STIM.conditions) = {'contrast0','contrast1','contrast2'};
cont_0 = STIM.contrast == 0 & STIM.fixedc == 0;
cont_1 = STIM.contrast == 0.1060 & STIM.fixedc == 0;
cont_2 = STIM.contrast == 0.2240 & STIM.fixedc == 0;
cont_5 = STIM.contrast == 0.4730 & STIM.fixedc == 0;

cont_1b = STIM.contrast == 0.1060 & STIM.fixedc == 0.1060;
cont_2b = STIM.contrast == 0.2240 & STIM.fixedc == 0.2240;
cont_5b = STIM.contrast == 0.4730 & STIM.fixedc == 0.4730;

% get average across logical-selected trials.

cont0 = mean(STIM.aMUA(:,:,cont_0),3);

cont1 = mean(STIM.aMUA(:,:,cont_1),3);
cont2 = mean(STIM.aMUA(:,:,cont_2),3);
cont5 = mean(STIM.aMUA(:,:,cont_5),3);

cont1b = mean(STIM.aMUA(:,:,cont_1b),3);
cont2b = mean(STIM.aMUA(:,:,cont_2b),3);
cont5b = mean(STIM.aMUA(:,:,cont_5b),3);

% baseline correct
bc_cont0 = BMbasecorrect(cont0);

bc_cont1 = BMbasecorrect(cont1);
bc_cont2 = BMbasecorrect(cont2);
bc_cont5 = BMbasecorrect(cont5);

bc_cont1b = BMbasecorrect(cont1b);
bc_cont2b = BMbasecorrect(cont2b);
bc_cont5b = BMbasecorrect(cont5b);

%plot all contacts over time

h2 = figure('position',[15,135,1200,500]);
subplot(1,4,1);
f_ShadedLinePlotbyDepth(cont0,1:24,refwin,1:24,1)
hold on
plot([0 0], ylim,'k')
plot([1200 1200], ylim,'k')
%ylim([-2 2]);
title('0 contrast in both eyes');
xlabel('time (ms)')
ylabel('contacts indexed down from surface')
hold off

subplot(1,4,2);
f_ShadedLinePlotbyDepth(cont1,1:24,refwin,1:24,1)
hold on
plot([0 0], ylim,'k')
plot([1200 1200], ylim,'k')
title('0.1 contrast in DE');
xlabel('time (ms)')
ylabel('contacts indexed down from surface')
hold off

subplot(1,4,3);
f_ShadedLinePlotbyDepth(cont2,1:24,refwin,1:24,1)
hold on
plot([0 0], ylim,'k')
plot([1200 1200], ylim,'k')
title('0.2 contrast in DE');
xlabel('time (ms)')
ylabel('contacts indexed down from surface')
hold off

subplot(1,4,4);
f_ShadedLinePlotbyDepth(cont5,1:24,refwin,1:24,1)
hold on
plot([0 0], ylim,'k')
plot([1200 1200], ylim,'k')
title('0.5 contrast in DE');
xlabel('time (ms)')
ylabel('contacts indexed down from surface')
hold off

sgtitle({'aMUA | Varying contrast to dominant eye | ',BRdatafile},'Interpreter','none');
%% Collapsing across time to generate a single value for each contrast level

c_cont0 = mean(bc_cont0(pre+90:post,:),1); %c stands for collapsed; 0 contrast in both eyes
% monocular stimulus
c_cont1 = mean(bc_cont1(pre+90:post,:),1); %0.1 contrast in DE, 0 in nonDE
c_cont2 = mean(bc_cont2(pre+90:post,:),1); %0.2 contrast in DE, 0 in nonDE
c_cont5 = mean(bc_cont5(pre+90:post,:),1); %0.5 contrast in DE, 0 in nonDE
% binocular stimuli
c_cont1b = mean(bc_cont1b(pre+90:post,:),1); %0.1 contrast in both eyes
c_cont2b = mean(bc_cont2b(pre+90:post,:),1); %0.2 contrast in both eyes
c_cont5b = mean(bc_cont5b(pre+90:post,:),1); %0.5 contrast in both eyes

E = [c_cont0;c_cont1;c_cont2;c_cont5]; % Monocular: 4x1 vector of avg'd values (1 value for each contrast lvl)
B = [c_cont0;c_cont1b;c_cont2b;c_cont5b]; % Binocular: 4x1 vector of avg'd values (1 value for each contrast lvl)
C = [E;B]; %combined into one matrix
contrasts = [0 .1 .2 .5];
%% bar plots

figure('Position', [60 211 1125 320]);
subplot(1,3,1);
%superbar(E(:,5),'BarFaceColor',[.75 .75 .75],'BarEdgeColor','k','BarLineWidth',0.8);
bar(E(:,5),'FaceColor',[.75 .75 .75],'EdgeColor','k','LineWidth',0.8);
title('Contact 5')
set(gca,'box','off');
ylim([0 1.5]);
xlabel('Contrast in Dominant Eye')
ylabel('AVG multi-unit activity (uV)');
xticklabels({'0','0.1','0.2','0.5'})

subplot(1,3,2);
%superbar(E(1:4,15),'BarFaceColor',[.75 .75 .75],'BarEdgeColor','k','BarLineWidth',0.8)
bar(E(:,15),'FaceColor',[.75 .75 .75],'EdgeColor','k','LineWidth',0.8);
ylim([0 1.5]);
title('Contact 15')
set(gca,'box','off');
xlabel('Contrast in Dominant Eye')
ylabel('AVG multi-unit activity (uV)');
xticklabels({'0','0.1','0.2','0.5'})

subplot(1,3,3);
%superbar(E(:,5),'BarFaceColor',[.75 .75 .75],'BarEdgeColor','k','BarLineWidth',0.8);
bar(E(:,20),'FaceColor',[.75 .75 .75],'EdgeColor','k','LineWidth',0.8);
title('Contact 20')
set(gca,'box','off');
ylim([0 1.5]);
xlabel('Contrast in Dominant Eye')
ylabel('AVG multi-unit activity (uV)');
xticklabels({'0','0.1','0.2','0.5'})

sgtitle({'aMUA averaged across time (40ms-stim offset) for each contact',BRdatafile},'Interpreter','none');
%% normalization 
deresp = (E(:,5));
biresp = (B(:,5));
bothresp = (C(:,5));
mn                = min(bothresp); 
mx                = max(bothresp); 
normde       = (deresp - mn)/(mx - mn); 
normbi        = (biresp - mn)/(mx - mn);

h5 = figure('Position', [60 211 1125 376.6]);
subplot(1,3,1);
plot(contrasts, normde,'-o', 'color','k');
hold on
plot(contrasts, normbi,'-o', 'color','b');
title('Contact 5')
ylim([0 1]);
xlabel('Contrast level')
ylabel('Normalized contrast response');
set(gca,'box','off');
legend('Monocular stimulus','Binocular stimulus','Location','SouthEast');
hold off

% contact 15

deresp = (E(:,15));
biresp = (B(:,15));
bothresp = (C(:,15));
mn                = min(bothresp); 
mx                = max(bothresp); 
normde       = (deresp - mn)/(mx - mn); 
normbi        = (biresp - mn)/(mx - mn);

subplot(1,3,2);
plot(contrasts, normde,'-o', 'color','k');
hold on
plot(contrasts, normbi,'-o', 'color','b');
title('Contact 15')
ylim([0 1]);
xlabel('Contrast level')
ylabel('Normalized contrast response');
set(gca,'box','off');
legend('Monocular stimulus','Binocular stimulus','Location','SouthEast');
hold off

% contact 20

deresp = (E(:,20));
biresp = (B(:,20));
bothresp = (C(:,20));
mn                = min(bothresp); 
mx                = max(bothresp); 
normde       = (deresp - mn)/(mx - mn); 
normbi        = (biresp - mn)/(mx - mn);

subplot(1,3,3);
plot(contrasts, normde,'-o', 'color','k');
hold on
plot(contrasts, normbi,'-o', 'color','b');
title('Contact 20')
ylim([0 1]);
xlabel('Contrast level')
ylabel('Normalized contrast response');
set(gca,'box','off');
legend('Monocular stimulus','Binocular stimulus','Location','SouthEast');
hold off

sgtitle({'Monocular vs binocular stimulation as a function of contrast',BRdatafile},'Interpreter','none');

%% Normalization with z-scores

normdata = nan(size(STIM.aMUA,1),size(STIM.aMUA,2),size(STIM.aMUA,3));

for t = 1:size(STIM.aMUA,3)
    for c = 1:size(STIM.aMUA,2)
        normdata(:,c,t) = (STIM.aMUA(:,c,t)-mean(STIM.aMUA(25:75,c,t))/std(STIM.aMUA(25:75,c,t)));
    end
end

AVG_normdata = mean(normdata,3);
AVG_normdata = mean(AVG_normdata(90:1200,:),1);



%% testing CRF example code
%{
Gr = 100; %multiplicative response gain factor (=highest response amplitude)
Gc = 50;  %normalization pool (determines x position)
q = 10;   %exponent that determines rise and saturation
s = 0;    %exponent that determines rise and saturation
b = 0;    %baseline offset w/o stim
clear Rc
clear mRc
for c=1:100
   %compute response Rc for each conrast level c:
   Rc(c) = Gr*[c^(q+s)]/[c^q + Gc^q]+b;
   %compute multiplicative gain prediction
   Grm = 80;
   mRc(c) = Grm*[c^(q+s)]/[c^q + Gc^q]+b;
   %compute contrast-gain prediction
   cGc = 55;
   cRc(c) = Gr*[c^(q+s)]/[c^q + cGc^q]+b;
   %compute additive gain prediction
   ab = -10;
   aRc(c) = Gr*[c^(q+s)]/[c^q + Gc^q]+ab;
end
figure(1),clf
plot([1:100],Rc,'k')
hold on
plot([1:100],mRc,'r')
plot([1:100],cRc,'b')
plot([1:100],aRc,'g')
ylim([0 100])
legend('monocular CRF','multiplicative gain','contrast gain',...
   'additive gain','Location','NorthWest')
%% Testing Kacie's code for CRF curve fitting
id = 5;

clear DE_example BIN_example 
DE_example = E(:,5); %not sure if this is the right input for the runCRFFit function below
BIN_example = B(:,5); %^
unqc = contrasts;
colors = parula(10); 

[de_prd,de_xprd]   = runCRFFit(DE_example',unqc); 
[bin_prd,bin_xprd] = runCRFFit(BIN_example',unqc);
 
 
figure, 
plot(de_xprd,de_prd,'Color',colors(1,:)); 
hold on; 
plot(bin_xprd,bin_prd,'Color',colors(5,:)); 
hold on; 
plot(unqc.*100',DE_example,'o','Color',colors(1,:)); 
hold on; 
plot(unqc.*100',BIN_example,'o','Color',colors(5,:)); 
set(gca,'xscale','log','box','off','tickdir','out','linewidth',1.5); 
title(gca,num2str(id)); 

%}