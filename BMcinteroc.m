%% BMcinteroc
% script to select an input file, load in .nev,
% .ns2, and ns6, pull out stim onset and tie those timepoints to the
% raw neural data. Each stim onset with animal fixation is a trial. Generate LFP, aMUA, and CSD --triggered to a reference
% window. Average across trials and plot each response by
% contact channel for the duration of the stimulus. 
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
    datadir  = '/users/bmitc/Box Sync/DATA/';
    %datadir = 'users/bmitc/Documents/MATLAB/data/';
end

addpath(genpath(npmkdir))
addpath(genpath(nbanalysisdir))
addpath(genpath(datadir))

BRdatafile = '160523_E_mcosinteroc002'; 
filename = [datadir BRdatafile];

%% Define layers by contact
supra = 10:13; % contact range for supragranular layer
gran = 14:19; % contact range for granular layer
infra = 20:24; % contact range for infragranular layer

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
elseif isequal(BRdatafile, '190415_B_cinteroc002') || isequal(BRdatafile, '190321_B_cinteroc001')...
        || isequal(BRdatafile,'161003_E_cinteroc002') || isequal(BRdatafile,'190210_B_cinteroc001')
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

switch sortdirection
    case 'ascending'
        MUA = MUA(:,idx);
        LFP = LFP(:,idx);
        sortedLabels = NeuralLabels(idx); 
    case 'descending'
        MUA = MUA(:,flipud(idx));
        LFP = LFP(:,flipud(idx));
        sortedLabels = NeuralLabels(flipud(idx)); 
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
avg.LFP = mean(STIM.LFP,3);
avg.aMUA = mean(STIM.aMUA,3);
avg.CSD = mean(STIM.CSD,3);

[bsl.LFP] = BMbasecorrect(avg.LFP);
[bsl.aMUA] = BMbasecorrect(avg.aMUA);
[bsl.CSD] = BMbasecorrect(avg.CSD);

%% Vectors of logicals for conditions of interest
contrast = unique(STIM.contrast);

clear i STIM.Mconditions
for i = 1:length(contrast)
STIM.Mconditions(i,:) = STIM.contrast == contrast(i) & STIM.fixedc == 0; 
end

clear i STIM.NDEconditions

for i = 1:length(contrast)
STIM.NDEconditions(i,:) = STIM.contrast == 0 & STIM.fixedc == contrast(i); 
end

clear i STIM.Bconditions
for i = 1:length(contrast)
STIM.Bconditions(i,:) = STIM.contrast == contrast(i) & STIM.fixedc == contrast(i); 
end

%% Data conversion for aMUA (Z score or percent change)

clear cMUA

cMUA = nan(size(STIM.aMUA,1),size(STIM.aMUA,2),size(STIM.aMUA,3));

clear t c
for t = 1:size(STIM.aMUA,3)
    for c = 1:size(STIM.aMUA,2)
%       cMUA(:,c,t) = (STIM.aMUA(:,c,t)-mean(STIM.aMUA(25:75,c,t)))./(std(STIM.aMUA(25:75,c,t))); %z score
        cMUA(:,c,t) = (STIM.aMUA(:,c,t)-mean(STIM.aMUA(25:75,c,t)))./(mean(STIM.aMUA(25:75,c,t)))*100; %percent change
    end
end

avg.cMUA = mean(cMUA,3);
bsl.cMUA = BMbasecorrect(avg.cMUA);

%% Averaged trials by condition

clear m Mon_cMUA NDE_cMUA Bin_cMUA 
for m = 1:size(STIM.Mconditions,1)
    Mon_cMUA(m).contrast = mean(cMUA(:,:,STIM.Mconditions(m,:)),3); %#ok<SAGROW>
    NDE_cMUA(m).contrast = mean(cMUA(:,:,STIM.NDEconditions(m,:)),3);
    Bin_cMUA(m).contrast = mean(cMUA(:,:,STIM.Bconditions(m,:)),3); %#ok<SAGROW>
end

clear m Mon_CSD BIN_CSD
for m = 1:size(STIM.Mconditions,1)
    Mon_CSD(m).contrast = mean(STIM.CSD(:,:,STIM.Mconditions(m,:)),3); %#ok<SAGROW>
    NDE_CSD(m).contrast = mean(STIM.CSD(:,:,STIM.NDEconditions(m,:)),3);
    Bin_CSD(m).contrast = mean(STIM.CSD(:,:,STIM.Bconditions(m,:)),3); %#ok<SAGROW>
end

%% collapsing across time for each condition

clear i coll_mon
for i=1:size(Mon_cMUA,2)
    coll_mon.full(i,:)  = mean(Mon_cMUA(i).contrast(80:offset,:),1);
end

clear i coll_bin
for i=1:size(Bin_cMUA,2)
    coll_bin.full(i,:)  = mean(Bin_cMUA(i).contrast(80:offset,:),1);
end

clear i coll_nde
for i=1:size(NDE_cMUA,2)
    coll_bin.full(i,:)  = mean(Bin_cMUA(i).contrast(80:offset,:),1);
end

clear i
for i=1:size(Mon_cMUA,2)
    coll_mon.transient(i,:)  = mean(Mon_cMUA(i).contrast(80:200,:),1);
end

clear i
for i=1:size(Bin_cMUA,2)
    coll_bin.transient(i,:)  = mean(Bin_cMUA(i).contrast(80:200,:),1);
end

clear i
for i=1:size(Mon_cMUA,2)
    coll_mon.sustained(i,:)  = mean(Mon_cMUA(i).contrast(201:offset,:),1);
end

clear i 
for i=1:size(Bin_cMUA,2)
    coll_bin.sustained(i,:)  = mean(Bin_cMUA(i).contrast(201:offset,:),1);
end

%% Plotting: all averaged, baseline corrected trials (SNAPSHOT)

refwin = pre:post; % reference window for line plotting
channels = 1:nct;  % how many channels (nct is a predefined variable with the exact number of channels

h1 = figure('position',[15,135,1200,500]);
clear i
avg_fields = fieldnames(avg);
for i = 1:length(avg_fields)
subplot(1,4,i)
f_ShadedLinePlotbyDepthMod((avg.(avg_fields{i})),0:(1/(numel(channels))):1,refwin, channels, 1);
hold on
plot([0 0], ylim,'k')
plot([offset offset], ylim,'k','linestyle','-.','linewidth',0.5)
title(avg_fields{i})
xlabel('time (ms)')
ylabel('contacts indexed down from surface')
hold off
end

bAVG_iCSD = filterCSD(bsl.CSD')';

h = subplot(1,4,4);
imagesc(refwin,channels,bAVG_iCSD');
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

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_snapshot',BRdatafile), '-jpg', '-transparent');

%% Varying contrast, monocular stimulation (Contrasts)

figure('position',[15,135,1200,500]);
subplot(1,length(contrast),numel(contrast))
hold on
contrastValue = max(contrast);
global scalingfactor
f_ShadedLinePlotbyDepth(mean(STIM.aMUA(:,:,STIM.Mconditions(numel(contrast),:)),3),0:(1/(numel(channels))):1,refwin,channels,1,1);
title('1 contrast in DE');
xlabel('time (ms)');
hold off

clear i 
for i = 1:length(contrast)-1
    subplot(1,length(contrast),i);
    f_ShadedLinePlotbyDepth_BAM(mean(STIM.aMUA(:,:,STIM.Mconditions(i,:)),3),0:(1/(numel(channels))):1,refwin,channels,1,1,false,scalingfactor);
   
plot([0 0], ylim,'k')
plot([offset offset], ylim,'k')
if i == 1
    title({contrast(i),' contrast in both eyes'});
else 
    title({contrast(i),' contrast in DE'});
end
xlabel('time (ms)')
ylabel('contacts indexed down from surface')
hold off
end

sgtitle({'aMUA | Varying contrast to dominant eye',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_contrasts-DE',BRdatafile), '-jpg', '-transparent');

%% Bar plot (Bar-contrasts)

figure('Position', [60 211 1100 300]);
selectchannels = [1:3:length(channels)];
%selectchannels = [15 16 17 18];
clear c
for c = 1:length(selectchannels)
    subplot(1,length(selectchannels),c)
    bar(coll_bin.full(:,selectchannels(c)),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
    hold on
    bar(coll_mon.full(:,selectchannels(c)),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
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
    'Responses collapsed across stimulus duration',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_bar-contrasts',BRdatafile), '-jpg', '-transparent');


%% Tightplot
figure('position',[185 150 887 450]);

[ha, pos] = tight_subplot(1,(numel(contrast)),[0.005 .03],[.10 .2],[.05 .05]); %channels, columns, [spacing], [bottom and top margin], [left and right margin]
clear c
for c = 1:length(contrast)
    
    axes(ha(c)); % ha is a variable that gets the axis of each subplot
    plot(fliplr(coll_mon.full(c,:)),channels,'b','linewidth',.5);
    hold on
    plot(fliplr(coll_bin.full(c,:)),channels,'.-r','linewidth',0.5);
    xlim([-10 80])
    yticklabels({flipud(1:length(channels))});
    grid off
    hold off
    
    xlabel('Percent change');
    
    if c == 1
        ylabel('contacts indexed down from surface');
        title({contrast(c),' contrast'});
    else 
        title({contrast(c),' contrast'});
    end   
   
end
sgtitle({'Monocular vs binocular MUA as a function of contrast level',BRdatafile},'Interpreter', 'none')
set(ha(1:numel(contrast)), 'box', 'off');
% legend('Monocular stimulus','Binocular stimulus','Location','eastoutside');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_tightplot',BRdatafile), '-jpg', '-transparent');

%% Summary line plots (contrast lines)

figure('position',[185 150 887 450]);

h = subplot(1,4,1);
imagesc(refwin,channels,bAVG_iCSD');
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
plot(fliplr(coll_bin.full(:,:)),channels);
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
plot(fliplr(coll_mon.full(:,:)),channels);
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
plot(fliplr(coll_bin.full(:,:))-fliplr(coll_mon.full(:,:)),channels);
hold on 
grid on
xlim([-5 20]);
set(gca,'box','off');
yticks(1:nct)
yticklabels(fliplr(1:nct))
ylim([1 nct])
xlabel('Percent change');
title('Subtraction (bin - mon)');
%legend(num2str(contrast),'Location','southoutside','orientation','horizontal');
hold off

sgtitle({'Monocular vs binocular aMUA averaged over stimulus duration',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_lineplots',BRdatafile), '-jpg', '-transparent');

%% Binning contacts into V1 layers

layers_MON.full = [mean(coll_mon.full(:,supra),2),mean(coll_mon.full(:,gran),2),mean(coll_mon.full(:,infra),2)];
layers_MON.transient = [mean(coll_mon.transient(:,supra),2),mean(coll_mon.transient(:,gran),2),mean(coll_mon.transient(:,infra),2)];
layers_MON.sustained = [mean(coll_mon.sustained(:,supra),2),mean(coll_mon.sustained(:,gran),2),mean(coll_mon.sustained(:,infra),2)];

layers_BIN.full = [mean(coll_bin.full(:,supra),2),mean(coll_bin.full(:,gran),2), mean(coll_bin.full(:,infra),2)];
layers_BIN.transient = [mean(coll_bin.transient(:,supra),2),mean(coll_bin.transient(:,gran),2),mean(coll_bin.transient(:,infra),2)];
layers_BIN.sustained = [mean(coll_bin.sustained(:,supra),2),mean(coll_bin.sustained(:,gran),2),mean(coll_bin.sustained(:,infra),2)];


%% Contrast lines across time (Timeplots)

figure('position',[213.6666666666667,149.6666666666667,724.6666666666666,425.3333333333334]);
subplot(2,3,1)
clear c
for c = 1:length(contrast)
plot(refwin,mean(Mon_cMUA(c).contrast(:,supra),2),'color','b')
hold on
plot(refwin,mean(Bin_cMUA(c).contrast(:,supra),2),'color','r')
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
for c = 1:length(contrast)
plot(refwin,mean(Mon_cMUA(c).contrast(:,gran),2),'color','b')
hold on
plot(refwin,mean(Bin_cMUA(c).contrast(:,gran),2),'color','r')
set(gca,'ylim',[-10 ylimit],'Box','off','TickDir','out')
end

xlabel('time (ms)');
title('Granular');

subplot(2,3,3)
clear c
for c = 1:length(contrast)
plot(refwin,mean(Mon_cMUA(c).contrast(:,infra),2),'color','b')
hold on
plot(refwin,mean(Bin_cMUA(c).contrast(:,infra),2),'color','r')
set(gca,'ylim',[-10 ylimit],'Box','off','TickDir','out')
end

xlabel('time (ms)');
title('Infragranular');

subplot(2,3,4)
clear c
for c = 1:length(contrast)
plot(refwin,smooth(mean(Bin_cMUA(c).contrast(:,supra),2)-(mean(Mon_cMUA(c).contrast(:,supra),2)),.1))
hold on
ylimit = max(abs(get(gcf,'ylim')));
set(gca,'ylim',[-10 ylimit/3],'Box','off','TickDir','out')
end
ylabel({'Percent difference'...
    '(bin - mon)'});
xlabel('time (ms)');

subplot(2,3,5)
clear c
for c = 1:length(contrast)
plot(refwin,smooth(mean(Bin_cMUA(c).contrast(:,gran),2)-(mean(Mon_cMUA(c).contrast(:,gran),2)),.1))
hold on
ylimit = max(abs(get(gcf,'ylim')));
set(gca,'ylim',[-10 ylimit/3],'Box','off','TickDir','out')
end

xlabel('time (ms)');

subplot(2,3,6)
clear c
for c = 1:length(contrast)
plot(refwin,smooth(mean(Bin_cMUA(c).contrast(:,infra),2)-(mean(Mon_cMUA(c).contrast(:,infra),2)),.1))
hold on
ylimit = max(abs(get(gcf,'ylim')));
set(gca,'ylim',[-10 ylimit/3],'Box','off','TickDir','out')
end

xlabel('time (ms)');

sgtitle({'V1 laminar contrast response profiles: monocular (blue) vs binocular (red)'...
   ,BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_timeplots',BRdatafile), '-jpg', '-transparent');

%% Bar Graphs of Binned Layers (Binned Layers)

rcontrast = round(contrast,2,'significant');

figure('Position', [148,73,633,487]);
subplot(3,3,1)
bar(layers_BIN.full(:,1),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(layers_MON.full(:,1),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 80]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('percent change');
title('Supragranular');
hold off

subplot(3,3,2)
bar(layers_BIN.full(:,2),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(layers_MON.full(:,2),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 80]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('percent change');
title('Granular');
hold off

subplot(3,3,3)
bar(layers_BIN.full(:,3),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(layers_MON.full(:,3),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 80]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('percent change');
title('Infragranular');
hold off

subplot(3,3,4)
bar((layers_BIN.full(:,1)-layers_MON.full(:,1)),0.8,'FaceColor',[0.20, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-10 30]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('percent difference');
%title('Supragranular');
hold off

subplot(3,3,5)
bar((layers_BIN.full(:,2)-layers_MON.full(:,2)),0.8,'FaceColor',[0.20, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-10 30]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('percent difference');
%title('Granular');
hold off

subplot(3,3,6)
bar((layers_BIN.full(:,3)-layers_MON.full(:,3)),0.8,'FaceColor',[0.20, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-10 30]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('percent difference');
hold off


subplot(3,3,7)
bar(((layers_BIN.full(:,1)-layers_MON.full(:,1))./(layers_MON.full(:,1))),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
plot(((layers_BIN.full(:,1)-layers_MON.full(:,1))./(layers_MON.full(:,1))),'-o','LineWidth',0.8);
set(gca,'box','off');
ylim([-1 2]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('fold change');
%title('Supragranular');
hold off

subplot(3,3,8)
bar(((layers_BIN.full(:,2)-layers_MON.full(:,2))./(layers_MON.full(:,2))),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
plot(((layers_BIN.full(:,2)-layers_MON.full(:,2))./(layers_MON.full(:,2))),'-o','LineWidth',0.8);
set(gca,'box','off');
ylim([-1 2]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('fold change');
%title('Granular');
hold off

subplot(3,3,9)
bar(((layers_BIN.full(:,3)-layers_MON.full(:,3))./(layers_MON.full(:,3))),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
plot(((layers_BIN.full(:,3)-layers_MON.full(:,3))./(layers_MON.full(:,3))),'-o','LineWidth',0.8);
set(gca,'box','off');
ylim([-1 2]);
xticklabels(rcontrast)
xlabel('contrast')
ylabel('fold change');
%title('Infragranular');
hold off

sgtitle({'Binned contacts by layer | aMUA responses',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_binned-layers',BRdatafile), '-jpg', '-transparent');

%% Transient vs sustained, Fold change Semilogx

trans = ((layers_BIN.transient(:,:)-layers_MON.transient(:,:))./(layers_MON.transient(:,:)));
sust = ((layers_BIN.sustained(:,:)-layers_MON.sustained(:,:))./(layers_MON.sustained(:,:)));
full = ((layers_BIN.full(:,:)-layers_MON.full(:,:))./(layers_MON.full(:,:)));

figure('position',[360,450.3,560,167.6]);
subplot(1, 3, 1)
format bank;
semilogx(contrast,trans(:,1),'-.k');
hold on
semilogx(contrast,sust(:,1),'k');
semilogx(contrast,full(:,1),'-o');
xticklabels('');
ylim([-2 2]);
xlabel('contrast level')
ylabel('Fold change');
title('Supragranular');
hold off

subplot(1, 3, 2)
format bank;
semilogx(contrast,trans(:,2),'-.k');
hold on
semilogx(contrast,sust(:,2),'k');
semilogx(contrast,full(:,2),'-o');
xticklabels('');
xlabel('contrast level')
ylim([-2 2]);
title('Granular');
hold off

subplot(1, 3, 3)
format bank;
semilogx(contrast,trans(:,3), '-.k');
hold on
semilogx(contrast,sust(:,3),'k');
semilogx(contrast,full(:,3),'-o');
xticklabels('');
ylim([-2 2]);
xlabel('contrast level')
title('Infragranular');
hold off
%legend('transient','sustained','full','Location','southoutside','orientation','vertical');

sgtitle('Fold change from monocular to binocular: transient vs sustained');

% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_layers-transvsust_2',BRdatafile), '-jpg', '-transparent');

%% Semilogx and Semilogy


%% Saving Workspace to D drive

Prompt = ('Would you like to save the workspace? (y/n)');
str = input(Prompt,'s');
if str == 'n' || str == 'N'
    error('Did not save workspace');
else 
    fprintf('\nSaving workspace...\n');
end

cd('D:\')
save(sprintf('%s',BRdatafile));

fprintf('Workspace saved');