%% BMznorm 
% script to select an input file, load in .nev,
% .ns2, and ns6, pull out stim onset and tie those timepoints to the
% raw neural data. Each stim onset with animal fixation is a trial. Generate LFP, aMUA, and CSD --triggered to a reference
% window. Average across trials. Baseline correct the data. Plot each
% contact channel for the duration of the stimulus. 
clear
%% Establish directories and set path

if strcmp(getenv('USER'),'maierav')                                      %retrieves environment variable 'USER' 
    npmkdir  = '/Users/alex 1/Desktop/LAB/Brock/OLD/NPMK-4.5.3.0/NPMK/'; %directory for Alex's machine
    nbanalysisdir   = '/Users/alex 1/Desktop/LAB/bootcamp/nbanalysis/';  %directory for Alex's machine
    datadir  = '/Users/alex 1/Desktop/LAB/';                             %directory for the stored data
else
    npmkdir  = '/users/bmitc/Documents/MATLAB/NPMK/';                    %neural processing matlab kit (NPMK)
    nbanalysisdir   = '/users/bmitc/Documents/MATLAB/nbanalysis/';       %directory with various tools for opening, loading, and processing 
    %datadir  = '/users/bmitc/Box Sync/DATA/';
    datadir = 'users/bmitc/Documents/MATLAB/data/';
end

addpath(genpath(npmkdir))
addpath(genpath(nbanalysisdir))
addpath(genpath(datadir))

BRdatafile = '160922_E_rfori001'; 
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
elseif isequal(BRdatafile, '190415_B_cinteroc002') || isequal(BRdatafile, '190321_B_cinteroc001')...
        || isequal(BRdatafile,'161003_E_cinteroc002')
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

pre   = -50; % 50ms before stim onset
post  = 300; % ~20 ms after stim offset

STIM.LFP  = trigData(LFP,STIM.onsetsdown,-pre,post); %pre variable is in absolute units 
STIM.CSD  = trigData(CSD,STIM.onsetsdown,-pre,post); 
STIM.aMUA = trigData(MUA,STIM.onsetsdown,-pre,post); 

%% Averaging across trials & baseline correct

avg.LFP = mean(STIM.LFP,3);
avg.aMUA = mean(STIM.aMUA,3);
avg.CSD = mean(STIM.CSD,3);

[bsl.LFP] = BMbasecorrect(avg.LFP);
[bsl.aMUA] = BMbasecorrect(avg.aMUA);
[bsl.CSD] = BMbasecorrect(avg.CSD);

%% Vectors of logicals for conditions of interest
contrast = unique(STIM.contrast);

clear i
for i = 1:length(contrast)
STIM.Mconditions(i,:) = STIM.contrast == contrast(i) & STIM.fixedc == 0; 
end

clear i
for i = 1:length(contrast)
STIM.Bconditions(i,:) = STIM.contrast == contrast(i) & STIM.fixedc == contrast(i); 
end

%% Data conversion for aMUA (Z score or percent change)

cMUA = nan(size(STIM.aMUA,1),size(STIM.aMUA,2),size(STIM.aMUA,3));

for t = 1:size(STIM.aMUA,3)
    for c = 1:size(STIM.aMUA,2)
%       cMUA(:,c,t) = (STIM.aMUA(:,c,t)-mean(STIM.aMUA(25:75,c,t)))./(std(STIM.aMUA(25:75,c,t))); %z score
        cMUA(:,c,t) = (STIM.aMUA(:,c,t)-mean(STIM.aMUA(25:75,c,t)))./(mean(STIM.aMUA(25:75,c,t)))*100; %percent change
    end
end

avg.cMUA = mean(cMUA,3);
bsl.cMUA = BMbasecorrect(avg.cMUA);


%% Averaged trials by condition
%  Now that I have my conditions as logicals in a structure and my aMUA in
%  z-scores, time to pull out AVG'd z-scored aMUA for each condition

clear m MzMUA BzMUA
for m = 1:size(STIM.Mconditions,1)
    Mon_cMUA(m).contrast = mean(cMUA(:,:,STIM.Mconditions(m,:)),3); %#ok<SAGROW>
    Bin_cMUA(m).contrast = mean(cMUA(:,:,STIM.Bconditions(m,:)),3); %#ok<SAGROW>
end

%% collapsing across time for each condition

clear i coll_mon
for i=1:size(cMUA,2)
    coll_mon.contrast(i,:)  = mean(cMUA(i).contrast(80:330,:),1);
end

clear i coll_bin
for i=1:size(cMUA,2)
    coll_bin.contrast(i,:)  = mean(cMUA(i).contrast(80:330,:),1);
end

%% Plotting

%% Plotting all averaged, baseline corrected trials (SNAPSHOT)

refwin = pre:post; % reference window for line plotting
channels = 1:nct;  % how many channels (nct is a predefined variable with the exact number of channels
offset = ((STIM.offsets(1)-STIM.onsets(1))/(30)); % stimulus offset as calculated from grating text file.

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
%plot([0 0], ylim,'k')
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
% h2 = figure('position',[15,135,1200,500]);
% clear i 
% for i = 1:length(contrast)
% subplot(1,length(contrast),i);
% f_ShadedLinePlotbyDepth(mean(STIM.aMUA(:,:,STIM.Mconditions(i,:)),3),1:24,refwin,channels,1)
% hold on
% plot([0 0], ylim,'k')
% plot([offset offset], ylim,'k')
% if i == 1
%     title({contrast(i),' contrast in both eyes'});
% else 
%     title({contrast(i),' contrast in DE'});
% end
% xlabel('time (ms)')
% ylabel('contacts indexed down from surface')
% hold off
% end
% 
% sgtitle({'aMUA | Varying contrast to dominant eye',BRdatafile},'Interpreter','none');
% 
% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_contrasts',BRdatafile), '-jpg', '-transparent');

h2 = figure('position',[15,135,1200,500]);
subplot(1,length(contrast),numel(contrast))
hold on
contrastValue = max(contrast);
f_ShadedLinePlotbyDepth(mean(STIM.aMUA(:,:,STIM.Mconditions(numel(contrast),:)),3),0:(1/(numel(channels))):1,refwin,channels,1);
scalingfactor = get(gca);
title('1 contrast in DE');
clear i 
for i = 1:length(contrast)-1
    subplot(1,length(contrast),i);
    useblack = false;
    f_ShadedLinePlotbyDepth_BAM(mean(STIM.aMUA(:,:,STIM.Mconditions(i,:)),3),1:24,refwin,channels,1,useblack,.3)
   
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
% 
% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_contrasts',BRdatafile), '-jpg', '-transparent');

%% Bar plot (Bar-contrasts)
format bank; rcontrast = round(contrast,2,'significant');

h4 = figure('Position', [60 211 1100 300]);
subplot(1,7,1)
bar(coll_bin.contrast(:,1),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(coll_mon.contrast(:,1),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels('')
xlabel('contrast level')
ylabel('percent change from baseline')
title('Contact 1');
hold off

subplot(1,7,2)
bar(coll_bin.contrast(:,6),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(coll_mon.contrast(:,6),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels('')
xlabel('contrast level')
%ylabel('percent change from baseline');
title('Contact 6');
hold off

subplot(1,7,3)
bar(coll_bin.contrast(:,11),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(coll_mon.contrast(:,11),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels('')
xlabel('contrast level')
%ylabel('percent change from baseline')
title('Contact 11');
%legend('Binocular stimulation','Monocoular stimulation','Location','southoutside');
hold off

subplot(1,7,4)
bar(coll_bin.contrast(:,16),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(coll_mon.contrast(:,16),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels('')
xlabel('contrast level')
%ylabel('percent change from baseline')
title('Contact 16');
hold off

subplot(1,7,5)
bar(coll_bin.contrast(:,21),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(coll_mon.contrast(:,21),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels('')
xlabel('contrast level')
%ylabel('percent change from baseline')
title('Contact 17');
hold off

subplot(1,7,6)
bar(coll_bin.contrast(:,21),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(coll_mon.contrast(:,21),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels('')
xlabel('contrast level')
%ylabel('percent change from baseline')
title('Contact 21');
hold off

subplot(1,7,7)
bar(coll_bin.contrast(:,24),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(coll_mon.contrast(:,24),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels('')
xlabel('contrast level')
%ylabel('percent change from baseline')
title('Contact 24');
hold off


%sgtitle({'Monocular vs Binocular response as a function of contrast',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_bar-contrasts',BRdatafile), '-jpg', '-transparent');


%% 3D surface plot (MESH)
format bank; rcontrast = round(contrast,2,'significant');

h7 = figure('Position', [49 111 1100 470]);
sp2 = subplot(1,3,2);
surf(contrast,fliplr(channels),coll_bin.contrast');
hold on
colormap(gca, 'jet'); % this makes the red color the sinks and the blue color the sources (convention)
cm2 = colorbar; 
set(gca,'tickdir','out'); 
climit = max(abs(get(gca,'CLim'))*1);
%climit = 1;
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
title('Binocular response')
xlabel('contrast level')
clrbar = colorbar; clrbar.Label.String = 'percent change'; 
set(clrbar.Label,'rotation',270,'fontsize',8,'VerticalAlignment','middle');
ylabel('contacts');
xticklabels(0:0.5:1);
yticklabels({'','18','8'});
% zaxis = get(gca,'zlim');
zaxis = [-10 100];
set(gca,'zlim',zaxis);
hold off

sp1 = subplot(1,3,1);
surf(contrast,fliplr(channels),coll_mon.contrast');
hold on
colormap(gca,'jet'); % this makes the red color the sinks and the blue color the sources (convention)
cm1 = colorbar; 
set(gca,'zlim',zaxis,'tickdir','out');
%climit = max(abs(get(gca,'CLim'))*1);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
%zlim = zaxis;
title('Monocular response')
xlabel('contrast level')
clrbar = colorbar; clrbar.Label.String = 'percent change'; 
set(clrbar.Label,'rotation',270,'fontsize',8,'VerticalAlignment','middle');
ylabel('contacts');
zlabel('percent change from baseline');
xticklabels(0:0.5:1);
yticklabels({'','18','8'});
hold off

sp3 = subplot(1,3,3);
surf(contrast,fliplr(channels),coll_bin.contrast'-coll_mon.contrast');
hold on
colormap(gca,'bone'); % this makes the red color the sinks and the blue color the sources (convention)
cm3 = colorbar; 
ax = set(gca,'tickdir','out'); 
climit = max(abs(get(gca,'CLim'))*1);
% climit = .4;
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
set(gca,'zlim',zaxis,'tickdir','out');
title('Monocular subtracted from binocular')
xlabel('contrast level')
clrbar = colorbar; clrbar.Label.String = 'percent difference'; 
set(clrbar.Label,'rotation',270,'fontsize',8,'VerticalAlignment','middle');
ylabel('contacts');
xticklabels(0:0.5:1);
yticklabels({'','18','8'});
zaxis = get(gca,'zlim');
set(gca,'zlim',zaxis);
hold off

sgtitle({'Monocular versus Binocular response as a function of contrast',BRdatafile},'Interpreter','none');

set(sp1, 'Position', [.05 .2 .20 .6])
set(sp2, 'Position', [.38 .2 .20 .6])
set(sp3, 'Position', [.70 .2 .20 .6])

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_surf-pc',BRdatafile), '-jpg', '-transparent');

%% 2D Surface plot (IMAGESC)

h8 = figure('Position', [49 111 1400 500]);

sp2 = subplot(1,3,2);
imagesc(1:length(contrast),channels,coll_bin.contrast');
hold on
colormap(colormap('jet')); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; 
set(gca,'tickdir','out');  
climit = max(abs(get(gca,'CLim'))*1);
% climit = 1;
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
title('Binocular contrast response');
xlabel('\fontsize{12}contrast level')
clrbar = colorbar; clrbar.Label.String = '% change'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('\fontsize{12}contacts indexed down from surface');
xticklabels(rcontrast);
hold off

sp1 = subplot(1,3,1);
imagesc(1:length(contrast),channels,coll_mon.contrast');
hold on
colormap(colormap('jet')); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; 
set(gca,'tickdir','out');  
%climit = max(abs(get(gca,'CLim'))*1);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
title('Monocular contrast response');
xlabel('\fontsize{12}contrast level')
clrbar = colorbar; clrbar.Label.String = '% change'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('\fontsize{12}contacts indexed down from surface');
xticklabels(rcontrast);
hold off

sp3 = subplot(1,3,3);
imagesc(1:length(contrast),channels,coll_bin.contrast'-coll_mon.contrast');
hold on
colormap(gca,'bone'); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; 
set(gca,'tickdir','out');  
climit = max(abs(get(gca,'CLim'))*1);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
title('Monocular subtracted from binocular');
xlabel('\fontsize{12}contrast level')
clrbar = colorbar; clrbar.Label.String = '% difference'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('\fontsize{12}contacts indexed down from surface');
xticklabels(rcontrast);
hold off

%sgtitle({'Monocular vs Binocular response as a function of contrast',BRdatafile},'Interpreter','none');

set(sp1, 'Position', [.05 .2 .20 .6])
set(sp2, 'Position', [.38 .2 .20 .6])
set(sp3, 'Position', [.70 .2 .20 .6])

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_imagesc',BRdatafile), '-jpg', '-transparent');

%% Experimental tightplot
h9 = figure('Position', [280,58,380,582]);

[ha, pos] = tight_subplot(1,6,[0.005 .03],[.1 .15],[.2 .2]); %channels, columns, [spacing], [top and bottom margin], [left and right margin]
for c = 1:length(contrast)
    
    axes(ha(c)); % ha is a variable that gets the axis of each subplot
    
    plot(channels,coll_mon.contrast(c,:),'k','linewidth',.5);
    hold on
    plot(channels,coll_bin.contrast(c,:),'.-b','linewidth',0.5);
    %fill([contrast fliplr(contrast)], [coll_mon.contrast(:,c) fliplr(coll_bin.contrast(:,c))], 'b')
    hold off
    ylim([-10 100]);
%     for cc =1:5:24
%     yticklabels(c);
%     end
    grid off
    
    if c < 6
    set(gca, 'XTick', '') %removing the units on the x axis (if i < 24) 
    set(gca, 'YTick', '') %removing the units on the x axis (if i < 24)
    p = gca; p.XAxis.Visible = 'off'; %remove the x axis
    end
    
    if c == 1
        title({'Contrast response curves across contacts',BRdatafile},'Interpreter', 'none')
    end
   
%     if c == 12
%         ylabel('\fontsize{12}contacts in order of depth');
%     end
    
    if c == 6
        xticks(1:4:24)
        xticklabels(1:4:24)
    end
end

set(ha(1:5), 'XTickLabel',''); 
set(ha(1:6), 'box', 'off');

xlabel('\fontsize{12}contacts in order of depth');
% legend('Monocular stimulus','Binocular stimulus','Location','eastoutside');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_tightplot',BRdatafile), '-jpg', '-transparent');

%% Experimental summary line plots

figure('position',[185 150 887 450]);
subplot(1,3,1)
hold on
plot(fliplr(coll_mon.contrast(:,:)),channels);
ylabel('contacts indexed down from surface');
yticks = 'on';
yticklabels({;'','20','15','10','5'});
xlim([-5 30]);
% xticks(-10 100);
%xticklabels(-5:10:40);
xlabel('Percent change');
title('Monocular');
% yticks(1:24)
% yticklabels(1:24);
hold off

subplot(1,3,2)
plot(fliplr(coll_bin.contrast(:,:)),channels);
hold on
ylabel('contacts indexed down from surface');
set(gca,'box','off');
yticklabels({;'','20','15','10','5'});
xlim([-5 30]);
xlabel('Percent change');
title('Binocular');
hold off

subplot(1,3,3)
plot(fliplr(coll_bin.contrast(:,:))-fliplr(coll_mon.contrast(:,:)),channels);
hold on 
ylabel('contacts indexed down from surface');
yticks = 'on';
yticklabels({;'','20','15','10','5'});
xlim([-5 30]);
set(gca,'box','off');
xlabel('Percent change');
title('Subtraction (bin - mon)');
hold off

% legend('0','0.05','0.14','0.37','1','Location','eastoutside','orientation','horizontal');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_lineplots',BRdatafile), '-jpg', '-transparent');

%% experimental summary lineplots 2

figure('position',[185 150 887 450]);
subplot(1,5,1)
plot(fliplr(coll_mon.contrast(1,:)),channels,'k','linewidth',.5);
hold on
plot(fliplr(coll_bin.contrast(1,:)),channels,'k','linewidth',.5);
ylabel('contacts indexed down from surface');
set(gca,'box','off');
yticklabels({;'','20','15','10','5'});
xlim([-5 30]);
xlabel('Percent change');
title('0 contrast');
hold off

subplot(1,5,2)
plot(fliplr(coll_mon.contrast(2,:)),channels,'-b','linewidth',.5);
hold on
plot(fliplr(coll_bin.contrast(2,:)),channels,'.-r','linewidth',.5);
ylabel('contacts indexed down from surface');
set(gca,'box','off');
yticklabels({;'','20','15','10','5'});
xlim([-5 30]);
xlabel('Percent change');
title('0.05 contrast');
hold off

subplot(1,5,3)
plot(fliplr(coll_mon.contrast(3,:)),channels,'-b','linewidth',.5);
hold on
plot(fliplr(coll_bin.contrast(3,:)),channels,'.-r','linewidth',.5);
ylabel('contacts indexed down from surface');
set(gca,'box','off');
yticklabels({;'','20','15','10','5'});
xlim([-5 30]);
xlabel('Percent change');
title('0.14 contrast');
hold off

subplot(1,5,4)
plot(fliplr(coll_mon.contrast(4,:)),channels,'-b','linewidth',.5);
hold on
plot(fliplr(coll_bin.contrast(4,:)),channels,'.-r','linewidth',.5);
ylabel('contacts indexed down from surface');
set(gca,'box','off');
yticklabels({;'','20','15','10','5'});
xlim([-5 30]);
xlabel('Percent change');
title('0.37 contrast');
hold off

subplot(1,5,5)
plot(fliplr(coll_mon.contrast(5,:)),channels,'-b','linewidth',.5);
hold on
plot(fliplr(coll_bin.contrast(5,:)),channels,'.-r','linewidth',.5);
ylabel('contacts indexed down from surface');
set(gca,'box','off');
yticklabels({;'','20','15','10','5'});
xlim([-5 30]);
xlabel('Percent change');
title('1 contrast');
hold off

legend('Monocular stimulation','Binocular stimulation','Location','eastoutside');

% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_lineplots2',BRdatafile), '-jpg', '-transparent');

%%
supra_mon = mean(coll_mon.contrast(:,1:10),2);
gran_mon = mean(coll_mon.contrast(:,11:17),2);
infra_mon = mean(coll_mon.contrast(:,18:24),2);

supra_bin = mean(coll_bin.contrast(:,1:10),2);
gran_bin = mean(coll_bin.contrast(:,11:17),2);
infra_bin = mean(coll_bin.contrast(:,18:24),2);

layers_MON = [supra_mon,gran_mon, infra_mon];
layers_BIN = [supra_bin,gran_bin, infra_bin];

figure('Position', [60 211 1100 300]);
subplot(1,3,1)
bar(layers_BIN(:,1),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(layers_MON(:,1),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('percent change from baseline');
title('Supragranular');
hold off

subplot(1,3,2)
bar(layers_BIN(:,2),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(layers_MON(:,2),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('percent change from baseline');
title('Granular');
hold off

subplot(1,3,3)
bar(layers_BIN(:,3),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(layers_MON(:,3),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('percent change from baseline');
title('Infragranular');
hold off

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_layers-uninformed',BRdatafile), '-jpg', '-transparent');

%% Semilog 

figure;
format bank;
subplot(2,3,1)
semilogx(layers_MON(:,1));
hold on
semilogx(layers_BIN(:,1));
xticklabels(rcontrast);
ylim([-5 30]);
title('Supragranular');
hold off

subplot(2,3,2)
semilogx(layers_MON(:,2));
hold on
semilogx(layers_BIN(:,2));
xticklabels(rcontrast);
ylim([-5 30]);
title('Granular');
hold off

subplot(2,3,3)
semilogx(layers_MON(:,3));
hold on
semilogx(layers_BIN(:,3));
xticklabels(rcontrast);
ylim([-5 30]);
title('Infragranular');
hold off

subplot(2,3,4)
semilogy(layers_MON(:,1));
hold on
semilogy(layers_BIN(:,1));
ylim([0 50]);
xlabel('contrast');
hold off

subplot(2,3,5)
semilogy(layers_MON(:,2));
hold on
semilogy(layers_BIN(:,2));
ylim([0 50]);
xlabel('contrast');
hold off

subplot(2,3,6)
semilogy(layers_MON(:,3));
hold on
semilogy(layers_BIN(:,3));
ylim([0 50]);
xlabel('contrast');
hold off

legend('Binocular stimulation','Monocoular stimulation','Location','southoutside');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_semilog',BRdatafile), '-jpg', '-transparent');