%% BMcinteroc
% script to select an input file, load in .nev,
% .ns2, and ns6, pull out stim onset and tie those timepoints to the
% raw neural data. Each stim onset with animal fixation is a trial. Generate LFP, aMUA, and CSD --triggered to a reference
% window. Average across trials and plot each response by
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

BRdatafile = '161003_E_cinteroc002'; 
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

offset = round(((STIM.offsets(1)-STIM.onsets(1))/(30)),0); % stimulus offset as calculated from grating text file.
pre   = -50; % 50ms before stim onset
post  = (round(10^-2*offset)/10^-2)+100; % ~100 ms after stim offset


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

for i = 1:length(contrast)
STIM.NDEconditions(i,:) = STIM.contrast == 0 & STIM.fixedc == contrast(i); 
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

clear m Mon_cMUA Bin_cMUA
for m = 1:size(STIM.Mconditions,1)
    Mon_cMUA(m).contrast = mean(cMUA(:,:,STIM.Mconditions(m,:)),3); %#ok<SAGROW>
    Bin_cMUA(m).contrast = mean(cMUA(:,:,STIM.Bconditions(m,:)),3); %#ok<SAGROW>
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

clear i coll_mon.transient
for i=1:size(Mon_cMUA,2)
    coll_mon.transient(i,:)  = mean(Mon_cMUA(i).contrast(80:200,:),1);
end

clear i coll_bin.transient
for i=1:size(Bin_cMUA,2)
    coll_bin.transient(i,:)  = mean(Bin_cMUA(i).contrast(80:200,:),1);
end

clear i coll_mon.sustained
for i=1:size(Mon_cMUA,2)
    coll_mon.sustained(i,:)  = mean(Mon_cMUA(i).contrast(201:offset,:),1);
end

clear i coll_bin.sustained
for i=1:size(Bin_cMUA,2)
    coll_bin.sustained(i,:)  = mean(Bin_cMUA(i).contrast(201:offset,:),1);
end
%% Plotting

%% Plotting all averaged, baseline corrected trials (SNAPSHOT)

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

subplot(1,4,4);
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
export_fig(sprintf('%s_contrasts',BRdatafile), '-jpg', '-transparent');

%% Bar plot (Bar-contrasts)

figure('Position', [60 211 1100 300]);
selectchannels = [3 6 9 12 15 18 21 24];
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

sgtitle({'Monocular vs Binocular MUA by contact as a function of contrast'...
    'Responses collapsed across full stimulus duration',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_bar-contrasts-sustained',BRdatafile), '-jpg', '-transparent');


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
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
title('Binocular response')
xlabel('contrast level')
clrbar = colorbar; clrbar.Label.String = 'percent change'; 
set(clrbar.Label,'rotation',270,'fontsize',8,'VerticalAlignment','middle');
ylabel('contacts');
xticklabels(0:0.5:1);
yticklabels({'','18','8'});
zaxis = [-10 100];
set(gca,'zlim',zaxis);
hold off

sp1 = subplot(1,3,1);
surf(contrast,fliplr(channels),coll_mon.contrast');
hold on
colormap(gca,'jet'); % this makes the red color the sinks and the blue color the sources (convention)
cm1 = colorbar; 
set(gca,'zlim',zaxis,'tickdir','out');
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
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

sgtitle({'Monocular versus Binocular MUA as a function of contrast',BRdatafile},'Interpreter','none');

set(sp1, 'Position', [.05 .2 .20 .6])
set(sp2, 'Position', [.38 .2 .20 .6])
set(sp3, 'Position', [.70 .2 .20 .6])

% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_surf-pc',BRdatafile), '-jpg', '-transparent');

%% 2D Surface plot (IMAGESC)

h8 = figure('Position', [49 111 1400 500]);

sp2 = subplot(1,3,2);
imagesc(1:length(contrast),channels,coll_bin.full');
hold on
colormap(colormap('jet')); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; 
set(gca,'tickdir','out');  
climit = max(abs(get(gca,'CLim'))*1);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
title('Binocular contrast response');
xlabel('\fontsize{12}contrast level')
clrbar = colorbar; clrbar.Label.String = '% change'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('\fontsize{12}contacts indexed down from surface');
xticklabels(rcontrast);
hold off

sp1 = subplot(1,3,1);
imagesc(1:length(contrast),channels,coll_mon.full');
hold on
colormap(colormap('jet')); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; 
set(gca,'tickdir','out');  
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
title('Monocular contrast response');
xlabel('\fontsize{12}contrast level')
clrbar = colorbar; clrbar.Label.String = '% change'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('\fontsize{12}contacts indexed down from surface');
xticklabels(rcontrast);
hold off

sp3 = subplot(1,3,3);
imagesc(1:length(contrast),channels,coll_bin.full'-coll_mon.full');
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

% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('%s_imagesc',BRdatafile), '-jpg', '-transparent');

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
    yticklabels({;'','20','15','10','5'});
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
%set(ha(2:numel(contrast)), 'YTickLabel',''); 
set(ha(1:numel(contrast)), 'box', 'off');
% legend('Monocular stimulus','Binocular stimulus','Location','eastoutside');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_tightplot',BRdatafile), '-jpg', '-transparent');

%% Summary line plots (contrast lines)

figure('position',[185 150 887 450]);

subplot(1,3,2)
plot(fliplr(coll_bin.full(:,:)),channels);
hold on
ylabel('contacts indexed down from surface');
set(gca,'box','off');
yticklabels({;'','20','15','10','5'});
xlabel('Percent change');
climit = max(abs(get(gca,'xlim'))*1);
xlim([-5 climit*1.2]);
title('Binocular');
hold off

subplot(1,3,1)
hold on
plot(fliplr(coll_mon.full(:,:)),channels);
ylabel('contacts indexed down from surface');
yticklabels({;'','20','15','10','5'});
xlim([-5 climit*1.2]);
xlabel('Percent change');
title('Monocular');
hold off

subplot(1,3,3)
plot(fliplr(coll_bin.full(:,:))-fliplr(coll_mon.full(:,:)),channels);
hold on 
ylabel('contacts indexed down from surface');
yticklabels({;'','20','15','10','5'});
xlim([-5 climit*1.2]);
set(gca,'box','off');
xlabel('Percent change');
title('Subtraction (bin - mon)');
hold off

sgtitle({'Monocular vs binocular aMUA as a function of contrast',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_lineplots',BRdatafile), '-jpg', '-transparent');

%% Percent difference from monocular to binocular stimulation

%% Binning contacts into V1 layers

supra = 1:10; % contact range for supragranular layer
gran = 11:16; % contact range for granular layer
infra = 17:24; % contact range for infragranular layer

layers_MON.full = [mean(coll_mon.full(:,supra),2),mean(coll_mon.full(:,gran),2),mean(coll_mon.full(:,infra),2)];
layers_MON.transient = [mean(coll_mon.transient(:,supra),2),mean(coll_mon.transient(:,gran),2),mean(coll_mon.transient(:,infra),2)];
layers_MON.sustained = [mean(coll_mon.sustained(:,supra),2),mean(coll_mon.sustained(:,gran),2),mean(coll_mon.sustained(:,infra),2)];

layers_BIN.full = [mean(coll_bin.full(:,supra),2),mean(coll_bin.full(:,gran),2), mean(coll_bin.full(:,infra),2)];
layers_BIN.transient = [mean(coll_bin.transient(:,supra),2),mean(coll_bin.transient(:,gran),2),mean(coll_bin.transient(:,infra),2)];
layers_BIN.sustained = [mean(coll_bin.sustained(:,supra),2),mean(coll_bin.sustained(:,gran),2),mean(coll_bin.sustained(:,infra),2)];

%% Bar Graphs of Binned Layers (Layers Informed)
figure('Position', [60 211 1100 300]);
subplot(1,3,1)
bar(layers_BIN.full(:,1),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(layers_MON.full(:,1),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('percent change from baseline');
title('Supragranular');
hold off

subplot(1,3,2)
bar(layers_BIN.full(:,2),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(layers_MON.full(:,2),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('percent change from baseline');
title('Granular');
hold off

subplot(1,3,3)
bar(layers_BIN.full(:,3),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(layers_MON.full(:,3),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-10 100]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('percent change from baseline');
title('Infragranular');
hold off

sgtitle({'Binned contacts by layer | aMUA responses',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_layers-informed',BRdatafile), '-jpg', '-transparent');

%% Binned by layer, relative change from Monocular to Binocular

figure('Position', [60 211 1100 300]);
subplot(1,3,1)
bar(((layers_BIN.full(:,1)-layers_MON.full(:,1))./(layers_MON.full(:,1))),0.8,'FaceColor',[0.20, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-0.2 4]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('fold change from monocular MUA');
title('Supragranular');
hold off

subplot(1,3,2)
bar(((layers_BIN.full(:,2)-layers_MON.full(:,2))./(layers_MON.full(:,2))),0.8,'FaceColor',[0.20, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-0.2 4]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('fold change from monocular MUA');
title('Granular');
hold off

subplot(1,3,3)
bar(((layers_BIN.full(:,3)-layers_MON.full(:,3))./(layers_MON.full(:,3))),0.8,'FaceColor',[0.20, 0.2, 0.2],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-0.2 4]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('fold change from monocular MUA');
title('Infragranular');
hold off

%sgtitle({'Binned contacts by layer | aMUA responses',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_layers-foldchange',BRdatafile), '-jpg', '-transparent');

%% Semilog 

figure;
format bank;
subplot(2,3,1)
semilogx(layers_MON.full(:,1));
hold on
semilogx(layers_BIN.full(:,1));
xticklabels('');
ylim([-10 100]);
title('Supragranular');
hold off

subplot(2,3,2)
semilogx(layers_MON.full(:,2));
hold on
semilogx(layers_BIN.full(:,2));
xticklabels('');
ylim([-10 100]);
title('Granular');
hold off

subplot(2,3,3)
semilogx(layers_MON.full(:,3));
hold on
semilogx(layers_BIN.full(:,3));
ylim([-10 100]);
xticklabels('');
title('Infragranular');
hold off

subplot(2,3,4)
semilogy(contrast,layers_MON.full(:,1));
hold on
semilogy(contrast,layers_BIN.full(:,1));
ylim([0 100]);
xlabel('contrast');
xticklabels(0:0.5:1);
hold off

subplot(2,3,5)
semilogy(contrast,abs(layers_MON.full(:,2)));
hold on
semilogy(contrast,abs(layers_BIN.full(:,2)));
ylim([0 100]);
xticklabels(0:0.5:1);
xlabel('contrast');
hold off

subplot(2,3,6)
semilogy(contrast,layers_MON.full(:,3));
hold on
semilogy(contrast,layers_BIN.full(:,3));
ylim([0 100]);
xticklabels(0:0.5:1);
xlabel('contrast');
hold off

sgtitle({'Binned contacts by layer | aMUA responses',BRdatafile},'Interpreter','none');

%legend('Binocular stimulation','Monocoular stimulation','Location','southoutside');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_semilog',BRdatafile), '-jpg', '-transparent');

%% Binned percent difference and change
trans_v_sust = [((layers_BIN.transient(6,:)-layers_MON.transient(6,:))./(layers_MON.transient(6,:)));((layers_BIN.sustained(6,:)-layers_MON.sustained(6,:))./(layers_MON.sustained(6,:)))]';
figure;
b = bar([1 2 3],trans_v_sust,'FaceColor','Flat');
b(1).FaceColor = [.85 .85 .85]; b(1).EdgeColor = 'k';
b(2).FaceColor = [.2 .2 .2]; b(2).EdgeColor = 'k';

hold on
ylabel('fold change');
set(gca,'box','off');
xticklabels({'Supragranular','Granular','Infragranular'});
title({'From Monocular to Binocular:'...
    'Relative fold change in MUA at full contrast'});
legend('Transient','Sustained','Location','northeast');
hold off

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_transVsust',BRdatafile), '-jpg', '-transparent');
