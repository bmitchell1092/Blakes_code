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
    datadir  = '/users/bmitc/Box Sync/DATA/';
end

addpath(genpath(npmkdir))
addpath(genpath(nbanalysisdir))
addpath(genpath(datadir))

BRdatafile = '160202_I_cinteroc004'; 
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
elseif isequal(BRdatafile,'161003_E_cinteroc002')
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
else 
    ebank = 'eB';
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

%% Averaging across trials

avg.LFP = mean(STIM.LFP,3);
avg.aMUA = mean(STIM.aMUA,3);
avg.CSD = mean(STIM.CSD,3);

%% Baseline correct 

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

%% Z-score Normalization for aMUA

zMUA = nan(size(STIM.aMUA,1),size(STIM.aMUA,2),size(STIM.aMUA,3));

for t = 1:size(STIM.aMUA,3)
    for c = 1:size(STIM.aMUA,2)
        zMUA(:,c,t) = (STIM.aMUA(:,c,t)-mean(STIM.aMUA(25:75,c,t)))./(std(STIM.aMUA(25:75,c,t)));
    end
end

avg.zMUA = mean(zMUA,3);
bsl.zMUA = BMbasecorrect(avg.zMUA);


%% Averaged trials by condition
%  Now that I have my conditions as logicals in a structure and my aMUA in
%  z-scores, time to pull out AVG'd z-scored aMUA for each condition

clear m MzMUA BzMUA
for m = 1:size(STIM.Mconditions,1)
    MzMUA(m).contrast = mean(zMUA(:,:,STIM.Mconditions(m,:)),3); %#ok<SAGROW>
    BzMUA(m).contrast = mean(zMUA(:,:,STIM.Bconditions(m,:)),3); %#ok<SAGROW>
end

%% collapsing across time for each condition
clear i coll_mon
for i=1:size(MzMUA,2)
    coll_mon.contrast(i,:)  = mean(MzMUA(i).contrast(80:330,:),1);
end

clear i coll_bin
for i=1:size(BzMUA,2)
    coll_bin.contrast(i,:)  = mean(BzMUA(i).contrast(80:330,:),1);
end

%% Plotting

%% Plotting all averaged, baseline corrected trials (SNAPSHOT)

refwin = pre:post; % reference window for line plotting
channels = 1:nct;  % how many channels (nct is a predefined variable with the exact number of channels
offset = ((STIM.offsets(1)-STIM.onsets(1))/(30));

h1 = figure('position',[15,135,1200,500]);
clear i
avg_fields = fieldnames(avg);
for i = 1:length(avg_fields)
subplot(1,4,i)
f_ShadedLinePlotbyDepth((avg.(avg_fields{i})),1:24,refwin,channels,1);
hold on
plot([0 0], ylim,'k')
plot([offset offset], ylim,'k','linestyle','-.','linewidth',0.5)
title(avg_fields{i})
xlabel('time (ms)')
ylabel('contacts indexed down from surface')
hold off
end

bAVG_iCSD = filterCSD(bsl.CSD')';

subplot(1,4,4)
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

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_snapshot',BRdatafile), '-jpg', '-transparent');

%% Varying contrast, monocular stimulation (Contrasts)
h2 = figure('position',[15,135,1200,500]);
clear i 
for i = 1:length(contrast)
subplot(1,length(contrast),i);
f_ShadedLinePlotbyDepth(mean(STIM.aMUA(:,:,STIM.Mconditions(i,:)),3),1:24,refwin,channels,1)
hold on
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
format bank; rcontrast = round(contrast,2,'significant');

h4 = figure('Position', [60 211 1100 300]);
subplot(1,4,1)
bar(coll_bin.contrast(:,5),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(coll_mon.contrast(:,5),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-1 5.5]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('Z-normalized MUA')
title('contact 5');
hold off

subplot(1,4,2)
bar(coll_bin.contrast(:,10),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(coll_mon.contrast(:,10),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-1 5.5]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('Z-normalized MUA')
title('contact 10');
hold off

subplot(1,4,3)
bar(coll_bin.contrast(:,15),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(coll_mon.contrast(:,15),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-1 5.5]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('Z-normalized MUA')
title('contact 15');
legend('Binocular stimulation','Monocoular stimulation','Location','southoutside');
hold off

subplot(1,4,4)
bar(coll_bin.contrast(:,20),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
bar(coll_mon.contrast(:,20),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
set(gca,'box','off');
ylim([-1 5.5]);
xticklabels(rcontrast)
xlabel('contrast level')
ylabel('Z-normalized MUA')
title('contact 20');
hold off

sgtitle({'Monocular vs Binocular response as a function of contrast',BRdatafile},'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_bar-contrasts',BRdatafile), '-jpg', '-transparent');


%% 3D surface plot (MESH)
format bank; rcontrast = round(contrast,2,'decimals');

h7 = figure('Position', [49 111 1100 470]);
sp2 = subplot(1,3,2);
surf(contrast,channels,coll_bin.contrast');
hold on
colormap(gca, 'jet'); % this makes the red color the sinks and the blue color the sources (convention)
cm2 = colorbar; 
set(gca,'tickdir','out'); 
climit = max(abs(get(gca,'CLim'))*1);
%climit = 1;
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
title('Binocular response')
xlabel('contrast level')
clrbar = colorbar; clrbar.Label.String = 'z-score'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts');
xticklabels(0:0.5:1);
% zaxis = get(gca,'zlim');
zaxis = [-1 5.5];
set(gca,'zlim',zaxis);
hold off

sp1 = subplot(1,3,1);
surf(contrast,channels,coll_mon.contrast');
hold on
colormap(gca,'jet'); % this makes the red color the sinks and the blue color the sources (convention)
cm1 = colorbar; 
set(gca,'zlim',zaxis,'tickdir','out');
%climit = max(abs(get(gca,'CLim'))*1);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
%zlim = zaxis;
title('Monocular response')
xlabel('contrast level')
clrbar = colorbar; clrbar.Label.String = 'z-score'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts');
xticklabels(0:0.5:1);
hold off

sp3 = subplot(1,3,3);
surf(contrast,channels,coll_bin.contrast'-coll_mon.contrast');
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
clrbar = colorbar; clrbar.Label.String = 'z-score difference'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('contacts');
xticklabels(0:0.5:1);
zaxis = get(gca,'zlim');
set(gca,'zlim',zaxis);
hold off

sgtitle({'Monocular versus Binocular response as a function of contrast',BRdatafile},'Interpreter','none');

set(sp1, 'Position', [.08 .2 .20 .6])
set(sp2, 'Position', [.38 .2 .20 .6])
set(sp3, 'Position', [.68 .2 .20 .6])

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_surf',BRdatafile), '-jpg', '-transparent');

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
title('2D Binocular contrast response');
xlabel('\fontsize{12}contrast level')
clrbar = colorbar; clrbar.Label.String = 'z-score'; 
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
title('2D Monocular contrast response');
xlabel('\fontsize{12}contrast level')
clrbar = colorbar; clrbar.Label.String = 'z-score'; 
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
clrbar = colorbar; clrbar.Label.String = 'z-score difference'; 
set(clrbar.Label,'rotation',270,'fontsize',10,'VerticalAlignment','middle');
ylabel('\fontsize{12}contacts indexed down from surface');
xticklabels(rcontrast);
hold off

sgtitle({'Monocular vs Binocular response as a function of contrast',BRdatafile},'Interpreter','none');

set(sp1, 'Position', [.10 .2 .20 .6])
set(sp2, 'Position', [.40 .2 .20 .6])
set(sp3, 'Position', [.70 .2 .20 .6])

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_imagesc',BRdatafile), '-jpg', '-transparent');

%% Experimental tightplot
h9 = figure('Position', [280,58,380,582]);

[ha, pos] = tight_subplot(24,1,[0.005 .03],[.1 .15],[.2 .2]); %channels, columns, [spacing], [top and bottom margin], [left and right margin]
for c = 1:24
    
    axes(ha(c)); % ha is a variable that gets the axis of each subplot
    
    plot(coll_mon.contrast(:,c),'k','linewidth',.5);
    hold on
    plot(coll_bin.contrast(:,c),'.-b','linewidth',0.5);
    %fill([contrast fliplr(contrast)], [coll_mon.contrast(:,c) fliplr(coll_bin.contrast(:,c))], 'b')
    hold off
    ylim([-0.5 5]);
    for cc =1:5:24
    yticklabels(c);
    end
    grid off
    
    if c < 24
    set(gca, 'XTick', '') %removing the units on the x axis (if i < 24) 
    % set(gca, 'YTick', '') %removing the units on the x axis (if i < 24)
    p = gca; p.XAxis.Visible = 'off'; %remove the x axis
    end
    
    if c == 1
        title({'Contrast response curves',BRdatafile},'Interpreter', 'none')
    end
   
    if c == 12
        ylabel('\fontsize{12}contacts in order of depth');
    end
    
    if c == 24
        xticklabels(contrast)
    end
end

set(ha(1:23), 'XTickLabel',''); %remove labels from first 23 plots
set(ha(1:24), 'box', 'off');

xlabel('\fontsize{12}contrast');
% legend('Monocular stimulus','Binocular stimulus','Location','eastoutside');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('%s_tightplot',BRdatafile), '-jpg', '-transparent');
