% July 2014
% simple filtering of BR data (NS6) into a LFP and MUA signal
% decimated to 1kHz
% steps

function preprocessBRconv(fileinput,savepath)

% get NPMK on path
currentFolder = pwd;
cd('/users/kaciedougherty/documents/code/blackrock') %!! Change to NMPK directory 
setupNPMK('/users/kaciedougherty/documents/code/blackrock')
cd(currentFolder)

% get filename
if nargin == 1 | nargin == 2
    [datapath,datafile,~] = fileparts(fileinput); % will ignore extension and query below
else
    cd('/Volumes/DROBO/DATA/NEUROPHYS')
    [fileinput, pathinput] = uigetfile('*.*');
    [datapath,datafile,~] = fileparts([pathinput '/' fileinput]);
    cd(currentFolder);
end

filename = fullfile(datapath,datafile);
if ~exist(savepath)
savepath = '/Volumes/DROBO/USERS/Kacie/Analysis/preprocessed/';
end


fprintf('\npreprocessing %s\n',filename)

% check if file exist and load NEV w/o spiking data
if exist(strcat(filename,'.nev'),'file') == 2;
    NEV = openNEV(strcat(filename,'.nev'),'noread','nomat','nosave');
else
    error('the following file does not exist\n%s.nev',filename);
end

% check for continious file NS#
queryuser = 0;
if exist(strcat(filename,'.ns6'),'file') == 2;
    extension = 'ns6';
elseif exist(strcat(filename,'.ns5'),'file') == 2;
    % read in header info and see if data was filtered at collection
    NS = openNSx(strcat(filename,'.ns5'),'noread','precision','double');
    if mode([NS.ElectrodesInfo.HighFreqCorner]) == 300 && mode([NS.ElectrodesInfo.HighFreqCorner]) == 500000
        % N5 is unfiltered
        extension = 'ns5';
    else
        % NS is filtered
        queryuser = 1;
    end
    clear NS
else
    queryuser = 1;
end

% if the continous data was not NS6 or unfiltered NS5, ask user if they want to proceede
if queryuser
    prompt = 'We typically run this code on an unfiltered 30kHz signal (.ns6 or .ns5 file).\nProvide file extension (ns#, no ".") to continue or any other value to abort:\n';
    answer = input(prompt,'s');
    answer = lower(answer);
    if  strfind(answer,'ns')
        extension = answer;
    else
        return
    end
end

% get event codes from NEV, then clear
EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventTimes = NEV.Data.SerialDigitalIO.TimeStampSec * 1000; %ms, to match 1kHz
clear NEV

% Read in NS Header
NS_Header = openNSx(strcat(filename,'.',extension),'noread');
% get basic info about recorded data
neural = ~strcmp('E',{NS_Header.ElectrodesInfo.ConnectorBank}); % bank E is the BNCs on the front of the NSP
N.electrodes = length(neural);
N.neural = sum( neural);
N.analog = sum(~neural);

% % only run if data is certain length AND there's TWO electrodes:
% if  N.electrodes > 40 &  (NS_Header.MetaTags.DataDurationSec/60) > 5 % only run if data is certain length & there's TWO electrodes:
%get labels 
NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label};
NeuralInfo = NS_Header.ElectrodesInfo(neural);
BNCLabels = {NS_Header.ElectrodesInfo(~neural).Label};
BNCInfo = NS_Header.ElectrodesInfo(~neural);
%check that D2A range is same for all files (to convert appropriately)
rnge = ((length(NS_Header.ElectrodesInfo(1).MinAnalogValue:NS_Header.ElectrodesInfo(1).MaxAnalogValue))./...
    (length(NS_Header.ElectrodesInfo(1).MinDigiValue:NS_Header.ElectrodesInfo(1).MaxDigiValue)));
% if abs(rnge - 0.25) > .001
%     error('check D2A range\n'); 
% end

% get sampeling frequnecy
Fs = NS_Header.MetaTags.SamplingFreq;
nyq = Fs/2;
r = Fs/1000; % 1000 is the sampling frequency we want after decimation
 
% counters
clear act nct
act = 0;
nct = 0;

% process data electrode by electrode
for e = 1:N.electrodes
    fprintf('\nreading electrode %u of %u\n',e,N.electrodes); 
    
    clear NS DAT
    
    electrode = sprintf('c:%u',e);
    NS = openNSx(strcat(filename,'.',extension),electrode,'read');
    DAT = double(NS.Data);
    NS.Data = [];
    
    
    if e == 1
        %preallocation
        N.samples = length(DAT); %samples in header diffrent from actual data length???
        clear BNC LFP MUA
        BNC = zeros(ceil(N.samples/r),N.analog);
        LFP = zeros(ceil(N.samples/r),N.neural);
        MUA = zeros(ceil(N.samples/r),N.neural);
    end
    
    if strcmp('E',NS.ElectrodesInfo.ConnectorBank)
        datatype = 'BNC';
    else
        datatype = 'neural';
  
    end
    
    switch datatype
        
        case 'BNC'
            % analog input to NSP breakout board (BNC connectros on front)
            act = act+1;
            
            if Fs > 1000
                fprintf('\n decimating BNC, ch %u of %u \n',act, N.analog)
                BNC(:,act) = decimate(DAT,r); % downsample to 1kHz, lowpass filtered previously to avoid aliasing
            else
                BNC(:,act) = DAT;
            end
            clear DAT
            
            % save and clear when all channels have been read in
            if act == N.analog
                fprintf('\nsaving BNC...\n')
                savefile = sprintf('%s/%s.bnc',savepath,datafile);
                save(savefile,'BNC','BNCLabels','BNCInfo','EventCodes','EventTimes','-mat')
                clear BNC 
            end
                
            
        case 'neural'
            
            % data was collected on preamp
            % filter and downsample neural data
            nct = nct+1;
            
            DAT = DAT./4; 
            
            % LFP
            lpc = 200; %low pass cutoff
            lWn = lpc/nyq;
            [bwb,bwa] = butter(4,lWn,'low');
            fprintf('\nfiltering LFP, ch %u of %u \n',nct, N.neural)
            fLFP = filtfilt(bwb,bwa,DAT);  %low pass filter
            
            if Fs > 1000
                fprintf('\n decimating LFP, ch %u of %u \n',nct, N.neural)
                LFP(:,nct) = decimate(fLFP,r); % downsample to 1kHz, lowpass filtered previously to avoid aliasing
            else
                LFP(:,nct) = fLFP;
            end
            clear fLFP 
            
            % save and clear when all channels have been read in
            if nct == N.neural
                fprintf('\nsaving LFP...\n')
                savefile = sprintf('%s/%s.lfp',savepath,datafile);
                save(savefile,'LFP','NeuralLabels','NeuralInfo','EventCodes','EventTimes','-mat')
                clear LFP
            end
            
                     
            % MUA, analog
            hpc = 750;  %high pass cutoff
            hWn = hpc/nyq;
            [bwb,bwa] = butter(4,hWn,'high');
            fprintf('\nhp filtering MUA, ch %u of %u \n',nct, N.neural)
            hpMUA = abs(filtfilt(bwb,bwa,DAT)); %high pass filter &rectify
            
%             lpc = 50; %low pass cutoff
%             lWn = lpc/nyq;
%             [bwb,bwa] = butter(4,lWn,'low');
%             fprintf('\nlp filtering MUA, ch %u of %u \n',nct, N.neural)
%             lpMUA = filtfilt(bwb,bwa,lpMUA);  %low pass filter to smooth
            
            if Fs > 1000
                fprintf('\ndownsampling MUA, ch %u of %u \n',nct, N.neural)
                tbcMUA = downsample(hpMUA,r); % downsample to 1kHz, lowpass filtered previously to avoid aliasing
                k = jnm_kernel('psp',10); 
                
                MUA(:,nct) = doConv(tbcMUA,k)*10^3; 
            else
                tbcMUA = hpMUA; 
                
                MUA(:,nct) = doConv(tbcMUA,k)*10^3;
            end
            clear lpMUA hpMUA
             

            
            %Followng Self et al. 2013 supplement
            % MUA, analog
            % high pass at 500 Hz
            hpc = 500;  %high pass cutoff
            hWn = hpc/nyq;
            [bwb,bwa] = butter(4,hWn,'high');
            fprintf('\nhp filtering MUA, ch %u of %u \n',nct, N.neural)
            hpMUA = filtfilt(bwb,bwa,DAT); %high pass filter 
            
            % low pass at 5000 Hz and rectifiy
            lpc = 5000;  %low pass cutoff
            lWn = lpc/nyq;
            [bwb,bwa] = butter(4,lWn,'low');
            hpMUA = abs(filtfilt(bwb,bwa,hpMUA)); %low pass filter &rectify
            
            
            % low pass filter at 200 Hz (as in Self et al. 2013)
            lpc = 50; %low pass cutoff
            lWn = lpc/nyq;
            [bwb,bwa] = butter(4,lWn,'low');
            fprintf('\nlp filtering MUA, ch %u of %u \n',nct, N.neural)
            lpMUA = filtfilt(bwb,bwa,hpMUA);  %low pass filter to smooth
            
            % then downsample:
            if Fs > 1000
                fprintf('\ndownsampling MUA, ch %u of %u \n',nct, N.neural)
                selfMUA(:,nct) = downsample(lpMUA,r); % downsample to 1kHz, lowpass filtered previously to avoid aliasing
            else
                selfMUA(:,nct) = lpMUA;
            end
            clear lpMUA hpMUA
            
        
            
            % save and clear when all channels have been read in
            if nct == N.neural
                fprintf('\nsaving MUA...\n')
                savefile = sprintf('%s/%s.mua',savepath,datafile);
                save(savefile,'MUA','selfMUA','NeuralLabels','NeuralInfo','EventCodes','EventTimes','-mat')
                clear MUA
            end
           
           
    end
    
end
    
%end










