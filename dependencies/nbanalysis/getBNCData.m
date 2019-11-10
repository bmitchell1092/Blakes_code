function [BNC, BNCLabels,fs] = getBNCData(theseelectrodes,filename,sampledown)

BNC       = [];
BNCLabels = [];

if nargin == 2
    sampledown = 0; % downsample to 1 kHz?
end

if ~exist(strcat(filename),'file') == 2;
    error('file does not exist\n');
end
extension = filename(end-2:end); 

% Read in NS Header
NS_Header    = openNSx(filename,'noread');
neural       = ~strcmp('E',{NS_Header.ElectrodesInfo.ConnectorBank}); % bank E is the BNCs on the front of the NSP
N.electrodes = length(neural);
N.analog     = sum(~neural);

%get labels
BNCLabels    = {NS_Header.ElectrodesInfo(~neural).Label};
BNCInfo      = NS_Header.ElectrodesInfo(~neural);

% get sampling frequency
Fs           = NS_Header.MetaTags.SamplingFreq;
nyq          = Fs/2;
r            = Fs/1000; % 1000 is the sampling frequency we want after decimation

% counters
clear act nct
act = 0;

% sort electrode contacts in ascending order:
for ch = 1:length(theseelectrodes)
    chname = theseelectrodes{ch};
    id = find(~cellfun('isempty',strfind(BNCLabels,chname)));
    if ~isempty(id)
        ids(ch) = id;
    end
end

for e = ids
    
    electrode = sprintf('c:%u',( e + sum( neural ) ) ); % Jake Edit
    NS        = openNSx(strcat(filename),electrode,'read','uV');
    DAT       = double(NS.Data);
    NS.Data   = [];
    
    % analog input to NSP breakout board (BNC connectros on front)
    act = act+1;
    
    if Fs > 1000 && sampledown == 1
        BNC(:,act) = decimate(DAT,r); % downsample to 1kHz, lowpass filtered previously to avoid aliasing
        fs = 1000; 
    else
        BNC(:,act) = DAT; 
        fs = Fs;
    end
    clear DAT
    
end