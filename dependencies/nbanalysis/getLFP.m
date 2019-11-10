function [xLFP,varargout] = getLFP(filename,extension,el,sortdirection)

if nargin < 4
    sortdirection = 'ascending';
end

% if ispc
%      drname = strcat('\\129.59.230.179\CerebrusData\',BRdatafile(1:8)); %sprintf('Z:\\%s',BRdatafile(1:8));
%  else
% %     if ~isempty(strfind(BRdatafile,'_I_'))
% %         rignum = '021';
% %     else
% %         rignum = '022';
% %     end
% %     d = datenum(BRdatafile(1:6),'yymmdd');
% %     if d < datenum('160831','yymmdd')
% %         drobo='Drobo';
% %     else
% %         drobo='Drobo2';
% %     end
% %     
%     drname = sprintf('/users/kaciedougherty/documents/neurophysdata/%s',BRdatafile(1:8));
% end
% THIS CODE DOES NOT DOWNSAMPLE OR FILTER DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%analyze analog data in response to stimulus presentations:
%load digital codes and neural data:
%check if file exist and load NEV w/o spiking data
if nargout > 1
    if exist(strcat(filename,'.nev'),'file') == 2
        NEV = openNEV(strcat(filename,'.nev'),'nomat','nosave');
    else
        error('the following file does not exist\n%s.nev',filename);
    end
    % get event codes from NEV, then clear
    EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
    EventTimes = floor(NEV.Data.SerialDigitalIO.TimeStampSec .* 1000); %ms, to match 1kHz
    clear NEV
    
    varargout{1} = EventCodes;
    varargout{2} = EventTimes;
    
end
%%
% Read in NS Header
NS_Header = openNSx(strcat(filename,'.',extension),'noread');
% get basic info about recorded data
neural = strcmp({NS_Header.ElectrodesInfo.ConnectorBank},el(2));
N.neural = sum( neural);

%get labels
NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label};

if isempty(NeuralLabels)
    xLFP = [];
    return
end

% get sampeling frequnecy
Fs = NS_Header.MetaTags.SamplingFreq;

% counters
clear nct
nct = 0;

% process data electrode by electrode
for e = 1:length(neural)
    
    if neural(e) == 1
        % good electrode
        nct = nct+1;
        
        clear NS DAT
        electrode = sprintf('c:%u',e);
        NS = openNSx(strcat(filename,'.',extension),electrode,'read','uV');
        if iscell(NS.Data)
            DAT = cell2mat(NS.Data);
        else
            DAT = NS.Data;
        end
        NS.Data = [];
        
        if nct == 1
            %preallocation
            N.samples = length(DAT); %samples in header diffrent from actual data length???
            clear LFP
            LFP = zeros(ceil(N.samples),N.neural);
        end
        
        LFP(:,nct) = DAT;
        clear DAT
        
    end
    
end

% sort electrode contacts in ascending order:
for ch = 1:length(NeuralLabels)
    
    chname = strcat(sprintf('%s',el),sprintf('%02d',ch));
    id = find(~cellfun('isempty',strfind(NeuralLabels,chname)));
    if ~isempty(id)
        ids(ch) = id;
    end
end

switch sortdirection
    case 'ascending'
        xLFP = LFP(:,ids);
    case 'descending'
        xLFP = LFP(:,fliplr(ids));
    otherwise
        error('need sort direction')
end



end



