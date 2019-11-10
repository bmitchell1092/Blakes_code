function brfs = readgDRFTBRFS(filename, startRow, endRow)
 
%% output is structure array describing stimuli on each trial

%% Check Input
[~,~,ext] = fileparts(filename); 
if ~any(strcmp(ext,{'.gDRFTBrfsGratings_L1','.gDRFTBrfsGratings_L2','.gDRFTBrfsGratings_R1','.gDRFTBrfsGratings_R2'}));
    error('wrong filetype for this function')
end

%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end
%us

%% Format string for each line of text:
% For more information, see the TEXTSCAN documentation.

    formatSpec = '%u%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';


%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
% for block=2:length(startRow)
%     frewind(fileID);
%     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
%     for col=1:length(dataArray)
%         dataArray{col} = [dataArray{col};dataArrayBlock{col}];
%     end
% end

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names

brfs.trial              = dataArray{:,1};
brfs.horzdva            = dataArray{:,2};
brfs.vertdva            = dataArray{:,3};
brfs.tilt               = dataArray{:,4};
brfs.sf                 = dataArray{:,5};
brfs.phase              = dataArray{:,6};
brfs.disparity          = dataArray{:,7};
brfs.diameter           = dataArray{:,8};
brfs.eye                = dataArray{:,9};
brfs.eye_contrast       = dataArray{:,10};
brfs.other_contrast     = dataArray{:,11};
brfs.oridist            = dataArray{:,12};
brfs.gaborfilter_on     = dataArray{:,13};
brfs.gabor_std          = dataArray{:,14};
brfs.total_stim_dur     = dataArray{:,15};
brfs.soa_dur            = dataArray{:,16};
brfs.isi_dur            = dataArray{:,17};
brfs.temporal_freq      = dataArray{:,18};
brfs.timestamp          = dataArray{:,19};

brfs.filename = filename;
brfs.startRow = startRow;
brfs.endRow = endRow;
