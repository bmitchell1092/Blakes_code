function perim = readgPerimetry(filename, startRow, endRow)
 
%% output is structure array describing stimuli on each trial

% MAC, DEC 2014
% Made with help from MATLAB import data
% modified by Kacie for analyzing dot mapping code 
% additional revisions Jan 2016 to add more features / corrections, make backwards compatable

%% Check Input
[~,BRdatafile,ext] = fileparts(filename); 
if ~any(strcmp(ext,{'.gPERIMETRY_di'}));
    error('wrong filetype for this function')
end


%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format string for each line of text:
% For more information, see the TEXTSCAN documentation.

    fields = {...
        'trial'...
        'horzdva'...
        'vertdva'...
        'xpos'...
        'ypos'...
        'contrast'...
        'diameter'...
        'eye'...
        'header'...
        'timestamp'...
        'fix_x',...
        'fix_y'};
    
    formatSpec = '%04u\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\r\n';
    

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this code.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Allocate imported array to structure column variable names
if length(fields) ~= size(dataArray,2)
    error('bad formatSpec or structure fields for %s',filename)
end

    st = 1;
        en = length(dataArray{1});

% CORRECTIONS for File-Specific Issues

for f = 1:length(fields)
    if isnumeric(dataArray{:, f})
        perim.(fields{f}) = double(dataArray{f}(st:en));
    else
        perim.(fields{f}) = dataArray{f}(st:en);
    end
end

ntrls          = max(perim.trial); % total trials
npres          = mode(histc(perim.trial,1:max(perim.trial))); % number of "gen" calls written / trial, may be diffrent from RECORD & what was actually shown
perim.pres   = repmat([1:npres]',ntrls,1);

perim.filename = filename;
perim.startRow = startRow;
perim.endRow = endRow;



