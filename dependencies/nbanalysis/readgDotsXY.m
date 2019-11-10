function dot = readgDotsXY(filename)
 
%% output is structure array describing stimuli on each trial

% MAC, DEC 2014
% Made with help from MATLAB import data
% modified by Kacie for analyzing dot mapping code 

%% Check Input
[~,BRdatafile,ext] = fileparts(filename); 
if ~any(strcmp(ext,{'.gDotsXY_di'}));
    error('wrong filetype for this function')
end

%% Initialize variables.
delimiter = '\t';
endRow = inf;

n = datenum(BRdatafile(1:6),'yymmdd');

% CORRECTIONS for File-Specific Issues
if n == datenum('01/28/2016','mm/dd/yyyy');
    startRow = 2;
elseif n < datenum('01/29/2016','mm/dd/yyyy');
    startRow = 12; % for files that have weird starts 
else
    startRow = 2;
end

%% Format string for each line of text:
% For more information, see the TEXTSCAN documentation.

fields = {...
    'trial'...
    'horzdva'...
    'vertdva'...
    'dot_x'...
    'dot_y'...
    'dot_eye'...
    'diameter'...
    };

if n < datenum('01/29/2016','mm/dd/yyyy');
    formatSpec = '%u\t%f\t%f\t%f\t%f\t%f\t%f\r\n';
elseif n < datenum('170310','yymmdd'); 
    formatSpec = '%u\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n';
    fields = horzcat(fields,'contrast','timestamp');
else
     formatSpec = '%u\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n';
    fields = horzcat(fields,'contrast','fix_x','fix_y','timestamp');
end


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

% specific file corrections
switch BRdatafile
    case '160226_I_dotmapping008'
        st = 1;
        en = length(dataArray{1}) -1;
    case {'170324_I_dotmapping002','170317_I_dotmapping004','170314_I_dotmapping002','170309_I_dotmapping009','170315_I_dotmapping004','170309_I_dotmapping005','170309_I_dotmapping007','170309_I_dotmapping008','170314_I_dotmapping001','170315_I_dotmapping001','170411_I_dotmapping002','170411_I_dotmapping003','170411_I_dotmapping004'}
        st = 5;
        en = length(dataArray{1});
    otherwise
        st = 1;
        en = length(dataArray{1});
end


for f = 1:length(fields)
    if isnumeric(dataArray{f})
        dot.(fields{f}) = double(dataArray{f}(st:en));
    else
        dot.(fields{f}) = dataArray{f}(st:en);
    end
end

dot.filename = filename;
