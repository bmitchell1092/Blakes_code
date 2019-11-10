function [IMS,rnsinfo] = readgRNS(infofname,rnsfname)

delimiter = '\t';
startRow = 2;
endRow = inf;
% %%
% infofname = '/Users/kaciedougherty/Documents/neurophysdata/170214_T/170214_T_rns008.gRNSinfo_di'; 
% rnsfname = '/Users/kaciedougherty/Documents/neurophysdata/170214_T/170214_T_rns008.gRNS_di'; 
%%
fields = {'trial','rsize','csize'};
formatSpec = '%u%u%u'; 
fileID = fopen(infofname,'r');
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
fclose(fileID);
st = 1;
        en = length(dataArray{1});
for f = 1:length(fields)
    if isnumeric(dataArray{:, f})
        rnsinfo.(fields{f}) = double(dataArray{f}(st:en));
    else
        rnsinfo.(fields{f}) = dataArray{f}(st:en);
    end
end

%%
fileID = fopen(rnsfname,'r');
formatSpec = '%f';
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1,'ReturnOnError', false);
fclose(fileID);
%%
trs = length(rnsinfo.csize); 
IMS = reshape(dataArray{1},[rnsinfo.rsize(1) rnsinfo.csize(1) trs]); 
