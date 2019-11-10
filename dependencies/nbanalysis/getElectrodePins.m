function [pins] = getElectrodePins(NEV,el,chans)

if nargin == 2
    chans = [1:24]; 
end

labels = {NEV.ElectrodesInfo.ElectrodeLabel};
for  ch = chans
    elabel = sprintf('%s%02u',el,chans(ch));
    eidx = find(cell2mat(cellfun(@(x) ~isempty(strfind(x',elabel)),labels,'UniformOutput',0)));
    pins(ch) = NEV.ElectrodesInfo(eidx).ElectrodeID;
end