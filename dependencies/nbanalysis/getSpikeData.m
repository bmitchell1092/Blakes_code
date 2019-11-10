function [spikes] = getSpikeData(signals,theseelectrodes,NEV)

spikes = []; 

if any(strcmp(signals,'spikes')) || any(strcmp(signals,'SPIKES'))
for ch = 1:size(NEV.ElectrodesInfo,2)
    nevlabel{ch} = NEV.ElectrodesInfo(ch).ElectrodeLabel';
end

for ch = 1:length(theseelectrodes)
    
    clear chname; 
    chname = theseelectrodes{ch};
    nevid = find(~cellfun('isempty',strfind(nevlabel,chname)));
    
    if ~isempty(nevid)
    spkid = find(NEV.Data.Spikes.Electrode == nevid); 
    h_spkt  = NEV.Data.Spikes.TimeStamp(spkid); 
    spikes{ch} = unique(h_spkt./NEV.MetaTags.SampleRes.*1000); 
    end
    
end
end