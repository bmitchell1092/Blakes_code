function [norm] = calcNormData(data,set)
% data, samples x CHANS (within) OR samples X conditions (up to 2)

% normalize data so that values fall between 0 and 1
% use "within" case to normalize data within a channel
% use "across" case to nroamlize data across two conditions (i.e., 1 will
% correspond to the max value among the two conditions) 

switch set
    
    case 'within'
        
        for chan = 1:size(data,2)
            
            maxval  = max(data(:,chan));
            minval  = min(data(:,chan));
            
            norm(:,chan) = (data(:,chan) - repmat(minval,size(data,1),1))./...
                           repmat(maxval-minval,size(data,1),1); 
        end
        
    case 'across'
        
        
            maxval  = max(max(data));
            minval  = min(min(data));
        
        for chan = 1:size(data,2)
            
            norm(:,chan) = (data(:,chan) - repmat(minval,size(data,1),1))./...
                           repmat(maxval-minval,size(data,1),1); 
        end
        
        
end