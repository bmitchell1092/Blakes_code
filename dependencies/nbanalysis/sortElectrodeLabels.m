
function [ID,ids] = sortElectrodeLabels(NeuralLabels)
% sort electrode contacts in ascending order:
            el = {'eA';'eB';'eC';'eD'};
            for e = 1:size(el,1)
                for ch = 1:length(NeuralLabels)
                    ch1 = strcat(sprintf('%s',el{e}),sprintf('%02d',ch));
                    id = find(~cellfun('isempty',strfind(NeuralLabels,ch1)));
                    if ~isempty(id)
                        ids(e,ch) = id; 
                    end
                end
            end
            
            ID = find(ids(:,1)>0);
            
            % ID 1,2,3,4 goes with A,B,C,D
            
            
end