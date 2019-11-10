function [rsVar,lims] = rsEqChans(var,int); 

totchan = size(var,2); 

ymx = int * (totchan +1); 

inter = linspace(-int,-ymx,totchan+1); 

% for i = 2:length(inter)
%     
%     meds(i) = median(inter(i):inter(i-1));
% end
% 
meds = inter; 
meds(1) = []; 

for chan = 1:totchan 
    
    mn = var(1,chan); 
    %mx = max(var(:,chan)); 
    
    %rng = mx - mn; 
    
    const = meds(chan) - mn; 
    
    rsVar(:,chan) = var(:,chan) + const; 
    
 clear rng mx mn const   
    
end

lims = [(inter(end) - 1.5*int) (inter(1) + 1.5*int)]; 
