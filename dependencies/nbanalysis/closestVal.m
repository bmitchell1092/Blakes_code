function [vecval,idx] = closestVal(vec,val)

%find index in vec that corresponds to the value in the array closest to val


for i = 1:length(val)
    diffvec = abs((vec - val(i)));
    if all(isnan(diffvec))
        vecval(i) = nan;
        idx(i) = nan;
    else
        [~,h_idx] = min(diffvec);
        vecval(i)  = vec(h_idx);
        idx(i) = h_idx;
    end
    
end

