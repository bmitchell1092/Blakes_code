function [sem] = calcSEM(data)


numdim = ndims(data); 

sem    = nanstd(data,0,numdim)./sqrt(size(data,numdim));  

end