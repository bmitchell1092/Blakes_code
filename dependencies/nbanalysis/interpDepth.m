function [zi, yi] = interpDepth(DAT,depth,npts)
 
% MAC, March 2017
 
% for computation
    % time must be along the 2nd demension
    % channels must be along the 1st demension
 
if size(DAT,1) > size(DAT,2)
    DAT = DAT'; 
end
 
if nargin < 3
    npts = 10;
end
 
 
x = 1:size(DAT,2);
[y, i] = sort(depth);
DAT = DAT(i,:); 
 
xi = x;
yi = y(1):abs(mode(diff(y)))/npts:y(end);
 
[X,Y] = meshgrid(x,y);
[XI,YI] = meshgrid(xi,yi);
 
z  = DAT;
zi = interp2(X,Y,z,XI,YI);