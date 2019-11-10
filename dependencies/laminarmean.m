function [lDAT, corticaldepth, N] = laminarmean(sDAT,depth)

% sDAT should be sessions/files x samples x channels with the most superficial channel as (:,:,1)
% depth should be sessions x channels, with depth(i,ch) = 0 indicating the
% channel at the bottom of the initial sink, positive values indicating
% super, and negative values indicating deep 
% can use NaN in either input to account for uneven cahnnel numbers 


d1 = size(sDAT,1);
d2 = size(sDAT,2);
corticaldepth = max(max(depth)):diff(depth(1,1:2)):min(min(depth));
d3 = length(corticaldepth); 
lDAT = NaN(d1,d2,d3); % sessions x samples x channes, will average across 1st dem before output

for i = 1:size(sDAT,1)
    try
    [~,I,II] = intersect(corticaldepth, depth(i,:),'stable');
    lDAT(i,:,I) = sDAT(i,:,II);
    catch err
        err
    end
end

N = sum(~isnan(squeeze(lDAT(:,1,:))));
lDAT = squeeze(nanmean(lDAT,1));

