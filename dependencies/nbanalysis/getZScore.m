function [zresp] = getZScore(data,baseline,maxs)


zresp = (nanmean(data,2) - nanmean(baseline))./(nanmean(maxs) - nanmean(baseline)) ;  

