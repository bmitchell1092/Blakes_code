function [BL] = BMbasecorrect(data)
%Basecorrect - baseline correct the indexed data
%   Detailed explanation goes here
basetp = data(0:100,:); %retrieves the reference range to which the data will be baseline corrected
mean_basetp = mean(basetp, 1); 
BL = data(:,:)- mean_basetp; 
end

