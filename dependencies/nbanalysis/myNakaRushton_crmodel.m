function [error] = myNakaRushton_crmodel(params)

% http://jn.physiology.org/content/114/4/2087.long
% Gr = 100; %multiplicative response gain factor (=highest response amplitude)
% Gc = 50;  %normalization pool (determines x position)
% q = 10;   %exponent that determines rise and saturation
% s = 0;    %exponent that determines rise and saturation
% b = 0;    %baseline offset w/o stim

global Data 

Gr = params(1); 
Gc = params(2); 
q  = params(3); 
b  = params(4); 
s = 0; % not varying s here...

contrasts = Data(1,:); 
for ci = 1:length(contrasts)
 
    c = contrasts(ci); 
    prd(ci) = Gr*[c^(q+s)]/[c^q + Gc^q]+b; % prediction 
    
end


sse = sum((Data(2,:) - prd).^2);
error = sse; 

% don't let negative values through (set the error high):
if (sum(params<0) >= 1) || ~isreal(sse)
    error = 9999999; 
end

