function [Gr Gc q s b fbest] = fitCRdata

% can increase number of starting points or 'maxiter' for possibly better fits
options = optimset('display', 'off', 'MaxIter', 5000);
nstarting_pts = 500; 

fbest = 999;

% params = [100 50 10 0 0]; Gr Gc q s b from online example 

for i = 1:nstarting_pts
    
    % initial starting paramaters for Gr Gc q b: 
    
    parInit(1:2) =  randi([1 200],2,1); %rand(1).*100;
    parInit(3)   = rand(1); %rand(1).*10; %
    parInit(4)   = randi([1 10],1,1);   
   
    [xret,fret,exitflag,output] = fminsearch(@myNakaRushton_crmodel, parInit, options);
    if i == 1
        xbest = xret; 
    end
    if fret < fbest
        fbest = fret;
        xbest = xret;
    end
    
end

Gr = xbest(1);    %multiplicative response gain factor (=highest response amplitude)
Gc = xbest(2);    %normalization pool (determines x position)
q  = xbest(3);    %exponent that determines rise and saturation
s  = 0 ;          %exponent that determines rise and saturation---not varying here
b  = xbest(4);    %baseline offset w/o stim

end