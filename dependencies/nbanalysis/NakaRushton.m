% NAKA RUSHTON FUNCTION
% from: http://docs.psychtoolbox.org/ComputeNakaRushton
% http://jn.physiology.org/content/114/4/2087.long

Gr = 100; %multiplicative response gain factor (=highest response amplitude)
Gc = 50;  %normalization pool (determines x position)
q = 10;   %exponent that determines rise and saturation
s = 0;    %exponent that determines rise and saturation
b = 0;    %baseline offset w/o stim

clear Rc
clear mRc
for c=1:100 
    %compute response Rc for each conrast level c:  
    Rc(c) = Gr*[c^(q+s)]/[c^q + Gc^q]+b; 
    
    %compute multiplicative gain prediction
    Grm = 80;
    mRc(c) = Grm*[c^(q+s)]/[c^q + Gc^q]+b;
    
    %compute contrast-gain prediction
    cGc = 55;
    cRc(c) = Gr*[c^(q+s)]/[c^q + cGc^q]+b;
    
    %compute additive gain prediction
    ab = -10;
    aRc(c) = Gr*[c^(q+s)]/[c^q + Gc^q]+ab;
end

figure(1),clf
plot([1:100],Rc,'k')
hold on
plot([1:100],mRc,'r')
plot([1:100],cRc,'b')
plot([1:100],aRc,'g')
ylim([0 100])
legend('monocular CRF','multiplicative gain','contrast gain',...
    'additive gain','Location','NorthWest')


