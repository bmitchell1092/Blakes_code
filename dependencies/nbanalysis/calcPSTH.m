function [edges,r,h] = calcPSTH(binsize,fs,spkvec,spktr,pre,post) 

% spike
%binsize = 4;    % in ms 
%fs      = 1000; % in Hz

ntrs    = length(unique(spktr)); 
trleng  = length(pre:post); % in ms

%Compute PSTH        
lastBin = binsize * ceil((trleng-1)*(1000/(fs*binsize)));
edges = 0 : binsize : lastBin;	
edges    = edges + pre; % edges of bins (relative to stim onset)
x = (mod(spkvec-1,trleng)+1)*(1000/fs);
r = (histc(x,edges)*1000) / (ntrs*binsize);

%Plot histogram 
figure,h = gca; h_color = getColor(5); 
axes(h);
ph=bar(edges(1:end-1),r(1:end-1),'histc');
set(ph,'edgecolor',h_color,'facecolor',h_color);