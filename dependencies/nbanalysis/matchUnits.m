% kilo sort date
clear
date = '180423';
monkey = 'I'; 
cd('/users/kaciedougherty/documents/neurophysdata/kiloout')

kfles    = dir(strcat(date,'*','ss.mat')); 
empArray = nan(length(kfles)*5,61);
clustIDs = nan(40,1); 
fileIDs = nan(40,1); 
count    = 0;

% dump waveforms across files into a matrix
for f = 1:length(kfles)
    clear ss
    load(sprintf('/users/kaciedougherty/documents/neurophysdata/kiloout/%s',kfles(f).name));
    if f == 1
        st = 1;
    else
        st = count + 1;
    end
    count                = count + size(ss.spikeWaves(1,:,ss.clusterMap(:,3) == 1),3);
    empArray(st:count,:) = squeeze(ss.spikeWaves(1,:,ss.clusterMap(:,3) == 1))'; 
    clustIDs(st:count)   = ss.clusterMap(ss.clusterMap(:,3) == 1,1); 
    fileIDs(st:count)    = f; 
end
empArray(count+1:end,:) = [];
clustIDs(count+1:end)   = []; 
fileIDs(count+1:end)    = []; 

% by eye groupings
one = [2 3 7 8 9 10 12 13 15];
two = [4 6 11 14 ];

% time to peak:
rmclust = [];
for i = 1:size(empArray,1)
    [fp] = findpeaks(empArray(i,:));
    if length(fp) > 5
        rmclust = [rmclust i]; 
    end
end
empArray(rmclust,:) = []; 
clustIDs(rmclust)   = []; 
fileIDs(rmclust)    = []; 

clear pks
pks   = nan(3,size(empArray,1));
first = [1:30];
second = [31:60];
for i = 1:size(empArray,1)
    clear fmx mn smxid
    [~,fmx]   = max(empArray(i,first));
    [~,mn]    = min(empArray(i,:));
    [~,smxid] = max(empArray(i,second));
    
    pks(1,i) = fmx; % first peak
    pks(2,i) = mn;  % min
    pks(3,i) = second(smxid);  % second peak
    pks(4,i) = mn -fmx; 
    pks(5,i) = abs((median([fmx:mn ]) ) - (median([mn:second(smxid)]) )); 
    
end

nclust = 2; 
idx = kmeans(pks',nclust); 

figure, 
for i = 1:nclust
subplot(1,nclust,i)
plot([1:61],empArray(idx == i,:)); 
end

 for i = 1:nclust
     
     id = find(idx == i); 
     groups(i).files = {kfles(fileIDs(id)).name}; 
     groups(i).clust = clustIDs(id); 
     groups(i).wave  = empArray(id,:); 
     
 end
 
 save(strcat('/users/kaciedougherty/documents/neurophysdata/',date,'_',monkey,'/groups.mat'),'groups'); 
 
 
 
 




