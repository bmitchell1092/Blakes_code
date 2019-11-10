function [onlat,peaklat,tridx] = getLatTrans(spkTPs,contrasts,SPK,pre,post,fs)

thrc = max(contrasts);  
HspkTPs = spkTPs(contrasts >= thrc,:); 

RESP = NaN(length(HspkTPs),1);

    trspks = [];
    trls = [];
    for r = 1:length(HspkTPs)
        st = HspkTPs(r,1);
        en = HspkTPs(r,2);
        if st ~= 0 && en ~= 0   
            refwin = st + pre: st + post;
            id = find(ismember(refwin,SPK(SPK>refwin(1) & SPK<refwin(end))));
            trls = [trls repmat(r,length(id),1)'];
            trspks = [trspks id];
        end
    end
    
    
    % psth:
    binsize = 5; % ms
    triallen = length([pre:post]);
 
    ntrials = length(HspkTPs);
    
    lastBin = binsize * ceil((triallen-1)*(1000/(fs*binsize)));
    edges = 0 : binsize : lastBin;
    tedges = edges +(pre./double((fs/1000)));
    x = (mod((trspks)-1,triallen)+1)*(1000/fs);
    resp = (histc(x,edges)) / (ntrials*binsize) .* (1000/binsize);
    
%     %Plot histogram
%     h = gca; h_color = [.2 .6 .4];
%     axes(h);
%     figure 
%     ph=bar(tedges(1:end-1),resp(1:end-1),'histc');
%     set(ph,'edgecolor',h_color,'facecolor',h_color);
%     set(gca,'XLim',[tedges(1) tedges(end)],'Box','off','TickDir','out');
%     ylabel('spks/sec'); xlabel('t (ms) from stimulus onset');
    
    %%
    % neuron latency:
    % threshold: 3*STDs above base firing rate
    
    m_bsl   = nanmean(resp(tedges<=0));
    std_bsl = nanstd(resp(tedges<=0));
    thr     =(2.75*m_bsl);
    faker = resp; faker(tedges<=0) = 0;
    onlat   = tedges(find(faker>thr,1,'first')); % first bin passing thr;     % Jiang et al. (2015): P (25.2, 2.0) M  (19.6, 2.2)
    peaklat = tedges(find(faker == max(faker))); % bin with max resp    % Jiang et al. (2015): P (72.9, 5.3) M  (43.5, 3.3)
    if isempty(onlat)
        fprintf('\nNote: No data points cross onset threshold\n');
    end

    % neuron transiency:
    % Jiang et al. (2015): P (22.95, 5.97) M  (39.53,7.5)
    tridx = 100 - (nanmean((resp(tedges>=100 & tedges <=200))) - nanmean(resp(tedges<0))) ./ (nanmean((resp(tedges>=0 & tedges <=100))) - nanmean(resp(tedges<0))).*100;
    
    fprintf('\nonset latency: %d ms\npeak latency: %d ms\ntransiency index %d\n\n',onlat,peaklat,tridx);
    
    %Plot histogram
    figure; h = gca; h_color = [.2 .6 .4];
    axes(h);
    ph=bar(tedges(1:end-1),resp(1:end-1),'histc');
    if~isempty(onlat)
        h1 = vline(onlat);
        set(h1,'Color','k');
    end
    h2 = vline(peaklat);  set(h2,'Color','k');
    set(ph,'edgecolor',h_color,'facecolor',h_color);
    set(gca,'XLim',[tedges(1) tedges(end)]);
    ylabel('spks/sec'); xlabel('t (ms) from stimulus onset');
    title(gca,sprintf('%d trls/lat %0.2f ms/peak lat %0.2f ms/trans idx %0.2f',length(spkTPs),onlat,peaklat,tridx));
    