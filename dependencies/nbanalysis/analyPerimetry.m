%% analyze drifting grating files (tfsf, rfori, cinteroc)

clearvars -except pre; warning off, %close all 
if ~ispc && exist('/volumes/drobo')
    addpath('/volumes/drobo/users/kacie/code/processdata_code');
    addpath('/volumes/drobo/users/kacie/code/mlanalysis');
    addpath(genpath('/volumes/drobo/lab software/neurophys analysis'));
elseif ~ispc && exist('/users/kaciedougherty/documents/code');
    addpath('/users/kaciedougherty/documents/code/mlanalysis');
     addpath(genpath('/users/kaciedougherty/documents/code/BLACKROCK'));
else
    addpath('/users/MLab/documents');
    addpath('/users/MLab/documents/mlanalysisonline');
    addpath('/users/MLab/documents/utils');
end

BRdatafile = '161027_I_perimetry011';
paradigm   = BRdatafile(10:end-3);paradigm = paradigm(1:4);

if ispc
    brdrname = strcat('\\129.59.230.179\CerebrusData\',BRdatafile(1:8));
    mldrname = strcat('\\129.59.231.215\MLData\',BRdatafile(1:8));
    
elseif exist('/users/kaciedougherty/documents/behavioraldata')
    brdrname = sprintf('/users/kaciedougherty/documents/behavioraldata/%s',BRdatafile(1:8));
    mldrname = brdrname;
else
    brdrname = sprintf('/Volumes/Drobo/DATA/NEUROPHYS/rig021/%s',BRdatafile(1:8));
    mldrname = brdrname;
    spath    = '/volumes/drobo/users/kacie/analysis/bcc';
end

ext = '.gPERIMETRY_di';

badobs = [1];

flag_RewardedTrialsOnly = false;
perim = readgPerimetry([mldrname filesep BRdatafile ext]);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load digital codes and neural data:
filename = fullfile(brdrname,BRdatafile);
if strcmp(BRdatafile(1:6),'151214')
    filename(end-19:end-14) = '151314';
end
% check if file exist and load NEV
if exist(strcat(filename,'.nev'),'file') == 2;
    NEV = openNEV(strcat(filename,'.nev'),'read','overwrite');
else
    error('the following file does not exist\n%s.nev',filename);
end
% get event codes from NEV
EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventTimes = floor(NEV.Data.SerialDigitalIO.TimeStampSec .* 1000); %ms, to match 1kHz
EventSampels = NEV.Data.SerialDigitalIO.TimeStamp;
[pEvC, pEvT] = parsEventCodesML(EventCodes,EventSampels);

% % get eye signals:
% [BNC] = getBNCData({'ainp2';'ainp3'},filename);
% eyex  = BNC(:,1); eyey = BNC(:,2); clear BNC

% sort/pick trials [before iterating units]
stimfeatures = {...
    'xpos'...
    'ypos'...
    'contrast'...
    'diameter'...
    'eye'...
    };

clear(stimfeatures{:})

% event codes for different behavioral outcomes: 
cd_correct_t     = [35 8 36 23 8 24 96 18 18 18]; correct_t =[];  
cd_correct_c     = [35 8 36 23 24 96 18 18 18];   correct_c = [];  

cd_nofix         = [35 18 18 18];               
cd_brokefix      = [35 8 97 18 18 18];

cd_notarget      = [35 8 36 23 18 18 18]; notarget = []; 
cd_broketarget   = [35 8 36 23 8 97 18 18 18]; broketarget = []; 

cd_brokedelay_c  = [35 8 36 97 18 18 18];  brokedelay_c = []; 
cd_broketarget_c = [35 8 36 23 97 18 18 18]; broketarget_c = []; 

obs = 0; clear TPs STIM
for t = 1:length(pEvC)

    
    stimon  =  ismember(pEvC{t},[23 25 27]);
    stimoff =  ismember(pEvC{t},[24 26 28]);
    
    st = pEvT{t}(stimon);
    
    if isempty(st)
        st = nan;
    end
    
    stim =  find(perim.trial == t); if any(diff(stim) ~= 1); error('check text file'); end
    
    
    for p = 1:length(st)
        
        % ignore early fixation breaks and no fixation trials:
        fixid = find(pEvC{t}==35); relcodes = pEvC{t}(fixid:end)' ;
        
        if ~isequal(relcodes,cd_nofix) && ~isequal(relcodes,cd_brokefix)
            obs = obs + 1;
        realtr(obs) = t;
            tk = 0;
            if isequal(relcodes,cd_correct_t)
                correct_t = [correct_t obs]; tk = 1;
        end
        if isequal(relcodes,cd_notarget)
            notarget = [notarget obs]; tk = 1; 
        end
        if isequal(relcodes,cd_broketarget)
            broketarget = [broketarget obs]; tk = 1;
        end
        
        
        if isequal(relcodes,cd_correct_c)
            correct_c = [correct_c obs]; tk = 1; 
        end
        if isequal(relcodes,cd_brokedelay_c)
            brokedelay_c = [brokedelay_c obs]; tk = 1;
        end
        if isequal(relcodes,cd_broketarget_c)
            broketarget_c = [broketarget_c obs]; tk = 1; 
        end
        
        if tk == 0 
            sprintf('this set of codes does not match any outcome: %u\n',relcodes)
            dfdf
        end
        
                if any(pEvC{t} == 96)
                    rew_tr(obs) = 1; % vector, 0 or 1--rewarded trial?
                    rxnt(obs) = abs(diff(pEvT{t}([find(pEvC{t} == 23) find(pEvC{t} == 23)+1])))./30;
                else
                    rew_tr(obs) = 0;
                    rxnt(obs) = nan;
                end
                if  any(t == badobs) || ...
                        (flag_RewardedTrialsOnly && ~any(pEvC{t} == 96))
                    TPs(obs,:) = [0 0];
                    
                    for f = 1:length(stimfeatures)
                        STIM.(stimfeatures{f})(obs,:) = NaN;
                    end
                else
                    
                    TPs(obs,:) = [st(p)];
                    
                    for f = 1:length(stimfeatures)
                        STIM.(stimfeatures{f})(obs,:) = perim.(stimfeatures{f})(stim(p));
            end
            
        end
    end
    
    end
end

% adjust position of stimuli based on fixation cross coordinates (i.e., find RF position of stimuli):
STIM.xpos = STIM.xpos -(perim.fix_x(1)); 
STIM.ypos = STIM.ypos -(perim.fix_y(1)); 

% target or non-target (catch)?
catch_tr = (STIM.xpos == 0 & STIM.ypos==0);

figure, h = plot(1:10,1:10,1:10,2:11,1:10,3:12,1:10,4:13,1:10,5:14,1:10,6:15,1:10,7:16,1:10,8:17); colors = get(h,'Color'); close %get colors for plotting
%% analyze data

ptc = length(correct_t)./sum(~catch_tr).*100; 
pnt = length(notarget)./sum(~catch_tr).*100;
ptb = length(broketarget)./sum(~catch_tr).*100;

pcc = length(correct_c)./sum(catch_tr).*100; 
pdc = length(brokedelay_c)./sum(catch_tr).*100; 
pct = length(broketarget_c)./sum(catch_tr).*100; 

perf = [ptc; pnt; ptb; pcc; pdc; pct]'; 
outcomes = {'target correct','no target','broke target','catch correct','broke delay','broke target (catch)'}; 
figure,
for i = 1:6 
    h = bar(i,perf(i),'FaceColor',colors{i});
    hold on; 
end
legend(outcomes); set(gca,'XTickLabel',[]); ylabel('% of target or catch trials');
setFigure(gcf,gca,[1 1 850 400]);


diameters = STIM.diameter; 
xs        = STIM.xpos;
ys        = STIM.ypos;

res = 0.05; % dva per pix in matrix
dd = max(diameters);
[X,Y] = meshgrid(min(xs)-dd:res:max(xs)+dd, min(ys)-dd:res:max(ys)+dd);
Z = NaN(size(X,1),size(X,2),length(xs));

for tr = 2:length(xs)
    xx = xs(tr);
    yy = ys(tr);
    dd = diameters(tr);
    fill = sqrt(abs(X-xx).^2 + abs(Y-yy).^2) < dd/2;
    if ~any(any(fill))
        error('check matrix')
    end
    trldat = Z(:,:,tr);
    [~,ecc] = cart2pol(xs(tr),ys(tr));
    trldat(fill) = rew_tr(tr);
    Z(:,:,tr) = trldat;
end
uZ = nanmean(Z,3);

figure,h = imagesc(X(1,:),Y(:,1),uZ);set(h,'AlphaData',~isnan(uZ));
set(gca,'YDir','normal','TickDir','out','FontSize',12)
xlabel('x position (dva)');
ylabel('y position (dva)');
cb = colorbar; set(get(cb,'ylabel'),'string','fraction correct'); 
caxis([0 1])
%pre = uZ; 
post = uZ; 
% hold on; 
% [rxs,rys] = getRFcoords;
% plot(rxs,rys,'*'); 

% number of trials per location :

res = 0.05; % dva per pix in matrix
dd = max(diameters);
[X,Y] = meshgrid(min(xs)-dd:res:max(xs)+dd, min(ys)-dd:res:max(ys)+dd);
Z = NaN(size(X,1),size(X,2),length(xs));

for tr = 2:length(xs)
    xx = xs(tr);
    yy = ys(tr);
    dd = diameters(tr);
    fill = sqrt(abs(X-xx).^2 + abs(Y-yy).^2) < dd/2;
    if ~any(any(fill))
        error('check matrix')
    end
    trldat = Z(:,:,tr);
    [~,ecc] = cart2pol(xs(tr),ys(tr));
    trldat(fill) = 1;
    Z(:,:,tr) = trldat;
end
uZ = nansum(Z,3);
uZ(uZ == 0) = nan; 
figure,h = imagesc(X(1,:),Y(:,1),uZ);set(h,'AlphaData',~isnan(uZ));
set(gca,'YDir','normal','TickDir','out','FontSize',12)
xlabel('x position (dva)');
ylabel('y position (dva)');
cb = colorbar; set(get(cb,'ylabel'),'string','n trials'); 
caxis([0 30]); 

%% plot reaction time data:
rxnt = rxnt(~catch_tr); 
rxnt = rxnt(rxnt>0);
mnrx = floor(min(rxnt)); mxrx = ceil(max(rxnt)); 
binsize = 10; 

figure,
histogram(rxnt,[mnrx:binsize:mxrx]); 
xlabel('time to target (ms)'); ylabel('n trials'); 
title(strcat(BRdatafile,'reaction time distribution'),'interpreter','none'); 
setFigure(gcf,gca,[1 1 600 300])

% %% eye movements on missed target trials:
% figure,
% pre = -100; post = 600; tvec = pre:post;
% for i = 1:length(notarget)
%     refwin = floor(TPs(notarget(i))./30) +pre : floor(TPs(notarget(i))./30) + post; 
%     plot(tvec,eyey(refwin)); 
%     %hold on; 
%     %plot(eyey(refwin)); 
%     hold on; 
% end

%%
newuZ = pre-post;
figure,h = imagesc(X(1,:),Y(:,1),newuZ);set(h,'AlphaData',~isnan(newuZ));
set(gca,'YDir','normal','TickDir','out','FontSize',12)
xlabel('x position (dva)');
ylabel('y position (dva)');
cb = colorbar; set(get(cb,'ylabel'),'string','difference'); 
caxis([0 0.45]); 