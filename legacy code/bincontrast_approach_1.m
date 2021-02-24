%% bincontrast.m
% Loads in all ditask units. IDX is the info struct. Constructs curves for
% each unit seperately, then bins the curves by ocularity/layer. Plotting included. 

clear
% choose dataset
dataset = 'diSTIM_allLevels';

% Setup directory for files of interest
if strcmp(getenv('username'),'mitchba2')
    didir = strcat('D:\dMUA\',dataset,'\');
elseif strcmp(getenv('username'),'bmitc')
    didir = strcat('C:\Users\bmitc\Documents\MATLAB\Data\',dataset,'\');
elseif strcmp(getenv('username'),'bmitc_000')
    didir = strcat('C:\Users\bmitc_000\Documents\MATLAB\Data\',dataset,'\');
end
list    = dir([didir '*_AUTO.mat']);


% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(didir)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% script choices
flag_save = 1;
baseline_correct = 1;
balanced = 1;

% Counts
N = 0;
uct = 0;

% bad curves
badCurves.mon = 0;
badCurves.bin = 0;
badCurves.pen = [];
        
% Penetration Loop
for pen = 1:length(list)
    tic
    
    % Load penetration data
    clear penetration
    penetration = list(pen).name(1:11);
    
    load([didir penetration '.mat'],'STIM')
    matobj = matfile([didir penetration '_AUTO.mat']);
    
    win_ms = matobj.win_ms;
    if ~isequal(win_ms,[40 140; 141 450; 50 250; -50 0])
        error('check RESP window')
    end
    N = N+1; % will have a running count of penetrations
    % Electrode Loop
    for e = 1:length(STIM.depths)
        uct = uct+1;
        goodfiles = unique(STIM.filen);
        
        resp  = squeeze(matobj.RESP(e,:,:)); % pulls out matobj RESP, (e x time x trial)
        resp = squeeze(bsxfun(@minus,resp(3,:), resp(4,:)))';% baseline corrects resp(3) by resp(4)
        
        X = diUnitTuning(resp,STIM,goodfiles); %get tuning info for the unit
        
        DE = X.dipref(1); % preferred eye
        NDE = X.dinull(1); % non-preferred ete
        PS = X.dipref(2); % preferred stimulus
        NS = X.dinull(2); % null stimulus
        
        % sort data so that they are [prefeye nulleye]
        clear eyes sortidx contrasts tilts
        eyes      = STIM.eyes;
        contrasts = STIM.contrast;
        tilts     = STIM.tilt;
        if X.dipref(1) == 2
            [eyes,sortidx] = sort(eyes,2,'ascend');
        else
            [eyes,sortidx] = sort(eyes,2,'descend');
        end
        for w = 1:length(eyes)
            contrasts(w,:) = contrasts(w,sortidx(w,:)); % sort contrasts in dominant eye and non-dominant eye
            tilts(w,:)     = tilts(w,sortidx(w,:));
        end; clear w
        
        
        STIM.monocular(find(STIM.adapted)+1) = 1; % not sure if I need this anymore
        
        % establish constant conditions
        I = STIM.ditask ...
            & STIM.adapted == 0 ... %is not adapted
            & STIM.rns == 0 ... %not random noise stimulus
            & STIM.cued == 0 ... %not cued or uncued
            & STIM.motion == 0 ... %not moving
            & ismember(STIM.filen,goodfiles); % things that should be included.
        
        % pull out the data for single electrode
        clear sdf sdftm resp
        resp = squeeze(matobj.RESP(e,:,:));
        
        if baseline_correct == true
            resp  = bsxfun(@minus,resp, resp(4,:)); % resp 5 is already the mean -.50 to 0
        end
        
        % Define stimulus levels for this unit
        stimcontrast = [0 X.dicontrasts]; % 0 and all contrast levels
        numC = length(stimcontrast);
        
        
        %% Pull out data by monocular condition, one contact at a time
        
        % pre-allocate data matrices
        clear moncond monSDF monSDFerror monRESP monRESPerror monTrlNum
        moncond     = {'DE_PS','NDE_PS','DE_NS','NDE_NS'};
        monRESP    = nan(numC,4,4); % contrast x timewindow x condition
        monRESPerror  = nan(numC,4,4); % contrast x timewindow x condition
        monTrlNum     = nan(numC,4);  % contrast x condition
        
        for mon = 1:size(moncond,2) % for each condition
            for c = 1:length(stimcontrast) % for each contrast level
                switch moncond{mon}
                    case 'DE_PS'
                        if c == 1
                            trls = STIM.blank; % zero contrast in both eyes
                        else
                            trls = I & STIM.monocular & DE... % is monocular and dominant eye
                                & contrasts(:,1) == stimcontrast(c)... % contrast in dom eye
                                & tilts(:,1) == X.dipref(2); % pref orientation in dom eye
                        end
                        
                    case 'NDE_PS'
                        if c == 1
                            trls = STIM.blank;
                        else
                            trls = I & STIM.monocular & NDE...  % is monocular and non-dominant eye
                                & contrasts(:,2) == stimcontrast(c) ... % contrast in non-dom eye
                                & tilts(:,2) == X.dipref(2); % pref orientation in non-dom eye
                        end
                        
                    case 'DE_NS'
                        if c == 1
                            trls = STIM.blank;
                        else
                            trls = I & STIM.monocular & DE... % is monocular and DE
                                & contrasts(:,1) == stimcontrast(c)... % contrast in dom eye
                                & tilts(:,1) == X.dinull(2); % pref orientation in dom eye
                        end
                        
                    case 'NDE_NS'
                        if c == 1
                            trls = STIM.blank;
                        else
                            trls = I & STIM.monocular & NDE...
                                & contrasts(:,2) == stimcontrast(c)... % contrast in dom eye
                                & tilts(:,2) == X.dinull(2); % pref orientation in dom eye
                        end
                end
                
                % pass if trial numbers are greater than 4
                if sum(trls) >= 5
                    monRESP(c,:,mon) = nanmean(resp(:,trls),2);
                    monRESP(1,:,mon) = nanmean(resp(4,:),2);
                    monRESPerror(c,:,mon)   = (nanstd(resp(:,trls),0,2))./(sqrt(sum(trls)));
                end
                monTrlNum(c,mon) = sum(trls); % stores trial count by contrast and condition
                Trls.mon(c,mon,pen) = sum(trls); % stores trial count by contrast, condition, and penetration
            end
        end
        clear trls mon c cond
        
        %% Binocular conditions
        clear binSDF binSDFerror binRESP binRESPerror binTrlNum
        bincond     = {'PS','NS'}; % PS = Preferred stimulus; NS = Null stimulus
        binRESP    = nan(numC,4,2); % contrast x timewindow x condition
        binRESPerror  = nan(numC,4,2); % contrast x time x condition
        binTrlNum     = nan(4,2); % contrast x condition
        
        for bin = 1:size(bincond,2) % for each binocular condition
            for c = 1:length(stimcontrast) % for each contrast level
                switch bincond{bin}
                    case 'PS'
                        if c == 1
                            trls = STIM.blank;
                        else
                            trls = I & STIM.dioptic... % should this be STIM.dioptic?
                                & contrasts(:,1) == stimcontrast(c)... % contrast in dom eye
                                & contrasts(:,2) == stimcontrast(c)... % contrast in null eye
                                & tilts(:,1) == X.dipref(2)... % pref orientation in dom eye
                                & tilts(:,2) == X.dipref(2); % pref orientation in null eye
                        end
                    case 'NS'
                        if c == 1
                            trls = STIM.blank;
                        else
                            trls = I & STIM.dioptic...
                                & contrasts(:,1) == stimcontrast(c)... % contrast in dom eye
                                & contrasts(:,2) == stimcontrast(c)... % contrast in null eye
                                & tilts(:,1) == X.dinull(2)... % null orientation in dom eye
                                & tilts(:,2) == X.dinull(2); % null orientation in null eye
                        end
                end
                
                if sum(trls) >= 5
                    binRESP(c,:,bin)    = nanmean(resp(:,trls),2);
                    binRESP(1,:,bin)    = nanmean(resp(4,:),2);
                    binRESPerror(c,:,bin)   = (nanstd(resp(:,trls),0,2))./(sqrt(sum(trls)));
                end
                
                binTrlNum(c,bin) = sum(trls);
                Trls.bin(c,bin,pen) = sum(trls);
            end
        end
        clear bin trls c
      
        
        %% SAVE UNIT in dMUA structure
        
        IDX(uct).penetration = STIM.penetration;
        IDX(uct).v1lim = STIM.v1lim;
        IDX(uct).depth = STIM.depths(e,:)';
        
        IDX(uct).prefeye    = DE;
        IDX(uct).prefori    = PS;
        IDX(uct).nulleye    = NDE;
        IDX(uct).nullori    = NS;
        IDX(uct).effects     = X.dianp; % p for main effect of each 'eye' 'tilt' 'contrast'
        
        IDX(uct).X      =   X;
        
        IDX(uct).occana       = X.occana;
        IDX(uct).oriana       = X.oriana;
        IDX(uct).diana        = X.diana;
        
        
        IDX(uct).occ   = X.occ';    % how much it prefers one eye over the other
        IDX(uct).ori   = X.ori';    % how much it prefers one orientation over the other
        IDX(uct).bio   = X.bio';    % How much it prefers both eyes over one
        
        IDX(uct).monLevels     = stimcontrast;
        IDX(uct).monTrials     = monTrlNum;
        IDX(uct).binTrials     = binTrlNum;
        IDX(uct).exactTrials   = [sum(monTrlNum(2:end,:),'all'),...
            sum(binTrlNum(2:end,:),'all')];
        IDX(uct).Total  = sum(IDX(uct).exactTrials,'all');
        
        %% Calculate Curves
        
        clear -global Data
        clear bin_all mon_all nBin nMon
        clear prd1 prd2 cv curves a K n b x
        
        % strings
        condition = {'mon','bin'};
        x = IDX(uct).monLevels*100; x(1) = x(1)+1;
        
        % data
        mon_all = squeeze(monRESP(:,:,1)); % only the DE_PS
        bin_all = squeeze(binRESP(:,:,1)); % only the BIN_PS
       
        
        global Data
        Data(1,:) = x;
        for tw = 1:2 % for transient and sustained
            
            mon = mon_all(:,tw);
            mon_ref = mon_all(:,1);
            bin = bin_all(:,tw);
            
            % normalize the response (on a unit-by-unit basis)
            mn        = min(mon_ref);
            mx        = max(mon_ref);
            nMon      = (mon - mn)./(mx - mn);
            nBin      = (bin - mn)./(mx - mn);
            
            for cond = 1:size(condition,2) % for each curve (bin and mon)
                switch condition{cond}
                    case 'mon'
                        Data(2,:) = nMon;
                    case 'bin'
                        Data(2,:) = nBin;
                end
                
                [a,K,n,~] = BMfitCRdata;
                predictions = 1:100;
                if cond == 1
                    if sum(isnan(mon_all)) == 0
                        for c = 1:length(predictions) % generate prediction for mon
                            MON(c,tw,uct) = a*[(c^n)/((c^n) + (K^n))]; % mon curve
                        end
                    else
                        MON(1:100,tw,uct) = NaN;
                        badCurves.mon = badCurves.mon + 0.5;
                    end
                elseif cond == 2
                    if sum(isnan(bin_all)) == 0
                        for c = 1:length(predictions) % generate prediction for bin
                            BIN(c,tw,uct) = a*[(c^n)/((c^n) + (K^n))]; % bin curve
                        end
                    else 
                        BIN(1:100,tw,uct) = NaN;
                        badCurves.bin = badCurves.bin + 0.5;
                        
                    end
                end  
            end
        end
        
        toc
    end
    
    clear I L lf bf mf layers STIM
end

clearvars -except sdfWin Trls MON BIN IDX flag_save N uct dataset badCurves
%% run additional analyses
% Ocularity analysis
clear occValues
for i = 1:length(IDX)
occValues(i,1) = IDX(i).occ(3);
end

% sort ocularity values by absolute distance from zero
[B,I] = sort(occValues,'ComparisonMethod','abs');

missing = ismissing(B);
tooLarge = abs(B) > 1.5;
I2 = I;
I2(missing | tooLarge) = [];

var3c = @(oldvar) mat2cell(oldvar(:), [fix(numel(oldvar)/3) *[1, 1], numel(oldvar)-2*fix(numel(oldvar)/3)], 1);     % Create New Matrix From Original Vector
occGroups = var3c(I2); % rows are: low med high

occLengths = [numel(occGroups{1,1}),numel(occGroups{2,1}),numel(occGroups{3,1})];

clear i var3c
%% Create OCC
clear o RESP *OCC*
for o = 1:size(occGroups,1)
    monRESP = nan(size(MON,1),size(MON,2),occLengths(o));
    binRESP = nan(size(BIN,1),size(BIN,2),occLengths(o));
    for u = 1:length(occGroups{o,1})
        monRESP(:,:,u) = MON(:,:,occGroups{o,1}(u));
        binRESP(:,:,u) = BIN(:,:,occGroups{o,1}(u));
    end
    OCC.MON(o).units = monRESP;
    OCC.BIN(o).units = binRESP;
end

for o = 1:3
    OCC.MON(o).avg = nanmean(OCC.MON(o).units,3);
    OCC.MON(o).err = nanstd(OCC.MON(o).units,[],3)./sqrt(size(OCC.MON(o).units,3));
    OCC.BIN(o).avg = nanmean(OCC.BIN(o).units,3);
    OCC.BIN(o).err = nanstd(OCC.BIN(o).units,[],3)./sqrt(size(OCC.BIN(o).units,3));
end

clear monRESP binRESP
%% Visualize 
% Figure save directory
if strcmp(getenv('username'),'bmitc')
    figDir = 'C:\Users\bmitc\OneDrive - Vanderbilt\Maier Lab\Figures\';
elseif strcmp(getenv('username'),'mitchba2')
    figDir = 'C:\Users\mitchba2\OneDrive - Vanderbilt\Maier Lab\Figures\';
elseif strcmp(getenv('username'),'bmitc_000')
    figDir = 'C:\Users\bmitc_000\OneDrive - Vanderbilt\Maier Lab\Figures\';
end

% toggle to save figures
flag_figsave = 1;

% Plot
tw = 2;
figure('position',[272,340,1119,375]);
for o = 1:3
    subplot(1,3,o)
    s1 = semilogx(1:100,OCC.MON(o).avg(:,tw),'color','k','linestyle','-','linewidth',2);  hold on
    s2 = semilogx(1:100,OCC.BIN(o).avg(:,tw),'color','b','linestyle','-','linewidth',2);  hold on
    
    if tw == 1
        set(s1,'linestyle','-'); set(s2,'linestyle','-')
        set(gca,'ylim',[0 1.2]);
    else
        set(s1,'linestyle','-'); set(s2,'linestyle','-')
        set(gca,'ylim',[0 .6]);
    end
    
    set(gca,'FontSize',12,'linewidth',2,'box','off',...
        'xlim',[1 100]);
    
    if o == 1
        occRange = abs([occValues(occGroups{o,1}(1)),occValues(occGroups{o,1}(end))]);
        title(sprintf('Lowest (n = %d)\n OccIndex [%.3g : %.3g]',occLengths(o),occRange(1),occRange(2)),'FontSize',16);
        xlabel('Stimulus contrast'); ylabel(sprintf('Relative response \n (normalized to monocular)'));
        xticklabels({'1','50','100'});
    elseif o == 2
        occRange = abs([occValues(occGroups{o,1}(1)),occValues(occGroups{o,1}(end))]);
        title(sprintf('Average (n = %d)\n OccIndex [%.3g : %.3g]',occLengths(o),occRange(1),occRange(2)),'FontSize',16);
        mCt = sum(Trls.mon(2:end,1,:),'all'); bCt = sum(Trls.bin(2:end,1,:),'all');
        legend(sprintf('MON (%d trials)',mCt),sprintf('BIN (%d trials)',bCt),'location','northwest'); legend boxoff
        xlabel([]); ylabel([]); xticklabels([]);
    else
        occRange = abs([occValues(occGroups{o,1}(1)),occValues(occGroups{o,1}(end))]);
        title(sprintf('Highest (n = %d)\n OccIndex [%.3g : %.3g]',occLengths(o),occRange(1),occRange(2)),'FontSize',16);
        xlabel([]); ylabel([]); xticklabels([]);
    end
end

if flag_figsave == 1
    cd(strcat(figDir,'ocularity\'));
    saveas(gcf, strcat('allUnits_approach1_trans', '.svg'));
    disp("Figure saved");
else
    disp("Figure was not saved");
end
%% SAVE
% Need to save workspace

if flag_save == true
    if strcmp(getenv('username'),'mitchba2')
        cd('C:/users/mitchba2/Documents/MATLAB/workspace/');
    elseif strcmp(getenv('username'),'bmitc')
        cd('C:/Users/bmitc/Documents/MATLAB/workspaces/');
    elseif strcmp(getenv('username'),'bmitc_000')
        cd('C:/Users/bmitc_000/Documents/MATLAB/workspaces/');
    end
save(sprintf('Dioptic_curves_%s',dataset),'IDX','MON','BIN','Trls','N','uct','dataset','badCurves');
fprintf('Workspace saved\n');
end

fprintf('Complete\n');
