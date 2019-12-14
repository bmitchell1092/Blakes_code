%% BMxsessions
% Loading in all session varaibles. Align sessions by a common contact reference point: bottom of layer 4. %
% Perform various calculations (subtractions, model computations,
% t-scores). 

clear

%% Establish directory for sessions
myFolder = 'D:\mcosinteroc2\';  % Specify the folder where the files live.
load('SessionDepths'); % found in dependencies folder

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.mat'); 
matFiles = dir(filePattern);
for k = 1:length(matFiles)
baseFileName{k} = matFiles(k).name;
fullFileName{k} = fullfile(myFolder, baseFileName{k});
end

%% Sessions

clear i 
for i = 1: length(fullFileName)
tmp = load(fullFileName{i},'STIM');

% STIM.DE/NDE/BIN/DI
fn = fieldnames(tmp.STIM.DE.aMUA.pc);

    clear f
    for f = 1:length(fn)
    sessions.DE.aMUA.pc.(fn{f})(:,:,:,i) = tmp.STIM.DE.aMUA.pc.(fn{f})(:,:,:);
    sessions.NDE.aMUA.pc.(fn{f})(:,:,:,i) = tmp.STIM.NDE.aMUA.pc.(fn{f})(:,:,:);
    sessions.BIN.aMUA.pc.(fn{f})(:,:,:,i) = tmp.STIM.BIN.aMUA.pc.(fn{f})(:,:,:);
    sessions.BIN.aMUA.pc_LSM.(fn{f})(:,:,:,i) = tmp.STIM.BIN.aMUA.pc_LSM.(fn{f})(:,:,:);
    sessions.BIN.aMUA.pc_QSM.(fn{f})(:,:,:,i) = tmp.STIM.BIN.aMUA.pc_QSM.(fn{f})(:,:,:);
    sessions.DI.aMUA.pc.(fn{f})(:,:,:,i) = tmp.STIM.DI.aMUA.pc.(fn{f})(:,:,:);
    sessions.DI.aMUA.pc_LSM.(fn{f})(:,:,:,i) = tmp.STIM.DI.aMUA.pc_LSM.(fn{f})(:,:,:);
    sessions.DI.aMUA.pc_QSM.(fn{f})(:,:,:,i) = tmp.STIM.DI.aMUA.pc_QSM.(fn{f})(:,:,:);
    end
    
    % STIM.calc
    sn = fieldnames(tmp.STIM.calc.aMUA.pc.subtract);
    clear f
    for f = 1:length(sn)
        clear j 
        snn = fieldnames(tmp.STIM.calc.aMUA.pc.subtract.(sn{f}));
            for j = 1:length(snn)
            sessions.calc.aMUA.pc.subtract.(sn{f}).(snn{j})(:,:,:,i) = ...
                tmp.STIM.calc.aMUA.pc.subtract.(sn{f}).(snn{j})(:,:,:);
            end
    end 
end

%% session AVG
% Step 2: Average and STD

clear f 
for f = 1:length(fn)
    % DE
    sAVG.DE.aMUA.pc.(fn{f})(1).data = nanmean(sessions.DE.aMUA.pc.(fn{f}),4);
    sAVG.DE.aMUA.pc.(fn{f})(1).error = nanstd(sessions.DE.aMUA.pc.(fn{f}),0,4) ...
        / (sqrt(size(sessions.DE.aMUA.pc.(fn{f}),4)));
    % NDE
    sAVG.NDE.aMUA.pc.(fn{f})(1).data = nanmean(sessions.NDE.aMUA.pc.(fn{f}),4);
    sAVG.NDE.aMUA.pc.(fn{f})(1).error = nanstd(sessions.NDE.aMUA.pc.(fn{f}),0,4) ...
        / (sqrt(size(sessions.NDE.aMUA.pc.(fn{f}),4)));
    % BIN & Models
    sAVG.BIN.aMUA.pc.(fn{f})(1).data = nanmean(sessions.BIN.aMUA.pc.(fn{f}),4);
    sAVG.BIN.aMUA.pc.(fn{f})(1).error = nanstd(sessions.BIN.aMUA.pc.(fn{f}),0,4) ...
        / (sqrt(size(sessions.DE.aMUA.pc.(fn{f}),4)));
    sAVG.BIN.aMUA.pc_LSM.(fn{f})(1).data = nanmean(sessions.BIN.aMUA.pc_LSM.(fn{f}),4);
    sAVG.BIN.aMUA.pc_LSM.(fn{f})(1).error = nanstd(sessions.BIN.aMUA.pc_LSM.(fn{f}),0,4) ...
        / (sqrt(size(sessions.BIN.aMUA.pc_LSM.(fn{f}),4)));
    sAVG.BIN.aMUA.pc_QSM.(fn{f})(1).data = nanmean(sessions.BIN.aMUA.pc_QSM.(fn{f}),4);
    sAVG.BIN.aMUA.pc_QSM.(fn{f})(1).error = nanstd(sessions.BIN.aMUA.pc_QSM.(fn{f}),0,4) ...
        / (sqrt(size(sessions.BIN.aMUA.pc_QSM.(fn{f}),4)));
    
   % DI & Models
    sAVG.DI.aMUA.pc.(fn{f})(1).data = nanmean(sessions.DI.aMUA.pc.(fn{f}),4);
    sAVG.DI.aMUA.pc.(fn{f})(1).error = nanstd(sessions.DI.aMUA.pc.(fn{f}),0,4) ...
        / (sqrt(size(sessions.DE.aMUA.pc.(fn{f}),4)));
    sAVG.DI.aMUA.pc_LSM.(fn{f})(1).data = nanmean(sessions.DI.aMUA.pc_LSM.(fn{f}),4);
    sAVG.DI.aMUA.pc_LSM.(fn{f})(1).error = nanstd(sessions.DI.aMUA.pc_LSM.(fn{f}),0,4) ...
        / (sqrt(size(sessions.DI.aMUA.pc_LSM.(fn{f}),4)));
    sAVG.DI.aMUA.pc_QSM.(fn{f})(1).data = nanmean(sessions.DI.aMUA.pc_QSM.(fn{f}),4);
    sAVG.DI.aMUA.pc_QSM.(fn{f})(1).error = nanstd(sessions.DI.aMUA.pc_QSM.(fn{f}),0,4) ...
        / (sqrt(size(sessions.DI.aMUA.pc_QSM.(fn{f}),4)));
end

%%
%%
%%
%% Variable: Calc, Subtraction plots

% Step 1: load Variable
clear i tmp subtractionDE
for i = 1: length(fullFileName)
tmp = load(fullFileName{i},'STIM');   % loaded in like contrast,channel,trial
                                      % I have
                                      % contrast,samples,channel,trial
sessions.calc.subtractionDE.transient(:,:,i) = tmp.STIM.calc.contacts.subtractionDE.transient(2:4,:);
sessions.calc.subtractionDE.sustained(:,:,i) = tmp.STIM.calc.contacts.subtractionDE.sustained(2:4,:);
end
%% Attempt 1
clear f
% for f = 1:1
    clear c
    for c = 1:4
        temp = permute(squeeze(sessions.DE.aMUA.pc.coll(c,:,:,:)),[3 1 2]); 
        [sAVG.DE.aMUA.pc.coll.aligned(c,:,:), corticaldepth, ~] = laminarmean(temp,BOL4_test);
    end
% end
% Success!!!!!!!!!!!!!!!!
%%
% Step 2: Permute matrix to work with align function
sessions_subtractionDE_transient = permute(sessions.calc.subtractionDE.transient,[3 1 2]);
sessions_subtractionDE_sustained = permute(sessions.calc.subtractionDE.sustained,[3 1 2]);

% Step 2.5: Load in alignment matrix
load('SessionDepths'); % found in dependencies folder

% Step 3: Alignment
[AVG.calc.subtraction.transient, ~, ~] = laminarmean(sessions_subtractionDE_transient,BOL4);
[AVG.calc.subtraction.sustained, corticaldepth, N] = laminarmean(sessions_subtractionDE_sustained,BOL4);

clear sessions_subtractionDE_transient sessions_subtractionDE_sustained

%% Variable: Coll

%% Variable: CSD

% Step 2: Permute matrix to work with align function

perm.sessions_CSD = permute(sessions.CSD,[3 1 2]);

% Step 3: Alignment

[AVG.CSD, ~, ~] = laminarmean(perm.sessions_CSD,BOL4);

%% Model calculations (stats)

% Dioptic Coll
clear i p z

for i = 1:size(sessions.calc.models.data.BIN.transient,1)
[~,stats.QSM.BIN.coll.transient(i),~,sessions.stats.QSM.BIN.coll.transient(i)] = ttest(sessions.calc.models.data.BIN.transient(i,:,:),sessions.calc.models.QSM.BIN.transient(i,:,:)); 
[~,stats.QSM.BIN.coll.sustained(i),~,sessions.stats.QSM.BIN.coll.sustained(i)] = ttest(sessions.calc.models.data.BIN.sustained(i,:,:),sessions.calc.models.QSM.BIN.sustained(i,:,:));
[~,stats.LSM.BIN.coll.transient(i),~,sessions.stats.LSM.BIN.coll.transient(i)] = ttest(sessions.calc.models.data.BIN.transient(i,:,:),sessions.calc.models.LSM.BIN.transient(i,:,:)); 
[~,stats.LSM.BIN.coll.sustained(i),~,sessions.stats.LSM.BIN.coll.sustained(i)] = ttest(sessions.calc.models.data.BIN.sustained(i,:,:),sessions.calc.models.LSM.BIN.sustained(i,:,:));
end

% Dioptic Layers
for i = 1:size(sessions.BIN.layers.transient,1)
    for l = 1:size(sessions.BIN.layers.transient,2)
[~,stats.QSM.BIN.layers.transient(i,l,:),~,sessions.stats.QSM.BIN.layers.transient(i,l)] = ttest(sessions.BIN.layers.transient(i,l,:),sessions.QSM.BIN.layers.transient(i,l,:));
[~,stats.QSM.BIN.layers.sustained(i,l,:),~,sessions.stats.QSM.BIN.layers.sustained(i,l)] = ttest(sessions.BIN.layers.sustained(i,l,:),sessions.QSM.BIN.layers.sustained(i,l,:));
[~,stats.LSM.BIN.layers.transient(i,l,:),~,sessions.stats.LSM.BIN.layers.transient(i,l)] = ttest(sessions.BIN.layers.transient(i,l,:),sessions.LSM.BIN.layers.transient(i,l,:));
[~,stats.LSM.BIN.layers.sustained(i,l,:),~,sessions.stats.LSM.BIN.layers.sustained(i,l)] = ttest(sessions.BIN.layers.sustained(i,l,:),sessions.LSM.BIN.layers.sustained(i,l,:));
    end
end


%% PLOT: Binned Layer Bar Graphs (DE vs BIN)
% Transient

contrast = [0 .22 .45 .9];
figure('Position', [148,73,633,487]);
clear i
for i = 1:3
subplot(3,3,i)
bar(AVG.BIN.layers.transient.data(:,i),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar(AVG.BIN.layers.transient.data(:,i),AVG.BIN.layers.transient.error(:,i),'o','marker','none','color','k');
bar(AVG.DE.layers.transient.data(:,i),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
errorbar(AVG.DE.layers.transient.data(:,i),AVG.DE.layers.transient.error(:,i),'o','marker','none','color','k');
set(gca,'box','off');
ylim([-5 50]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent change');
if i == 1
    title('Supragranular');
elseif i == 2
    title('Granular');
else 
    title('Infragranular');
end
hold off
end

clear i
for i = 1:3
subplot(3,3,i+3)
bar((AVG.BIN.layers.transient.data(:,i)-AVG.DE.layers.transient.data(:,i)),0.8,'FaceColor',[0.3, .3, 0.3],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-5 20]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent difference');
hold off
end

for i = 1:3
subplot(3,3,i+6)
bar((AVG.BIN.layers.transient.data(:,i)-AVG.DE.layers.transient.data(:,i)) ...
    ./ (AVG.DE.layers.transient.data(:,i)),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-.2 1]);
xticklabels(contrast)
xlabel('contrast')
ylabel('fold change');
hold off
end

sgtitle(sprintf('Transient: Binocular vs monocular (DE) stimulation | %d sessions',length(fullFileName)),'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_layers_DEvsBIN_transient'), '-jpg', '-transparent');

%% PLOT: Binned Layer Bar Graphs Dichoptic (DE+NDE and BINvsQSMvsLSM)
% 2 = [.22], 3 = [.45], 4 = [.90]
cDE = 2;
cNDE = 3;
% 1 = DE22NDE45, 2 = DE22NDE90, 3 = DE45NDE22, 4 = DE45NDE90
% 5 = DE90NDE22, 6 = DE90NDE45
di = 1;

switch cDE
    case 2
        DELevel = '.22';
    case 3
        DELevel = '.45';
    case 4
        DELevel = '.90';
end

switch cNDE
    case 2
        NDELevel = '.22';
    case 3
        NDELevel = '.45';
    case 4
        NDELevel = '.90';
end

QSMdifields = fieldnames(QSM.DI.coll);
LSMdifields = fieldnames(LSM.DI.coll);
labels = {0, .22, .45, .9,[],0,.22,.45,.9}; format bank;
figure('Position', [148,73,1200,582]);
clear i
for i = 1:3
subplot(3,3,i)
bar([AVG.DE.layers.transient.data(cDE,i); NaN; AVG.DE.layers.sustained.data(cDE,i)],0.8,'grouped','FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar([AVG.DE.layers.transient.data(cDE,i); [0]; AVG.DE.layers.sustained.data(cDE,i)], [AVG.DE.layers.transient.error(cDE,i); [0]; AVG.DE.layers.sustained.error(cDE,i)],'o','marker','none','color','k');
bar([AVG.NDE.layers.transient.data(cNDE,i); NaN; AVG.NDE.layers.sustained.data(cNDE,i)],0.4,'grouped','FaceColor',[.22, 0.23, 0.22],'EdgeColor','k','LineWidth',0.8);
errorbar([AVG.NDE.layers.transient.data(cNDE,i); [0]; AVG.NDE.layers.sustained.data(cNDE,i)], [AVG.NDE.layers.transient.error(cNDE,i); [0]; AVG.NDE.layers.sustained.error(cNDE,i)],'o','marker','none','color','k');
set(gca,'box','off');
ylim([-5 60]);
xticklabels(labels);
xlabel('contrast')
ylabel('percent change');
if i == 1
    title('Supragranular');
elseif i == 2
    title('Granular');
else 
    title('Infragranular');
end
hold off
end

clear i
for i = 1:3
subplot(3,3,i+3)
bar([LSM.BIN.layers.transient(:,i);NaN;LSM.BIN.layers.sustained(:,i)],0.8,'FaceColor',[1, 1, 1],'linestyle','--','EdgeColor','k','LineWidth',0.8);
hold on
bar([AVG.BIN.layers.transient.data(:,i);NaN;AVG.BIN.layers.sustained.data(:,i)],0.6,'grouped','FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
errorbar([AVG.BIN.layers.transient.data(:,i); [0]; AVG.BIN.layers.sustained.data(:,i)], [AVG.BIN.layers.transient.error(:,i); 0; AVG.BIN.layers.sustained.error(:,1)],'o','marker','none','color','k');
bar([QSM.DI.layers.transient(:,i);NaN;QSM.DI.layers.sustained(:,i)],0.4,'FaceColor',[.85, .325, .098],'linestyle',':','EdgeColor','k','LineWidth',1);
set(gca,'box','off');
ylim([-5 60]);
xticklabels(labels);
xlabel('contrast')
ylabel('percent change');
hold off
end

clear i
for i = 1:3
subplot(3,3,i+6)
bar([LSM.BIN.layers.transient(:,i)-AVG.BIN.layers.transient.data(:,i);NaN;LSM.BIN.layers.sustained(:,i)-AVG.BIN.layers.sustained.data(:,i)],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',0.8,'linestyle','--');
hold on
bar([QSM.BIN.layers.transient(:,i)-AVG.BIN.layers.transient.data(:,i);NaN;QSM.BIN.layers.sustained(:,i)-AVG.BIN.layers.sustained.data(:,i)],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',1,'linestyle',':');
set(gca,'box','off');
ylim([-5 20]);
xticklabels(labels);
xlabel('contrast')
ylabel('percent difference');
hold off
end

sgtitle(sprintf('Quadratic Summation Model prediction for Binocular response\n Same contrast in both eyes | %d sessions',length(fullFileName)),'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_DEvsNDEvsModelvsBIN_T&S'), '-jpg', '-transparent');

%% PLOT: Binned Layer Bar Graphs (DE+NDE and BINvsQSMvsLSM)
% Transient and Sustained

labels = {0, .22, .45, .9,[],0,.22,.45,.9}; format bank;
figure('Position', [148,73,1200,582]);
clear i
for i = 1:3
subplot(3,3,i)
bar([AVG.DE.layers.transient.data(:,i); NaN; AVG.DE.layers.sustained.data(:,i)],0.8,'grouped','FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar([AVG.DE.layers.transient.data(:,i); [0]; AVG.DE.layers.sustained.data(:,i)], [AVG.DE.layers.transient.error(:,i); [0]; AVG.DE.layers.sustained.error(:,i)],'o','marker','none','color','k');
bar([AVG.NDE.layers.transient.data(:,i); NaN; AVG.NDE.layers.sustained.data(:,i)],0.4,'grouped','FaceColor',[.22, 0.23, 0.22],'EdgeColor','k','LineWidth',0.8);
errorbar([AVG.NDE.layers.transient.data(:,i); [0]; AVG.NDE.layers.sustained.data(:,i)], [AVG.NDE.layers.transient.error(:,i); [0]; AVG.NDE.layers.sustained.error(:,i)],'o','marker','none','color','k');
set(gca,'box','off');
ylim([-5 60]);
xticklabels(labels);
xlabel('contrast')
ylabel('percent change');
if i == 1
    title('Supragranular');
elseif i == 2
    title('Granular');
else 
    title('Infragranular');
end
hold off
end

clear i
for i = 1:3
subplot(3,3,i+3)
bar([LSM.BIN.layers.transient(:,i);NaN;LSM.BIN.layers.sustained(:,i)],0.8,'FaceColor',[1, 1, 1],'linestyle','--','EdgeColor','k','LineWidth',0.8);
hold on
bar([AVG.BIN.layers.transient.data(:,i);NaN;AVG.BIN.layers.sustained.data(:,i)],0.6,'grouped','FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
errorbar([AVG.BIN.layers.transient.data(:,i); [0]; AVG.BIN.layers.sustained.data(:,i)], [AVG.BIN.layers.transient.error(:,i); 0; AVG.BIN.layers.sustained.error(:,1)],'o','marker','none','color','k');
bar([QSM.BIN.layers.transient(:,i);NaN;QSM.BIN.layers.sustained(:,i)],0.4,'FaceColor',[.85, .325, .098],'linestyle',':','EdgeColor','k','LineWidth',1);
set(gca,'box','off');
ylim([-5 60]);
xticklabels(labels);
xlabel('contrast')
ylabel('percent change');
hold off
end

clear i
for i = 1:3
subplot(3,3,i+6)
bar([((LSM.BIN.layers.transient(:,i)-AVG.BIN.layers.transient.data(:,i))./(AVG.BIN.layers.transient.data(:,i)))*100;NaN;((LSM.BIN.layers.sustained(:,i)-AVG.BIN.layers.sustained.data(:,i))./(AVG.BIN.layers.sustained.data(:,i)))*100],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',0.8,'linestyle','--');
hold on
bar([((QSM.BIN.layers.transient(:,i)-AVG.BIN.layers.transient.data(:,i))./(AVG.BIN.layers.transient.data(:,i)))*100;NaN;((QSM.BIN.layers.sustained(:,i)-AVG.BIN.layers.sustained.data(:,i))./(AVG.BIN.layers.sustained.data(:,i)))*100],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',1,'linestyle',':');
set(gca,'box','off');
ylim([-15 70]);
xticklabels(labels);
xlabel('contrast')
ylabel('percent difference');
hold off
end

sgtitle(sprintf('Quadratic Summation Model prediction for Binocular response\n Same contrast in both eyes | %d sessions',length(fullFileName)),'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_DEvsNDEvsModelvsBIN_T&S'), '-jpg', '-transparent');

%% Just models vs data: Percent difference 
% 2 = [.22], 3 = [.45], 4 = [.90]

% decision
subtract = 0;
% contrast labels
labels = {'22|45', '22|90', '45|22', '45|90','90|22','90|45',[],'22|45','22|90','45|22','45|90','90|22','90|45'}; format bank;

figure('Position', [148,73,1200,582]);
clear i
for i = 1:3
subplot(2,3,i)
bar([LSM.DI.layers.transient(:,i);NaN;LSM.DI.layers.sustained(:,i)],0.8,'FaceColor',[1, 1, 1],'linestyle','--','EdgeColor','k','LineWidth',0.8);
hold on
bar([AVG.DI.layers.transient.data(:,i);NaN;AVG.DI.layers.sustained.data(:,i)],0.6,'grouped','FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
errorbar([AVG.DI.layers.transient.data(:,i); [0]; AVG.DI.layers.sustained.data(:,i)], [AVG.DI.layers.transient.error(:,i); 0; AVG.DI.layers.sustained.error(:,i)],'o','marker','none','color','k');
bar([QSM.DI.layers.transient(:,i);NaN;QSM.DI.layers.sustained(:,i)],0.4,'FaceColor',[.85, .325, .098],'linestyle',':','EdgeColor','k','LineWidth',1);
set(gca,'box','off');
ylim([-5 60]);
xticklabels(labels);
xlabel('contrast (DE | NDE)')
ylabel('percent change');
hold off
if i == 1
    title('Supragranular');
elseif i == 2
    title('Granular');
else 
    title('Infragranular');
end
xtickangle(45)
hold off
end

if subtract == true
    clear i
    for i = 1:3
    subplot(2,3,i+3)
    bar([LSM.DI.layers.transient(:,i)-AVG.DI.layers.transient.data(:,i);NaN;LSM.DI.layers.sustained(:,i)-AVG.DI.layers.sustained.data(:,i)],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',0.8,'linestyle','--');
    hold on
    bar([QSM.DI.layers.transient(:,i)-AVG.DI.layers.transient.data(:,i);NaN;QSM.DI.layers.sustained(:,i)-AVG.DI.layers.sustained.data(:,i)],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',1,'linestyle',':');
    set(gca,'box','off');
    ylim([-5 20]);
    xticklabels(labels);
    xlabel('contrast (DE | NDE)')
    ylabel('subtraction');
    xtickangle(45)
    hold off
    lgd = legend('LSM-BIN','QSM-BIN','location','northwest');
    lgd.FontSize = 6;
    end

else
    clear i
    for i = 1:3
    subplot(2,3,i+3)
    bar([((LSM.DI.layers.transient(:,i)-AVG.DI.layers.transient.data(:,i))./(AVG.DI.layers.transient.data(:,i)))*100;NaN;((LSM.DI.layers.sustained(:,i)-AVG.DI.layers.sustained.data(:,i))./(AVG.DI.layers.sustained.data(:,i)))*100],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',0.8,'linestyle','--');
    hold on
    bar([((QSM.DI.layers.transient(:,i)-AVG.DI.layers.transient.data(:,i))./(AVG.DI.layers.transient.data(:,i)))*100;NaN;((QSM.DI.layers.sustained(:,i)-AVG.DI.layers.sustained.data(:,i))./(AVG.DI.layers.sustained.data(:,i)))*100],0.8,'grouped','FaceColor',[1, 1, 1],'EdgeColor','k','LineWidth',1,'linestyle',':');
    set(gca,'box','off');
    ylim([-20 60]);
    xticklabels(labels);
    xlabel('contrast (DE | NDE)')
    ylabel('% difference from data');
    xtickangle(45)
    hold off
    lgd = legend('LSM','QSM','location','northwest');
    lgd.FontSize = 6;
    end
end

sgtitle(sprintf('Model predictions for Binocular response\n All Dichoptic presentations | Averaged across %d sessions',length(fullFileName)),'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_Dichoptic-all-subtract'), '-jpg', '-transparent');
%% PLOT: Collapsed lineplots (Same contrast, QSM prediction)
cIndex = 3;

switch cIndex
    case 2
        cLevel = '.22';
    case 3
        cLevel = '.45';
    case 4
        cLevel = '.90';
end

figure('position',[151,58.33333333333333,834.6666666666666,574.6666666666666]);
subplot(2,4,1)
plot(AVG.DE.coll.transient(cIndex,:),corticaldepth,'Color',[0, 0.4470, 0.7410]);
hold on 
plot(AVG.BIN.coll.transient(cIndex,:),corticaldepth,'Color', 'r');
grid on
xlim([-5 70]);
%hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title(sprintf('%s DE',cLevel));
legend('DE','BIN','Location','northeast','orientation','vertical');
hold off

subplot(2,4,2)
plot(AVG.NDE.coll.transient(cIndex,:),corticaldepth,'Color',[.3 .2 .3]);
hold on 
grid on
xlim([-5 70]);
hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title(sprintf('%s NDE',cLevel));
legend('NDE','Location','northeast','orientation','vertical');
hold off

subplot(2,4,3)
plot(AVG.BIN.coll.transient(cIndex,:),corticaldepth,'color','r');
hold on 
plot(QSM.BIN.coll.transient(cIndex,:),corticaldepth,'-.','color','k');
plot(LSM.BIN.coll.transient(cIndex,:),corticaldepth,'linewidth',.6,'linestyle','--','color',[.35 .4 .3]);
grid on
xlim([-5 70]);
hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('QSM vs Binocular response');
legend('BIN','QSM','LSM','Location','northeast','orientation','vertical');
hold off

h = subplot(2,4,4);
bAVG_iCSD = filterCSD(AVG.CSD')';
imagesc(-50:600,1:37,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(30); set(v,'color','k','linestyle','-.','linewidth',.5);
set(gca,'tickdir','in','ytick','');
climit = max(abs(get(gca,'CLim'))*.7);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
plot([100 100], ylim,'k','linestyle','-.','linewidth',0.5)
hline(28,':')
title('CSD (All Trials)')
xlabel('time (ms)')
clrbar = colorbar; clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',8,'VerticalAlignment','middle');
ylabel('<-- cortical depth -->');
hold off

subplot(2,4,5)
plot(AVG.DE.coll.sustained(cIndex,:),corticaldepth,'Color',[0, 0.4470, 0.7410]);
hold on 
grid on
xlim([-5 70]);
hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title(sprintf('%s DE',cLevel));
legend('DE','Location','northeast','orientation','horizontal');
hold off

subplot(2,4,6)
plot(AVG.NDE.coll.sustained(cIndex,:),corticaldepth,'Color',[.3 .2 .3]);
hold on 
grid on
xlim([-5 70]);
hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title(sprintf('%s NDE',cLevel));
legend('NDE','Location','northeast','orientation','vertical');
hold off

subplot(2,4,7)
plot(AVG.BIN.coll.sustained(cIndex,:),corticaldepth,'color','r');
hold on 
plot(QSM.BIN.coll.sustained(cIndex,:),corticaldepth,'-.','color','k');
plot(LSM.BIN.coll.sustained(cIndex,:),corticaldepth,'linewidth',.6,'linestyle','--','color',[.35 .4 .3]);
grid on
xlim([-5 70]);
hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('QSM vs Binocular response');
legend('BIN','QSM','LSM','Location','northeast','orientation','vertical');
hold off

h = subplot(2,4,8);
bAVG_iCSD = filterCSD(AVG.CSD')';
imagesc(-50:600,1:37,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(101); set(v,'color','k','linestyle','-.','linewidth',.5);
set(gca,'tickdir','in','ytick','');  
climit = max(abs(get(gca,'CLim'))*.7);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
plot([450 450], ylim,'k','linestyle','-.','linewidth',0.5)
hline(28,'-.')
title('CSD (All Trials)')
xlabel('time (ms)')
clrbar = colorbar; clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',8,'VerticalAlignment','middle');
ylabel('<-- cortical depth -->');
hold off

sgtitle(sprintf('Session-averaged (N = %d) Quadratic Summation Model prediction \n %s contrast in DE eye | %s contrast in NDE eye',length(fullFileName),cLevel,cLevel));

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_QSMcoll_%s_botheyes_1',cLevel), '-jpg', '-transparent');

%% PLOT: Collapsed lineplots (Dichoptic contrasts, QSM prediction)
% 2 = [.22], 3 = [.45], 4 = [.90]
cDE = 4;
cNDE = 3;
% 1 = DE22NDE45, 2 = DE22NDE90, 3 = DE45NDE22, 4 = DE45NDE90
% 5 = DE90NDE22, 6 = DE90NDE45
di = 6;

switch cDE
    case 2
        DELevel = '.22';
    case 3
        DELevel = '.45';
    case 4
        DELevel = '.90';
end

switch cNDE
    case 2
        NDELevel = '.22';
    case 3
        NDELevel = '.45';
    case 4
        NDELevel = '.90';
end

QSMdifields = fieldnames(QSM.DI.coll);
LSMdifields = fieldnames(LSM.DI.coll);

figure('position',[151,58.33333333333333,834.6666666666666,574.6666666666666]);
subplot(2,4,1)
plot(AVG.DE.coll.transient(cDE,:),corticaldepth,'Color',[0, 0.4470, 0.7410]);
hold on 
grid on
xlim([-5 70]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title(sprintf('%s DE',DELevel));
legend('DE','Location','northeast','orientation','horizontal');
hold off

subplot(2,4,2)
plot(AVG.NDE.coll.transient(cNDE,:),corticaldepth,'Color',[.3 .2 .3]);
hold on 
grid on
xlim([-5 70]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title(sprintf('%s NDE',NDELevel));
legend('NDE','Location','northeast','orientation','vertical');
hold off

subplot(2,4,3)
plot(AVG.DI.coll.transient(di,:),corticaldepth,'color','r');
hold on 
plot(QSM.DI.coll.(QSMdifields{di}).transient,corticaldepth,'-.','color','k');
plot(LSM.DI.coll.(LSMdifields{di}).transient,corticaldepth,'linewidth',.6,'linestyle','--','color',[.35 .4 .3]);
grid on
xlim([-5 70]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('QSM vs Binocular response');
legend('BIN','QSM','LSM','Location','northeast','orientation','vertical');
hold off

h = subplot(2,4,4);
bAVG_iCSD = filterCSD(AVG.CSD')';
imagesc(-50:600,1:37,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(30); set(v,'color','k','linestyle','-.','linewidth',.5);
set(gca,'tickdir','in','ytick','');
climit = max(abs(get(gca,'CLim'))*.7);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
plot([100 100], ylim,'k','linestyle','-.','linewidth',0.5)
hline(28,'-.')
title('CSD (All Trials)')
xlabel('time (ms)')
clrbar = colorbar; clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',8,'VerticalAlignment','middle');
ylabel('<-- cortical depth -->');
hold off

subplot(2,4,5)
plot(AVG.DE.coll.sustained(cDE,:),corticaldepth,'Color',[0, 0.4470, 0.7410]);
hold on 
grid on
xlim([-5 70]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title(sprintf('%s DE',DELevel));
legend('DE','Location','northeast','orientation','horizontal');
hold off

subplot(2,4,6)
plot(AVG.NDE.coll.sustained(cNDE,:),corticaldepth,'Color',[.3 .2 .3]);
hold on 
grid on
xlim([-5 70]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
legend('NDE','Location','northeast','orientation','vertical');
title(sprintf('%s NDE',NDELevel));
hold off

subplot(2,4,7)
plot(AVG.DI.coll.sustained(di,:),corticaldepth,'color','r');
hold on 
plot(QSM.DI.coll.(QSMdifields{di}).sustained,corticaldepth,'-.','color','k');
plot(LSM.DI.coll.(LSMdifields{di}).sustained,corticaldepth,'linewidth',.6,'linestyle','--','color',[.35 .4 .3]);
grid on
xlim([-5 70]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('QSM vs Binocular response');
legend('BIN','QSM','LSM','Location','northeast','orientation','vertical');
hold off

h = subplot(2,4,8);
bAVG_iCSD = filterCSD(AVG.CSD')';
imagesc(-50:600,1:37,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(101); set(v,'color','k','linestyle','-.','linewidth',.5);
set(gca,'tickdir','in','ytick','');  
climit = max(abs(get(gca,'CLim'))*.7);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
plot([450 450], ylim,'k','linestyle','-.','linewidth',0.5)
hline(28,'-.')
title('CSD (All Trials)')
xlabel('time (ms)')
clrbar = colorbar; clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',8,'VerticalAlignment','middle');
ylabel('<-- cortical depth -->');
hold off

sgtitle(sprintf('Session-averaged (N = %d) Quadratic Summation Model prediction \n %s contrast in DE eye | %s contrast in NDE eye',length(fullFileName),DELevel,NDELevel));

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_QSMcoll_%sDE%sNDE',DELevel,NDELevel), '-jpg', '-transparent');

%% PLOT: Collapsed lineplots (Same contrast, for presentation)
cIndex = 2;

switch cIndex
    case 2
        cLevel = '.22';
    case 3
        cLevel = '.45';
    case 4
        cLevel = '.90';
end

figure('position',[208.3333333333333,134.3333333333333,834.6666666666667,409.9999999999999]);
subplot(1,3,1)
plot(AVG.DE.coll.transient(2,:),corticaldepth,'Color',[0, 0.4470, 0.7410]);
hold on 
plot(AVG.BIN.coll.transient(2,:),corticaldepth,'Color', 'r');
grid on
xlim([-5 70]);
%hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('22% contrast');
legend('One eye','Both eyes','Location','northeast','orientation','vertical');
hold off

subplot(1,3,2)
plot(AVG.DE.coll.transient(3,:),corticaldepth,'Color',[0, 0.4470, 0.7410]);
hold on 
plot(AVG.BIN.coll.transient(3,:),corticaldepth,'Color','r');
grid on
xlim([-5 70]);
%hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('45% contrast');
legend('One eye','Both eyes','Location','northeast','orientation','vertical');
hold off

subplot(1,3,3)
plot(AVG.DE.coll.transient(4,:),corticaldepth,'Color',[0, 0.4470, 0.7410]);
hold on 
plot(AVG.BIN.coll.transient(4,:),corticaldepth,'Color','r');
grid on
xlim([-5 70]);
%hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('90% contrast');
legend('One eye','Both eyes','Location','northeast','orientation','vertical');
hold off

%sgtitle(sprintf('Session-averaged (N = %d) Quadratic Summation Model prediction \n %s contrast in DE eye | %s contrast in NDE eye',length(fullFileName),cLevel,cLevel));

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_pres_both-eyes'), '-jpg', '-transparent');


%%
h = subplot(2,4,4);
bAVG_iCSD = filterCSD(AVG.CSD')';
imagesc(-50:600,1:37,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(30); set(v,'color','k','linestyle','-.','linewidth',.5);
set(gca,'tickdir','in','ytick','');
climit = max(abs(get(gca,'CLim'))*.7);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
plot([100 100], ylim,'k','linestyle','-.','linewidth',0.5)
hline(28,':')
title('CSD (All Trials)')
xlabel('time (ms)')
clrbar = colorbar; clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',8,'VerticalAlignment','middle');
ylabel('<-- cortical depth -->');
hold off

subplot(2,4,5)
plot(AVG.DE.coll.sustained(cIndex,:),corticaldepth,'Color',[0, 0.4470, 0.7410]);
hold on 
grid on
xlim([-5 70]);
hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title(sprintf('%s DE',cLevel));
legend('DE','Location','northeast','orientation','horizontal');
hold off

subplot(2,4,6)
plot(AVG.NDE.coll.sustained(cIndex,:),corticaldepth,'Color',[.3 .2 .3]);
hold on 
grid on
xlim([-5 70]);
hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title(sprintf('%s NDE',cLevel));
legend('NDE','Location','northeast','orientation','vertical');
hold off

subplot(2,4,7)
plot(AVG.BIN.coll.sustained(cIndex,:),corticaldepth,'color','r');
hold on 
plot(QSM.BIN.coll.sustained(cIndex,:),corticaldepth,'-.','color','k');
plot(LSM.BIN.coll.sustained(cIndex,:),corticaldepth,'linewidth',.6,'linestyle','--','color',[.35 .4 .3]);
grid on
xlim([-5 70]);
hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('QSM vs Binocular response');
legend('BIN','QSM','LSM','Location','northeast','orientation','vertical');
hold off

h = subplot(2,4,8);
bAVG_iCSD = filterCSD(AVG.CSD')';
imagesc(-50:600,1:37,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); % this makes the red color the sinks and the blue color the sources (convention)
colorbar; v = vline(101); set(v,'color','k','linestyle','-.','linewidth',.5);
set(gca,'tickdir','in','ytick','');  
climit = max(abs(get(gca,'CLim'))*.7);
set(gca,'CLim',[-climit climit],'Box','off','TickDir','out')
plot([450 450], ylim,'k','linestyle','-.','linewidth',0.5)
hline(28,'-.')
title('CSD (All Trials)')
xlabel('time (ms)')
clrbar = colorbar; clrbar.Label.String = 'nA/mm^3'; 
set(clrbar.Label,'rotation',270,'fontsize',8,'VerticalAlignment','middle');
ylabel('<-- cortical depth -->');
hold off

sgtitle(sprintf('Session-averaged (N = %d) Quadratic Summation Model prediction \n %s contrast in DE eye | %s contrast in NDE eye',length(fullFileName),cLevel,cLevel));

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_QSMcoll_%s_botheyes_1',cLevel), '-jpg', '-transparent');

%% Save Workspace

cd('D:\')
save(sprintf('session-wide'),'AVG','QSM','BOL4','LSM','sessions','corticaldepth');

cd('C:/users/bmitc/Documents/MATLAB/workspaces/'); 
save(sprintf('session-wide'),'AVG','QSM','BOL4','LSM','sessions','corticaldepth','stats');

fprintf('Workspace saved');