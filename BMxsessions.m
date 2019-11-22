%% BMxsessions
% Align sessions by a common contact reference point: bottom of layer 4. 

clear

%% Establish directory for sessions
myFolder = 'D:\mcosinteroc\';  % Specify the folder where the files live.

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

%% Variable: Layers

% Step 1: Load in variables
clear i 
for i = 1: length(fullFileName)
tmp = load(fullFileName{i},'STIM');
sessions.BIN.layers.full(:,:,i) = tmp.STIM.BIN.layers.full(:,:);
sessions.BIN.layers.transient(:,:,i) = tmp.STIM.BIN.layers.transient(:,:);
sessions.BIN.layers.sustained(:,:,i) = tmp.STIM.BIN.layers.sustained(:,:);
sessions.DE.layers.full(:,:,i) = tmp.STIM.DE.layers.full(:,:);
sessions.DE.layers.transient(:,:,i) = tmp.STIM.DE.layers.transient(:,:);
sessions.DE.layers.sustained(:,:,i) = tmp.STIM.DE.layers.sustained(:,:);
sessions.NDE.layers.full(:,:,i) = tmp.STIM.NDE.layers.full(:,:);
sessions.NDE.layers.transient(:,:,i) = tmp.STIM.NDE.layers.transient(:,:);
sessions.NDE.layers.sustained(:,:,i) = tmp.STIM.NDE.layers.sustained(:,:);
sessions.DI.layers.full(:,:,i) = tmp.STIM.DI.exclusive.layers.full(:,:);
sessions.DI.layers.transient(:,:,i) = tmp.STIM.DI.exclusive.layers.transient(:,:);
sessions.DI.layers.sustained(:,:,i) = tmp.STIM.DI.exclusive.layers.sustained(:,:);
end

% Step 2: Calculation (mean and std error)

% fieldnames for FOR loop
session_binfields = fieldnames(sessions.BIN.layers);
session_defields = fieldnames(sessions.DE.layers);
session_ndefields = fieldnames(sessions.NDE.layers);
session_difields = fieldnames(sessions.DI.layers);

clear i AVG
for i = 1:length(session_binfields)
    AVG.DE.layers.(session_defields{i}).data = nanmean(sessions.DE.layers.(session_defields{i}),3);
    AVG.DE.layers.(session_defields{i}).error = nanstd(sessions.DE.layers.(session_defields{i}),0,3) ...
        / (sqrt(size(sessions.DE.layers.(session_defields{i}),3)));
    AVG.BIN.layers.(session_binfields{i}).data = nanmean(sessions.BIN.layers.(session_binfields{i}),3);
    AVG.BIN.layers.(session_binfields{i}).error = nanstd(sessions.BIN.layers.(session_binfields{i}),0,3) ...
        / (sqrt(size(sessions.BIN.layers.(session_binfields{i}),3)));
    AVG.NDE.layers.(session_ndefields{i}).data = nanmean(sessions.NDE.layers.(session_ndefields{i}),3);
    AVG.NDE.layers.(session_ndefields{i}).error = nanstd(sessions.NDE.layers.(session_ndefields{i}),0,3) ...
        / (sqrt(size(sessions.NDE.layers.(session_ndefields{i}),3)));
    AVG.DI.layers.(session_difields{i}).data = nanmean(sessions.DI.layers.(session_difields{i}),3);
    AVG.DI.layers.(session_difields{i}).error = nanstd(sessions.DI.layers.(session_difields{i}),0,3) ...
        / (sqrt(size(sessions.DI.layers.(session_difields{i}),3)));
end

%% Variable: Calc, Subtraction plots

% Step 1: load Variable
clear i tmp subtractionDE
for i = 1: length(fullFileName)
tmp = load(fullFileName{i},'STIM');
sessions.calc.subtractionDE.transient(:,:,i) = tmp.STIM.calc.contacts.subtractionDE.transient(2:4,:);
sessions.calc.subtractionDE.sustained(:,:,i) = tmp.STIM.calc.contacts.subtractionDE.sustained(2:4,:);
end

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

% Step 1: Load variable
clear i tmp subtractionDE
for i = 1: length(fullFileName)
tmp = load(fullFileName{i},'STIM');
sessions.DE.coll.transient(:,:,i) = tmp.STIM.DE.coll.transient(:,:);
sessions.DE.coll.sustained(:,:,i) = tmp.STIM.DE.coll.sustained(:,:);
sessions.NDE.coll.transient(:,:,i) = tmp.STIM.NDE.coll.transient(:,:);
sessions.NDE.coll.sustained(:,:,i) = tmp.STIM.NDE.coll.sustained(:,:);
sessions.BIN.coll.transient(:,:,i) = tmp.STIM.BIN.coll.transient(:,:);
sessions.BIN.coll.sustained(:,:,i) = tmp.STIM.BIN.coll.sustained(:,:);
sessions.DI.coll.transient(:,:,i) = tmp.STIM.DI.exclusive.coll.transient(:,:);
sessions.DI.coll.sustained(:,:,i) = tmp.STIM.DI.exclusive.coll.sustained(:,:);
end

% Step 2: Permuate matrix to work with align function

perm.sessions_collDE_transient = permute(sessions.DE.coll.transient,[3 1 2]);
perm.sessions_collDE_sustained = permute(sessions.DE.coll.sustained,[3 1 2]);
perm.sessions_collNDE_transient = permute(sessions.NDE.coll.transient,[3 1 2]);
perm.sessions_collNDE_sustained = permute(sessions.NDE.coll.sustained,[3 1 2]);
perm.sessions_collBIN_transient = permute(sessions.BIN.coll.transient,[3 1 2]);
perm.sessions_collBIN_sustained = permute(sessions.BIN.coll.sustained,[3 1 2]);
perm.sessions_collDI_transient = permute(sessions.DI.coll.transient,[3 1 2]);
perm.sessions_collDI_sustained = permute(sessions.DI.coll.sustained,[3 1 2]);

% Step 3: Alignment

[AVG.DE.coll.transient, ~, ~] = laminarmean(perm.sessions_collDE_transient,BOL4);
[AVG.DE.coll.sustained, ~, ~] = laminarmean(perm.sessions_collDE_sustained,BOL4);
[AVG.NDE.coll.transient, ~, ~] = laminarmean(perm.sessions_collNDE_transient,BOL4);
[AVG.NDE.coll.sustained, ~, ~] = laminarmean(perm.sessions_collNDE_sustained,BOL4);
[AVG.BIN.coll.transient, ~, ~] = laminarmean(perm.sessions_collBIN_transient,BOL4);
[AVG.BIN.coll.sustained, ~, ~] = laminarmean(perm.sessions_collBIN_sustained,BOL4);
[AVG.DI.coll.transient, ~, ~] = laminarmean(perm.sessions_collDI_transient,BOL4);
[AVG.DI.coll.sustained, corticaldepth, N] = laminarmean(perm.sessions_collDI_sustained,BOL4);
clear perm

%% Variable: CSD

% Step 1: Load variable
clear i tmp subtractionDE
for i = 1: length(fullFileName)
tmp = load(fullFileName{i},'STIM');
sessions.CSD(:,:,i) = tmp.STIM.bsl.CSD(1:651,:);
end

% Step 2: Permute matrix to work with align function

perm.sessions_CSD = permute(sessions.CSD,[3 1 2]);

% Step 3: Alignment

[AVG.CSD, ~, ~] = laminarmean(perm.sessions_CSD,BOL4);


%% Model calculations
% Quadratic summation Model: sqrt((cL)^2 + (cR)^2)

QSM.BIN.layers.full = sqrt((AVG.NDE.layers.full.data(:,:).^2) + (AVG.DE.layers.full.data(:,:).^2));
QSM.BIN.layers.transient = sqrt((AVG.NDE.layers.transient.data(:,:).^2) + (AVG.DE.layers.transient.data(:,:).^2));
QSM.BIN.layers.sustained = sqrt((AVG.NDE.layers.sustained.data(:,:).^2) + (AVG.DE.layers.sustained.data(:,:).^2));
QSM.BIN.coll.transient = sqrt((AVG.NDE.coll.transient(:,:).^2) + (AVG.DE.coll.transient(:,:).^2));
QSM.BIN.coll.sustained = sqrt((AVG.NDE.coll.sustained(:,:).^2) + (AVG.DE.coll.sustained(:,:).^2));
QSM.DI.coll.DE22_NDE45.transient = sqrt((AVG.DE.coll.transient(2,:).^2) + (AVG.NDE.coll.transient(3,:).^2));
QSM.DI.coll.DE22_NDE45.sustained = sqrt((AVG.DE.coll.sustained(2,:).^2) + (AVG.NDE.coll.sustained(3,:).^2));
QSM.DI.coll.DE22_NDE90.transient = sqrt((AVG.DE.coll.transient(2,:).^2) + (AVG.NDE.coll.transient(4,:).^2));
QSM.DI.coll.DE22_NDE90.sustained = sqrt((AVG.DE.coll.sustained(2,:).^2) + (AVG.NDE.coll.sustained(4,:).^2));
QSM.DI.coll.DE45_NDE22.transient = sqrt((AVG.DE.coll.transient(3,:).^2) + (AVG.NDE.coll.transient(2,:).^2));
QSM.DI.coll.DE45_NDE22.sustained = sqrt((AVG.DE.coll.sustained(3,:).^2) + (AVG.NDE.coll.sustained(2,:).^2));
QSM.DI.coll.DE45_NDE90.transient = sqrt((AVG.DE.coll.transient(3,:).^2) + (AVG.NDE.coll.transient(4,:).^2));
QSM.DI.coll.DE45_NDE90.sustained = sqrt((AVG.DE.coll.sustained(3,:).^2) + (AVG.NDE.coll.sustained(4,:).^2));
QSM.DI.coll.DE90_NDE22.transient = sqrt((AVG.DE.coll.transient(4,:).^2) + (AVG.NDE.coll.transient(2,:).^2));
QSM.DI.coll.DE90_NDE22.sustained = sqrt((AVG.DE.coll.sustained(4,:).^2) + (AVG.NDE.coll.sustained(2,:).^2));
QSM.DI.coll.DE90_NDE45.transient = sqrt((AVG.DE.coll.transient(4,:).^2) + (AVG.NDE.coll.transient(3,:).^2));
QSM.DI.coll.DE90_NDE45.sustained = sqrt((AVG.DE.coll.sustained(4,:).^2) + (AVG.NDE.coll.sustained(3,:).^2));


% Linear summation model (LSM): (cL + cR) 
LSM.BIN.layers.full = AVG.NDE.layers.full.data(:,:) + AVG.DE.layers.full.data(:,:);
LSM.BIN.layers.transient = AVG.NDE.layers.transient.data(:,:) + AVG.DE.layers.transient.data(:,:);
LSM.BIN.layers.sustained = AVG.NDE.layers.sustained.data(:,:) + AVG.DE.layers.sustained.data(:,:);
LSM.BIN.coll.transient = AVG.NDE.coll.transient(:,:) + AVG.DE.coll.transient(:,:);
LSM.BIN.coll.sustained = AVG.NDE.coll.sustained(:,:) + AVG.DE.coll.sustained(:,:);
LSM.DI.coll.DE22_NDE45.transient = AVG.DE.coll.transient(2,:) + AVG.NDE.coll.transient(3,:);
LSM.DI.coll.DE22_NDE45.sustained = AVG.DE.coll.sustained(2,:) + AVG.NDE.coll.sustained(3,:);
LSM.DI.coll.DE22_NDE90.transient = AVG.DE.coll.transient(2,:) + AVG.NDE.coll.transient(4,:);
LSM.DI.coll.DE22_NDE90.sustained = AVG.DE.coll.sustained(2,:) + AVG.NDE.coll.sustained(4,:);
LSM.DI.coll.DE45_NDE22.transient = AVG.DE.coll.transient(3,:) + AVG.NDE.coll.transient(2,:);
LSM.DI.coll.DE45_NDE22.sustained = AVG.DE.coll.sustained(3,:) + AVG.NDE.coll.sustained(2,:);
LSM.DI.coll.DE45_NDE90.transient = AVG.DE.coll.transient(3,:) + AVG.NDE.coll.transient(4,:);
LSM.DI.coll.DE45_NDE90.sustained = AVG.DE.coll.sustained(3,:) + AVG.NDE.coll.sustained(4,:);
LSM.DI.coll.DE90_NDE22.transient = AVG.DE.coll.transient(4,:) + AVG.NDE.coll.transient(2,:);
LSM.DI.coll.DE90_NDE22.sustained = AVG.DE.coll.sustained(4,:) + AVG.NDE.coll.sustained(2,:);
LSM.DI.coll.DE90_NDE45.transient = AVG.DE.coll.transient(4,:) + AVG.NDE.coll.transient(3,:);
LSM.DI.coll.DE90_NDE45.sustained = AVG.DE.coll.sustained(4,:) + AVG.NDE.coll.sustained(3,:);

%% PLOT: Binned Layer Bar Graphs (DE vs BIN)
% Transient

contrast = [0 .22 .45 .9];
figure('Position', [148,73,633,487]);
clear i
for i = 1:3
subplot(3,3,i)
bar(AVG.BIN.layers.transient.data(2:4,i),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar(AVG.BIN.layers.transient.data(2:4,i),AVG.BIN.layers.transient.error(2:4,i),'o','marker','none','color','k');
%bar(AVG.DE.layers.transient.data(:,i),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
%errorbar(AVG.DE.layers.transient.data(:,i),AVG.DE.layers.transient.error(:,i),'o','marker','none','color','k');
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

%% PLOT: Binned Layer Bar Graphs (DE vs NDE)

contrast = [0 .22 .45 .9];
figure('Position', [148,73,633,487]);
clear i
for i = 1:3
subplot(3,3,i)
bar(AVG.DE.layers.transient.data(:,i),0.8,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar(AVG.DE.layers.transient.data(:,i),AVG.DE.layers.transient.error(:,i),'o','marker','none','color','k');
bar(AVG.NDE.layers.transient.data(:,i),0.4,'FaceColor',[.2, .2, .2],'EdgeColor','k','LineWidth',0.8);
errorbar(AVG.NDE.layers.transient.data(:,i),AVG.NDE.layers.transient.error(:,i),'o','marker','none','color','k');
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
bar((AVG.DE.layers.transient.data(:,i)-AVG.NDE.layers.transient.data(:,i)),0.8,'FaceColor',[0.3, .3, 0.3],'EdgeColor','k','LineWidth',0.8);
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
bar((AVG.DE.layers.transient.data(:,i)-AVG.NDE.layers.transient.data(:,i)) ...
    ./ (AVG.NDE.layers.transient.data(:,i)),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-.2 1]);
xticklabels(contrast)
xlabel('contrast')
ylabel('fold change');
hold off
end

sgtitle(sprintf('Transient: DE vs NDE stimulation | %d sessions',length(fullFileName)),'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_layers_DEvsNDE_transient'), '-jpg', '-transparent');

%% PLOT: Binned Layer Bar Graphs (DE+NDE and BIN+QSM)
% Transient and Sustained

labels = {0, .22, .45, .9,[],0,.22,.45,.9}; format bank;
figure('Position', [148,73,1200,582]);
clear i
for i = 1:3
subplot(2,3,i)
bar([AVG.DE.layers.transient.data(:,i); NaN; AVG.DE.layers.sustained.data(:,i)],0.8,'grouped','FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar([AVG.DE.layers.transient.data(:,i); [0]; AVG.DE.layers.sustained.data(:,i)], [AVG.DE.layers.transient.error(:,i); [0]; AVG.DE.layers.sustained.error(:,i)],'o','marker','none','color','k');
bar([AVG.NDE.layers.transient.data(:,i); NaN; AVG.NDE.layers.sustained.data(:,i)],0.4,'grouped','FaceColor',[.22, 0.23, 0.22],'EdgeColor','k','LineWidth',0.8);
errorbar([AVG.NDE.layers.transient.data(:,i); [0]; AVG.NDE.layers.sustained.data(:,i)], [AVG.NDE.layers.transient.error(:,i); [0]; AVG.NDE.layers.sustained.error(:,i)],'o','marker','none','color','k');
set(gca,'box','off');
ylim([-5 50]);
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
subplot(2,3,i+3)
bar([AVG.BIN.layers.transient.data(:,i);NaN;AVG.BIN.layers.sustained.data(:,i)],0.8,'grouped','FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar([AVG.BIN.layers.transient.data(:,i); [0]; AVG.BIN.layers.sustained.data(:,i)], [AVG.BIN.layers.transient.error(:,i); 0; AVG.BIN.layers.sustained.error(:,1)],'o','marker','none','color','k');
bar([QSM.BIN.layers.transient(:,i);NaN;QSM.BIN.layers.sustained(:,i)],0.4,'FaceColor',[0.9, .3, 0.1],'linestyle','-.','EdgeColor','w','LineWidth',0.8);
set(gca,'box','off');
ylim([-5 50]);
xticklabels(labels);
xlabel('contrast')
ylabel('percent change');
hold off
end

sgtitle(sprintf('Quadratic Summation Model prediction for Binocular response\n Same contrast in both eyes | %d sessions',length(fullFileName)),'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_DEvsNDEvsModelvsBIN_T&S'), '-jpg', '-transparent');

%% PLOT: Collapsed lineplot (.22 DE .90 NDE + QSM predictions)
% Transient and Sustained

figure('position',[151,58.33333333333333,834.6666666666666,574.6666666666666]);
subplot(2,4,1)
plot(AVG.DE.coll.transient(2,:),corticaldepth,'Color',[0, 0.4470, 0.7410]);
hold on 
grid on
xlim([-5 50]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('.22 DE');
legend('DE','Location','northeast','orientation','horizontal');
hold off

subplot(2,4,2)
plot(AVG.NDE.coll.transient(4,:),corticaldepth,'Color',[.3 .2 .3]);
hold on 
grid on
xlim([-5 50]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('.90 NDE');
legend('NDE','Location','northeast','orientation','vertical');
hold off

subplot(2,4,3)
plot(AVG.DI.coll.transient(2,:),corticaldepth,'color','r');
hold on 
plot(QSM.DI.coll.DE22_NDE90.transient,corticaldepth,'-.','color','k');
grid on
xlim([-5 50]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('QSM vs Binocular response');
legend('BIN','QSM','Location','northeast','orientation','vertical');
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
plot(AVG.DE.coll.sustained(2,:),corticaldepth,'Color',[0, 0.4470, 0.7410]);
hold on 
grid on
xlim([-5 50]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('.22 DE');
hold off

subplot(2,4,6)
plot(AVG.NDE.coll.sustained(4,:),corticaldepth,'Color',[.3 .2 .3]);
hold on 
grid on
xlim([-5 50]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('.90 NDE');
hold off

subplot(2,4,7)
plot(AVG.DI.coll.sustained(2,:),corticaldepth,'color','r');
hold on 
plot(QSM.DI.coll.DE22_NDE90.sustained,corticaldepth,'-.','color','k');
grid on
xlim([-5 50]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('QSM vs Binocular response');
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

sgtitle(sprintf('Session-averaged (N = %d) Quadratic Summation Model prediction \n Low contrast in DE | High contrast in NDE',length(fullFileName)));

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_coll_DE22_NDE90'), '-jpg', '-transparent');

%% DE .90 NDE .22 + QSM prediction
figure('position',[151,58.33333333333333,834.6666666666666,574.6666666666666]);

subplot(2,4,1)
plot(AVG.DE.coll.transient(4,:),corticaldepth,'Color',[0, 0.4470, 0.7410]);
hold on 
grid on
xlim([-5 50]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('.90 DE');
legend('DE','Location','northeast','orientation','horizontal');
hold off

subplot(2,4,2)
plot(AVG.NDE.coll.transient(2,:),corticaldepth,'Color',[.3 .2 .3]);
hold on 
grid on
xlim([-5 50]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('.22 NDE');
legend('NDE','Location','northeast','orientation','vertical');
hold off

subplot(2,4,3)
plot(AVG.DI.coll.transient(5,:),corticaldepth,'color','r');
hold on 
plot(QSM.DI.coll.DE90_NDE22.transient,corticaldepth,'-.','color','k');
grid on
xlim([-5 50]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('QSM vs Binocular response');
legend('BIN','QSM','Location','northeast','orientation','vertical');
hold off

h = subplot(2,4,4);
bAVG_iCSD = filterCSD(AVG.CSD')';
imagesc(-50:600,1:37,bAVG_iCSD');
hold on
colormap(flipud(colormap('jet'))); 
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
plot(AVG.DE.coll.sustained(4,:),corticaldepth,'Color',[0, 0.4470, 0.7410]);
hold on 
grid on
xlim([-5 50]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('.90 DE');
hold off

subplot(2,4,6)
plot(AVG.NDE.coll.sustained(2,:),corticaldepth,'Color',[.3 .2 .3]);
hold on 
grid on
xlim([-5 50]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('.22 NDE');
hold off

subplot(2,4,7)
plot(AVG.DI.coll.sustained(5,:),corticaldepth,'color','r');
hold on 
plot(QSM.DI.coll.DE90_NDE22.sustained,corticaldepth,'-.','color','k');
grid on
xlim([-5 50]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
title('QSM vs Binocular response');
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

sgtitle(sprintf('Session-averaged (N = %d) Quadratic Summation Model prediction \n High contrast in DE | Low contrast in NDE',length(fullFileName)));

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_coll_DE90_NDE22'), '-jpg', '-transparent');

%% PLOT: Collapsed lineplots (Same contrast, QSM prediction)
cIndex = 4;
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
grid on
xlim([-5 70]);
hline(0,':','BOL4')
ylim([-9 27])
xlabel('Percent change');
title(sprintf('%s DE',cLevel));
legend('DE','Location','northeast','orientation','horizontal');
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
export_fig(sprintf('AVG_QSMcoll_%s_botheyes',cLevel), '-jpg', '-transparent');

%% PLOT: Collapsed lineplots (Dichoptic contrasts, QSM prediction)
% 2 = [.22], 3 = [.45], 4 = [.90]
cDE = 3;
cNDE = 4;
% 1 = DE22NDE45, 2 = DE22NDE90, 3 = DE45NDE22, 4 = DE45NDE90
% 5 = DE90NDE22, 6 = DE90NDE45
di = 4;

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
legend('BIN','QSM','Location','northeast','orientation','vertical');
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
hold off

subplot(2,4,6)
plot(AVG.NDE.coll.sustained(cNDE,:),corticaldepth,'Color',[.3 .2 .3]);
hold on 
grid on
xlim([-5 70]);
hline(0,'-.','BOL4')
ylim([-9 27])
xlabel('Percent change');
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


%% Save Workspace

cd('D:\')
save(sprintf('session-wide'),'AVG','QSM','BOL4','LSM','sessions','corticaldepth');

cd('C:/users/bmitc/Documents/MATLAB/workspaces/'); 
save(sprintf('session-wide'),'AVG','QSM','BOL4','LSM','sessions','corticaldepth');

fprintf('Workspace saved');