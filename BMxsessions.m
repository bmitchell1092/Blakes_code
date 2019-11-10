%% BMxsessions
% Align sessions by a common contact reference point: bottom of layer 4. 

clear

%% Step 1: Load in variables from sessions
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

% Variable: Binned Layers
clear i 

for i = 1: length(fullFileName)
tmp = load(fullFileName{i},'STIM');
sessions.BIN.layers.full(:,:,i) = tmp.STIM.BIN.layers.full(2:4,:);
sessions.BIN.layers.transient(:,:,i) = tmp.STIM.BIN.layers.transient(2:4,:);
sessions.BIN.layers.sustained(:,:,i) = tmp.STIM.BIN.layers.sustained(2:4,:);
sessions.DE.layers.full(:,:,i) = tmp.STIM.DE.layers.full(2:4,:);
sessions.DE.layers.transient(:,:,i) = tmp.STIM.DE.layers.transient(2:4,:);
sessions.DE.layers.sustained(:,:,i) = tmp.STIM.DE.layers.sustained(2:4,:);
end

% Calculation: mean and error
session_binfields = fieldnames(sessions.BIN.layers);
session_defields = fieldnames(sessions.DE.layers);

clear i AVG
for i = 1:length(session_binfields)
    AVG.DE.layers.(session_defields{i}).data = mean(sessions.DE.layers.(session_defields{i}),3);
    AVG.DE.layers.(session_defields{i}).error = std(sessions.DE.layers.(session_defields{i}),0,3) ...
        / (sqrt(size(sessions.DE.layers.(session_defields{i}),3)));
    AVG.BIN.layers.(session_binfields{i}).data = mean(sessions.BIN.layers.(session_binfields{i}),3);
    AVG.BIN.layers.(session_binfields{i}).error = std(sessions.BIN.layers.(session_binfields{i}),0,3) ...
        / (sqrt(size(sessions.BIN.layers.(session_binfields{i}),3)));
end
    
% AVG.BIN.layers.full.data = mean(sessions.BIN.layers.full,3);
% AVG.BIN.layers.full.error = std(sessions.BIN.layers.full,0,3) / (sqrt(size(sessions.BIN.layers.full,3)));
% AVG.BIN.layers.transient.data = mean(sessions.BIN.layers.transient,3);
% AVG.BIN.layers.transient.error = std(sessions.BIN.layers.transient,0,3) / (sqrt(size(sessions.BIN.layers.transient,3)));
% AVG.BIN.layers.sustained.data = mean(sessions.BIN.layers.sustained,3);
% AVG.BIN.layers.sustained.error = std(sessions.BIN.layers.sustained,0,3) / (sqrt(size(sessions.BIN.layers.sustained,3)));
% 
% AVG.DE.layers.full.data = mean(sessions.DE.layers.full,3);
% AVG.DE.layers.full.error = std(sessions.DE.layers.full,0,3) / (sqrt(size(sessions.DE.layers.full,3)));
% AVG.DE.layers.transient.data = mean(sessions.DE.layers.transient,3);
% AVG.DE.layers.transient.error = std(sessions.DE.layers.transient,0,3)/(sqrt(size(sessions.DE.layers.transient,3)));
% AVG.DE.layers.sustained.data = mean(sessions.DE.layers.sustained,3);
% AVG.DE.layers.sustained.error = std(sessions.DE.layers.sustained,0,3)/(sqrt(size(sessions.DE.layers.sustained,3)));


%% PLOT (Averaged Bar Graphs)
% Transient
contrast = [.22 .45 .9];
figure('Position', [148,73,633,487]);
subplot(3,3,1)
bar(AVG.BIN.layers.transient.data(:,1),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar(AVG.BIN.layers.transient.data(:,1),AVG.BIN.layers.transient.error(:,1),'o','marker','none','color','k');
bar(AVG.DE.layers.transient.data(:,1),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
errorbar(AVG.DE.layers.transient.data(:,1),AVG.DE.layers.transient.error(:,1),'o','marker','none','color','k');
set(gca,'box','off');
ylim([-5 50]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent change');
title('Supragranular');
hold off

subplot(3,3,2)
bar(AVG.BIN.layers.transient.data(:,2),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar(AVG.BIN.layers.transient.data(:,2),AVG.BIN.layers.transient.error(:,2),'o','marker','none','color','k');
bar(AVG.DE.layers.transient.data(:,2),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
errorbar(AVG.DE.layers.transient.data(:,2),AVG.DE.layers.transient.error(:,2),'o','marker','none','color','k');
set(gca,'box','off');
ylim([-5 50]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent change');
title('Granular');
hold off

subplot(3,3,3)
bar(AVG.BIN.layers.transient.data(:,3),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar(AVG.BIN.layers.transient.data(:,3),AVG.BIN.layers.transient.error(:,3),'o','marker','none','color','k');
bar(AVG.DE.layers.transient.data(:,3),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
errorbar(AVG.DE.layers.transient.data(:,3),AVG.DE.layers.transient.error(:,3),'o','marker','none','color','k');
set(gca,'box','off');
ylim([-5 50]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent change');
title('Infragranular');
hold off

subplot(3,3,4)
bar((AVG.BIN.layers.transient.data(:,1)-AVG.DE.layers.transient.data(:,1)),0.8,'FaceColor',[0.3, .3, 0.3],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-5 20]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent difference');
%title('Supragranular');
hold off

subplot(3,3,5)
bar((AVG.BIN.layers.transient.data(:,2)-AVG.DE.layers.transient.data(:,2)),0.8,'FaceColor',[0.3, .3, 0.3],'EdgeColor','k','LineWidth',0.8);hold on
set(gca,'box','off');
ylim([-5 20]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent difference');
%title('Granular');
hold off

subplot(3,3,6)
bar((AVG.BIN.layers.transient.data(:,3)-AVG.DE.layers.transient.data(:,3)),0.8,'FaceColor',[0.3, .3, 0.3],'EdgeColor','k','LineWidth',0.8);hold on
set(gca,'box','off');
ylim([-5 20]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent difference');
hold off

subplot(3,3,7)
bar((AVG.BIN.layers.transient.data(:,1)-AVG.DE.layers.transient.data(:,1)) ...
    ./ (AVG.DE.layers.transient.data(:,1)),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-.2 1]);
xticklabels(contrast)
xlabel('contrast')
ylabel('fold change');
%title('Supragranular');
hold off

subplot(3,3,8)
bar((AVG.BIN.layers.transient.data(:,2)-AVG.DE.layers.transient.data(:,2)) ...
    ./ (AVG.DE.layers.transient.data(:,2)),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-.2 1]);
xticklabels(contrast)
xlabel('contrast')
ylabel('fold change');
%title('Granular');
hold off

subplot(3,3,9)
bar((AVG.BIN.layers.transient.data(:,3)-AVG.DE.layers.transient.data(:,3)) ...
    ./ (AVG.DE.layers.transient.data(:,3)),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-.2 1]);
xticklabels(contrast)
xlabel('contrast')
ylabel('fold change');
%title('Infragranular');
hold off

sgtitle(sprintf('Transient: Binocular vs monocular stimulation | %d sessions',length(fullFileName)),'Interpreter','none');

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_binned-layers_transient-new'), '-jpg', '-transparent');

%% PLOT (Averaged Bar Graphs)
% Sustained response

figure('Position', [148,73,633,487]);
subplot(3,3,1)
bar(AVG.BIN.layers.sustained.data(:,1),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar(AVG.BIN.layers.sustained.data(:,1),AVG.BIN.layers.sustained.error(:,1),'o','marker','none','color','k');
bar(AVG.DE.layers.sustained.data(:,1),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
errorbar(AVG.DE.layers.sustained.data(:,1),AVG.DE.layers.sustained.error(:,1),'o','marker','none','color','k');
set(gca,'box','off');
ylim([-5 50]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent change');
title('Supragranular');
hold off

subplot(3,3,2)
bar(AVG.BIN.layers.sustained.data(:,2),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar(AVG.BIN.layers.sustained.data(:,2),AVG.BIN.layers.sustained.error(:,2),'o','marker','none','color','k');
bar(AVG.DE.layers.sustained.data(:,2),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
errorbar(AVG.DE.layers.sustained.data(:,2),AVG.DE.layers.sustained.error(:,2),'o','marker','none','color','k');
set(gca,'box','off');
ylim([-5 50]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent change');
title('Granular');
hold off

subplot(3,3,3)
bar(AVG.BIN.layers.sustained.data(:,3),0.8,'FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar(AVG.BIN.layers.sustained.data(:,3),AVG.BIN.layers.sustained.error(:,3),'o','marker','none','color','k');
bar(AVG.DE.layers.sustained.data(:,3),0.4,'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
errorbar(AVG.DE.layers.sustained.data(:,3),AVG.DE.layers.sustained.error(:,3),'o','marker','none','color','k');
set(gca,'box','off');
ylim([-5 50]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent change');
title('Infragranular');
hold off

subplot(3,3,4)
bar((AVG.BIN.layers.sustained.data(:,1)-AVG.DE.layers.sustained.data(:,1)),0.8,'FaceColor',[0.3, .3, 0.3],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-5 20]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent difference');
%title('Supragranular');
hold off

subplot(3,3,5)
bar((AVG.BIN.layers.sustained.data(:,2)-AVG.DE.layers.sustained.data(:,2)),0.8,'FaceColor',[0.3, .3, 0.3],'EdgeColor','k','LineWidth',0.8);hold on
set(gca,'box','off');
ylim([-5 20]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent difference');
%title('Granular');
hold off

subplot(3,3,6)
bar((AVG.BIN.layers.sustained.data(:,3)-AVG.DE.layers.sustained.data(:,3)),0.8,'FaceColor',[0.3, .3, 0.3],'EdgeColor','k','LineWidth',0.8);hold on
set(gca,'box','off');
ylim([-5 20]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent difference');
hold off

subplot(3,3,7)
bar((AVG.BIN.layers.sustained.data(:,1)-AVG.DE.layers.sustained.data(:,1)) ...
    ./ (AVG.DE.layers.sustained.data(:,1)),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-.2 1]);
xticklabels(contrast)
xlabel('contrast')
ylabel('fold change');
%title('Supragranular');
hold off

subplot(3,3,8)
bar((AVG.BIN.layers.sustained.data(:,2)-AVG.DE.layers.sustained.data(:,2)) ...
    ./ (AVG.DE.layers.sustained.data(:,2)),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-.2 1]);
xticklabels(contrast)
xlabel('contrast')
ylabel('fold change');
%title('Granular');
hold off

subplot(3,3,9)
bar((AVG.BIN.layers.sustained.data(:,3)-AVG.DE.layers.sustained.data(:,3)) ...
    ./ (AVG.DE.layers.sustained.data(:,3)),0.8,'FaceColor',[0.7, 0.7, 0.7],'EdgeColor','k','LineWidth',0.8);
hold on
set(gca,'box','off');
ylim([-.2 1]);
xticklabels(contrast)
xlabel('contrast')
ylabel('fold change');
%title('Infragranular');
hold off

sgtitle(sprintf('Sustained: Binocular vs monocular stimulation | %d sessions',length(fullFileName)),'Interpreter','none');


cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_binned-layers_sustained'), '-jpg', '-transparent');

%% Alignment of collapsed contacts
%% Step 1: Load variable

% Variable: Collapsed subtractionDE
clear i tmp subtractionDE
for i = 1: length(fullFileName)
tmp = load(fullFileName{i},'STIM');
sessions.calc.subtractionDE.transient(:,:,i) = tmp.STIM.calc.contacts.subtractionDE.transient(2:4,:);
sessions.calc.subtractionDE.sustained(:,:,i) = tmp.STIM.calc.contacts.subtractionDE.sustained(2:4,:);
end

channels = size(sessions.calc.subtractionDE.transient,3);

% Step 2: Permuate matrix to work with align function
sessions_subtractionDE_transient = permute(sessions.calc.subtractionDE.transient,[3 1 2]);
sessions_subtractionDE_sustained = permute(sessions.calc.subtractionDE.sustained,[3 1 2]);

%% Alignment
% Step 3: load depth matrix (pre-defined)
depth = nan(length(fullFileName),size(sessions_subtractionDE_transient,3));

% Step 4: Alignment to bottom of layer 4
[lDAT_transient, ~, ~] = laminarmean(sessions_subtractionDE_transient,BOL4);

% aligned to bottom of layer 4:
[lDAT_sustained, corticaldepth, N] = laminarmean(sessions_subtractionDE_sustained,BOL4);


%% Plotting aligned data

figure;
subplot(1,2,1)
plot(lDAT_transient,corticaldepth);
hold on 
grid on
xlim([-5 20]);
hline(0,'-.','BOL4')
set(gca,'box','off');
%yticks(1:length(lDAT_transient))
% yticklabels(fliplr(1:nct))
% ylim([1 length(lDAT_transient)])
xlabel('Percent change');
title('Transient: Aligned');
%legend(num2str(contrast),'Location','southoutside','orientation','horizontal');
hold off

subplot(1,2,2)
plot(lDAT_sustained,corticaldepth);
hold on 
grid on
xlim([-5 20]);
hline(0,'-.','BOL4')
set(gca,'box','off');
% yticks(1:nct)
% yticklabels(fliplr(1:nct))
% ylim([1 nct])
xlabel('Percent change');
title('Sustained: Aligned');
hold off

%legend({num2str(contrast(:))},'Location','northoutside','orientation','vertical');

sgtitle(sprintf('Session-averaged (N = %d) subtraction plot (BIN - DE)\nTransient vs Sustained',length(fullFileName)));

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_subtraction-plots'), '-jpg', '-transparent');



%%
figure;

not_aligned_trans = permute(sessions_subtractionDE_transient,[2 3 1]);
%not_aligned_trans = permute(not_aligned_trans,[3 2 1]);
subplot(1,1,1)
hold on
clear i
for i=1:size(not_aligned_trans,3)
plot(fliplr(not_aligned_trans(:,:,i)),1:32)
end
grid on
xlim([-15 55]);
set(gca,'box','off');
%yticks(1:32)
yticklabels(fliplr(0:5:35))
%ylim([1 32])
xlabel('Percent change');
title('Transient: All sessions, NOT ALIGNED');
%legend(num2str(contrast),'Location','southoutside','orientation','horizontal');
hold off

cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
export_fig(sprintf('AVG_subtraction-plots-test-all'), '-jpg', '-transparent');
%%
subplot(1,2,2)
plot(lDAT_transient,corticaldepth);
hold on 
grid on
xlim([-5 20]);
hline(0,'-.')
set(gca,'box','off');
% yticks(1:length(lDAT_transient))
% yticklabels(fliplr(1:nct))
% ylim([1 length(lDAT_transient)])
xlabel('Percent change');
title('Transient: Aligned');
hold off

legend({num2str(contrast(:))},'Location','northoutside','orientation','horizontal');

sgtitle('Session-averaged subtraction plot (Binocular - monocular)');

% cd('C:\Users\bmitc\OneDrive\4. Vanderbilt\Maier Lab\Figures\')
% export_fig(sprintf('AVG_subtraction-plots-test'), '-jpg', '-transparent');


%% Experimental

figure('Position', [148,73,1200,582]);
subplot(3,3,1)
hBar = bar([AVG.BIN.layers.transient.data(:,1) AVG.BIN.layers.sustained.data(:,1)],0.8,'grouped','FaceColor',[0.8500, 0.3250, 0.0980],'EdgeColor','k','LineWidth',0.8);
hold on
errorbar([AVG.BIN.layers.transient.data(:,1) AVG.BIN.layers.sustained.data(:,1)], [AVG.BIN.layers.transient.error(:,1) AVG.BIN.layers.sustained.error(:,1)],'o','marker','none','color','k');
bar([AVG.DE.layers.transient.data(:,1) AVG.DE.layers.sustained.data(:,1)],0.4,'grouped','FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.8);
errorbar([AVG.DE.layers.transient.data(:,1) AVG.DE.layers.sustained.data(:,1)], [AVG.DE.layers.transient.error(:,1) AVG.DE.layers.sustained.error(:,1)],'o','marker','none','color','k');
set(gca,'box','off');
ylim([-5 50]);
xticklabels(contrast)
xlabel('contrast')
ylabel('percent change');
title('Supragranular');
hold off