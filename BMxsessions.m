%% BMxsessions
% Loading in all session varaibles. Align sessions by a common contact reference point: bottom of layer 4. %
% Perform various calculations (subtractions, model computations,
% t-scores). 

clear

%% Establish directory for sessions
myFolder = 'D:\mcosinteroc2\';  % Specify the folder where the files live.
%load('SessionDepths'); % found in dependencies folder

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
BOL4_2(i,:) = [max(tmp.STIM.laminae.gran)-1:-1:-(32-max(tmp.STIM.laminae.gran))];

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
%% Alignment across sessions
%% Attempt 1
sn = fieldnames(sessions);
pn = fieldnames(sessions.BIN.aMUA);
% DE & NDE

clear e f c
for e = 1:2  % for DE, NDE
    for f = 1:2:3 % for pc.all and pc.coll
        for c = 1:4 %for contrasts 0, 22, 45, 90
            temp = permute(squeeze(sessions.(sn{e}).aMUA.pc.(fn{f})(c,:,:,:)),[3 1 2]); 
            [sAVG.(sn{e}).aMUA.pc.(fn{f}).aligned(c,:,:), corticaldepth, ~] = laminarmean(temp,BOL4_2);
        end
    end
end
% Success!!!!!!!!!!!!!!!!

% BIN
clear p f c
for p = 1:length(pn) % for pc, pc_LSM, and pc_QSM
    for f = 1:2:3 % for 'all' and 'coll'
        clear c
        for c = 1:4 %for contrasts 0, 22, 45, 90
            temp = permute(squeeze(sessions.BIN.aMUA.(pn{p}).(fn{f})(c,:,:,:)),[3 1 2]); 
            [sAVG.BIN.aMUA.(pn{p}).(fn{f}).aligned(c,:,:), corticaldepth, ~] = laminarmean(temp,BOL4_2);
        end
    end
end

% DI
clear p f c
for p = 1:length(pn) % for pc, pc_LSM, and pc_QSM
    for f = 1:2:3 % for 'all' and 'coll'
        clear c
        for c = 1:6 %for contrasts 0, 22, 45, 90
            temp = permute(squeeze(sessions.DI.aMUA.(pn{p}).(fn{f})(c,:,:,:)),[3 1 2]); 
            [sAVG.DI.aMUA.(pn{p}).(fn{f}).aligned(c,:,:), corticaldepth, ~] = laminarmean(temp,BOL4_2);
        end
    end
end
clear sn pn temp p f c e


%% Save Workspace

cd('D:\')
save(sprintf('session-wide_newcode'),'sAVG','BOL4_2','corticaldepth','sessions');

cd('C:/users/bmitc/Documents/MATLAB/workspaces/'); 
save(sprintf('session-wide_newcode'),'sAVG','BOL4_2','corticaldepth','sessions');

fprintf('Workspace saved');