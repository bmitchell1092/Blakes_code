function cBHV = concatBHV(filelist)

% combines BHV files for analysis
% can also be used to simply convert .bhv to matlab var by giving length(filelist) == 1
% returns empty if length(filelist) == 1 and file has no recorded trials
% fall/winter 2014
% MAC

if ~iscell(filelist)
    filelist = {filelist};
end

cBHV = []; ct = 0;

for i = 1:length(filelist)
    fullfile = filelist{i};
    [pathname,filename,ext] = fileparts(fullfile);
    if isempty(pathname)
        % assume file is in CD
        fullfile = [pwd filesep filename ext];
    end
    % load file
    if strcmp('.bhv',ext)
        % need to convert to mat
        bhv  = m_bhv_read(fullfile);
    else
        clear bhv
        load(fullfile)
    end
    % check size
    ntr = bhv.NumTrials;
    if ntr > 1
        ct = ct + 1;
        field = fieldnames(bhv);
        fidx = find(~cellfun(@isempty,strfind(field,'numtrials')));
        fieldtrln = NaN(length(field),1);
        for f = fidx:length(field)-5
            [m, n] = size(bhv.(field{f}));
            if m == ntr || n == ntr;
                fieldtrln(f) = ntr;
            elseif m < 2 && n < 2
                fieldtrln(f) = NaN;
            else
                fieldtrln(f) = max([m n]);
            end
        end
        if ~all(fieldtrln(~isnan(fieldtrln))  == ntr)
            % trials missed because of file crash
            ntr = min(fieldtrln);
            bhv.numtrials = ntr;
            
            for f = fidx:length(field)-5
                [m, n] = size(bhv.(field{f}));
                if m > 2 || n > 2
                    if m > n
                        bhv.(field{f}) = bhv.(field{f})(1:ntr,:);
                    elseif  n > m
                        bhv.(field{f}) = bhv.(field{f})(:,1:ntr);
                    end
                end
            end
        end
    else
        continue
    end
    
    if ct == 1
        % just make bhv = cbhv with some small additions
        field = fieldnames(bhv);
        for f = 1:length(field)
            if ischar(bhv.(field{f}))
                bhv.(field{f}) = {bhv.(field{f})};
            end
        end
        clear cbhv
        cBHV = bhv;
        cBHV.TrialsPerFiles(ct)=  ntr;
        cBHV.FileNumberByTrial =  repmat(ct,ntr,1);
        
    else
        % concatenate
        field = fieldnames(bhv);
        fidx = find(~cellfun(@isempty,strfind(field,'NumTrials')));
        
        % see if new bhv has any fields not in cBHV and add
        addfields = setdiff(field,fieldnames(cBHV));
        for a = 1:length(addfields);
            objclass = class(bhv.(addfields{a}));
            cBHV.(addfields{a}) = eval([objclass '.empty']);
        end
        % so the same for subfield "VariableChanges"
        addfields = setdiff(fieldnames(bhv.VariableChanges),fieldnames(cBHV.VariableChanges));
        for a = 1:length(addfields);
            for ii=1:i
                cBHV.VariableChanges(ii).(addfields{a}).trial =[];
                cBHV.VariableChanges(ii).(addfields{a}).value =[];
            end
        end
        addfields = setdiff(fieldnames(cBHV.VariableChanges),fieldnames(bhv.VariableChanges));
        for a = 1:length(addfields);
            bhv.VariableChanges.(addfields{a}).trial =[];
            bhv.VariableChanges.(addfields{a}).value =[];
        end
        
        % convert char to cell
        for f = 1:length(field)
            if ischar(bhv.(field{f})) || strcmp('BlockOrder',field{f})
                bhv.(field{f}) = {bhv.(field{f})};
            end
        end
        for f = 1:length(field)
            [m, n] = size(bhv.(field{f}));
            if m == ntr && f > fidx
                cBHV.(field{f}) = cat(1, cBHV.(field{f}),bhv.(field{f}));
            elseif n == ntr && f > fidx
                cBHV.(field{f}) = cat(2, cBHV.(field{f}),bhv.(field{f}));
            elseif m == 1 && n == 1
                cBHV.(field{f}) = cat(2, cBHV.(field{f}),bhv.(field{f}));
            elseif m == 1
                cBHV.(field{f}) = cat(2, cBHV.(field{f}),bhv.(field{f}));
            elseif n == 1
                cBHV.(field{f}) = cat(1, cBHV.(field{f}),bhv.(field{f}));
            else
                cBHV.(field{f}) = cat(2, {cBHV.(field{f})},{bhv.(field{f})});
            end
        end
        cBHV.TrialsPerFiles(ct)=  ntr;
        cBHV.FileNumberByTrial =  cat(1,cBHV.FileNumberByTrial, repmat(ct,ntr,1));
        
        
        
    end
end


