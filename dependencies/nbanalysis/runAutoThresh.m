clear all
warning off

if ~ispc && ~exist('/users/kaciedougherty/documents/code')
    addpath('/volumes/drobo/users/kacie/code/processdata_code');
    % addpath('/volumes/drobo/users/kacie/code/mlanalysis');
    addpath('/volumes/drobo/users/kacie/code/nbanalysis');
    addpath(genpath('/volumes/drobo/lab software/neurophys analysis'));
elseif ~ispc
    addpath(genpath('/users/kaciedougherty/documents/code/nbanalysis'))
    addpath('/users/kaciedougherty/documents/neurophysdata')
    addpath(genpath('/users/kaciedougherty/documents/code/BLACKROCK'))
    addpath(genpath(('/Users/kaciedougherty/Documents/Code/KiloSortUtils')))
else
    addpath('/users/MLab/documents');
    addpath('/users/MLab/documents/mlanalysisonline');
    addpath('/users/MLab/documents/utils');
end
%%

totchan  = 24;
[~,filelist] =  getUnitInfo('160609_I_cinterocdrft013_11');
ftype       = 'bw';
if strcmp(ftype,'bw')
    for f = 1:length(filelist)
        
        [unit] = getUnitInfo(filelist{f});
        if isfield(unit,'bw')
            bw_filelist{f} = strcat(unit.bw,filelist{f}(end-2:end));
        else
            bw_filelist{f} = '';
        end
    end
    clear filelist;
    filelist = bw_filelist;
    
end
[~,cinteroclist] = getUnitInfo('160609_I_cinterocdrft013_11');


for f = 17:length(filelist)
    if isempty(filelist{f})
        continue 
    else
    fprintf('file %u of %u',f,length(filelist));
    
    if strcmp(ftype,'cinteroc')
        BRdatafile = filelist{f}(1:24);
    elseif strcmp(ftype,'bw')
        BRdatafile = filelist{f}(1:21);
    end
    
    brdrname   = sprintf('/users/kaciedougherty/documents/neurophysdata/%s',BRdatafile(1:8));
    filename = strcat(brdrname,'/',BRdatafile,'.ns6');
    
    if ~exist(filename)
        filename = strcat('/volumes/Toshiba External/',BRdatafile(1:8),'/',BRdatafile,'.ns6');
    end
    
    if ~exist(filename)
        continue
    else
        
        [ppNEV]  = offlineBRAutoSort(filename(1:end-4));
        
        fprintf('saving file %u of %u',f,length(filelist));
        save(strcat(filename(1:end-4),'.ppNEV'),'ppNEV','-v7.3');
    end
    end
    
end
