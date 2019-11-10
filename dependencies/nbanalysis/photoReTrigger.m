function [newTP,trigger,STIM] = photoReTrigger(TP,filename,ypos,trigger,STIM)

% TP = obs x [st en] in SAMPELS
% filename = full address, no extension, must contain path to BHV and NS6
% ypos = location of stimulus on screen, can be scaler or same length at TP
% Aug 4 2017; MAC

photolabel = 'ainp1';
nr         = 1; % see note at line 111
newTP      = nan(size(TP));

if nargin < 4
    trigger = [];
end
if nargin < 3
    ypos = 0;
end

if length(ypos) == 1
    x=ypos;
    ypos = zeros(length(TP),1);
    ypos(:) = x; clear x
end
ypos(isnan(ypos))=0;

% correct for certain bhv file issues
[brpath,BRdatafile,~] = fileparts(filename);

if strfind('drobo',brpath)
    
    
    switch BRdatafile
        case {'151205_E_rfori002','151205_E_rfori003','151205_E_cosinteroc001'}
            bhvfile =  '/Volumes/Drobo/Data/NEUROPHYS/rig022/151205_E/151205_E_rfori005.bhv';
        case '161005_E_rfori001'
            bhvfile = '/Volumes/Drobo2/Data/NEUROPHYS/rig022/161005_E/161005_E_rfori002.bhv';
        case '161006_E_evp001'
            bhvfile = '/Volumes/Drobo2/Data/NEUROPHYS/rig022/161006_E/161006_E_rfori001.bhv';
        case '170724_I_evp003'
            bhvfile = '/Volumes/Drobo2/Data/NEUROPHYS/rig021/170724_I/170724_I_rfori002.bhv';
        otherwise
            bhvfile = [filename '.bhv'];
    end
    
else
    switch BRdatafile
        case {'151205_E_rfori002','151205_E_rfori003','151205_E_cosinteroc001'}
            bhvfile =  '/volumes/LaCie Porsche/151205_E_rfori005.bhv';
        case '161005_E_rfori001'
            bhvfile =  '/volumes/LaCie Porsche/161005_E_rfori002.bhv';
        case '161006_E_evp001'
            bhvfile =  '/volumes/LaCie Porsche/161006_E_rfori001.bhv';
        case '170724_I_evp003'
            bhvfile =  '/volumes/LaCie Porsche/170724_I_rfori002.bhv';
        case '170724_I_mcosinteroc001'
            bhvfile =  '/volumes/LaCie Porsche/170724_I_rfori002.bhv';
        case '180808_I_evp001'
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/180808_I/180808_I_posdisparitydrft004.bhv';
        case {'160609_I_bwflicker001','160609_I_cinterocdrft013','160609_I_colorflicker003'}
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/160609_I/160609_I_cinterocdrft023.bhv';
        case '180810_I_evp001'
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/180810_I/180810_I_conedrft003.bhv';
        case {'180830_I_evp001','180830_I_evp002'}
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/180830_I/180830_I_color003.bhv';
        case '180906_I_evp001'
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/180906_I/180906_I_color001.bhv';
        case '181212_B_evp001'
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/181212_B/181212_B_cinterocdrft003.bhv';
        case '181217_B_evp001'
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/181217_B/181217_B_cinterocdrft002.bhv';
        case '190118_B_evp001'
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/190118_B/190118_B_cinterocdrft001.bhv';
        case '190119_B_evp001'
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/190119_B/190119_B_cinterocdrft002.bhv';
         case {'190120_B_evp001','190120_B_cinterocdrft001'}
        bhvfile = '/Volumes/LaCie Porsche/neurophys_data_for_ana/190120_B/190120_B_color001.bhv';
        case '160423_E_rfori002'
            bhvfile = '/Users/kaciedougherty/Documents/neurophysdata/160423_E_mcosinteroc002.bhv';
        case '190123_B_evp001'
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/190123_B/190123_B_cinterocdrft001.bhv';
        case '190209_B_evp002'
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/190209_B/190209_B_cinterocdrft001.bhv';
        case '190210_B_evp001'
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/190210_B/190210_B_cinterocdrft001.bhv';
        case '190319_B_evp002'
            bhvfile = '/volumes/LaCie Porsche/neurophys_data_for_ana/190319_B/190319_B_cinterocdrft001.bhv';
        otherwise
            bhvfile = [filename '.bhv'];
    end
    
    
end

% check for BHV and NS6
bhv = concatBHV(bhvfile);
ns6file = [filename '.ns6'];
if ~exist(bhvfile,'file') || ~exist(ns6file,'file') || isempty(bhv)
    if strcmp(BRdatafile(end-5:end-3),'evp')
        % grab any other bhv
        bhvlist = dir([brpath '/*.bhv']); bhv = [];
        for l = 1:length(bhvlist)
            try
                bhv = concatBHV([brpath filesep bhvlist(l).name]);
            end
            if ~isempty(bhv)
                break
            end
        end
        if ~isempty(bhv)
            newTP   = [];
            trigger = 'BHV or NS6 missing from disk';
            return
        end
    else
        newTP   = [];
        trigger = 'BHV or NS6 missing from disk';
        return
    end
end

% extract visual stim info from BHV and ypos
Fs        = 30000; % always using NS6
refresh   =  ceil((1/bhv.ActualVideoRefreshRate) * Fs); % time for a screen refresh in sampels
yconv     = 0.5 + (ypos*bhv.PixelsPerDegree/bhv.ScreenYresolution); % accounts for CRT scan path


% determin trigger type if not specified as input
if isempty(trigger)
    if  any(any(strcmp(bhv.TaskObject,'gen(gReferenceDriftingGrating_di)')))
        trigger = 'drifting';
    elseif (strcmp(bhv.TaskObject{end-3},'gen(gGrating_di)') && strcmp(bhv.TaskObject{end-2}([1:3 7:13]),'sqr[0 0 0]') )
        trigger = 'custom';
    elseif strcmp(bhv.PhotoDiodePosition{1},'none')
        trigger = 'none';
        newTP   = [];
        return
    else
        trigger = 'default';
    end
end


% constant re-triggering method
if any(strcmp({'constant'},trigger))
    
    % method specified as "constant"
    % assume constent offset relative to event code
    % good for when no photodiode signal is availiable for triggering
    
    newTP(:,1) = round(TP(:,1) + nr*refresh - refresh*yconv);
    newTP(:,2) = round(TP(:,2) + nr*refresh - refresh*yconv);
    return % stop exicuting here
    
end


% if not using constant mthod, look at BNC signal
ns_header  = openNSx(ns6file,'noread');
ns_labels  = cellfun(@(x) x(1:5),{ns_header.ElectrodesInfo.Label},'UniformOutput',0);
photoidx   = find(strcmp(ns_labels,photolabel));
if isempty(photoidx)
    newTP   = [];
    trigger = sprintf('no channel "%s"',photolabel);
    return
end

channel    = sprintf('c:%u', photoidx);

%  get photodiode signal around event
clear NS
ns_header  = openNSx(ns6file,'noread');
data       = [];
if ~ns_header.RawData.PausedFile
    NS                 = openNSx(ns6file,channel,'read','precision','double','sample');
    data               = NS.Data;
else
    segfiles           = dir([ns6file(1:end-4) '-s' '*']);
    for s = 1:length(segfiles)
        NS   = openNSx([brpath '/' segfiles(s).name],channel,'read','precision','double','sample');
        data = [data NS.Data];
    end
end

switch BRdatafile(8)
    case 'B'
        thresh = 1.7;
    case 'I'
        thresh = 3;
end

if any(strcmp(BRdatafile(1:8),'181207_B'))
    thresh = 1.4; 
end

for tr = 1:size(TP,1)
    clear trstart trend refwin dat detected tf fewer
    
    trstart      = floor(TP(tr,1) - (refresh*5));
    trend        = floor(TP(tr,2) + (refresh*5));
    refwin       = trstart:trend; 
    
    if trend > size(data,2)
           newTP(tr,:)  = [TP(tr,1) nan];
           continue
    end
    
    
    dat          = data(refwin);
    dat          = abs((dat - mean(dat)) / std(dat));
    
    if ~isfield(STIM,'temporal_freq')
        detected  = refwin(find(dat > thresh,1,'first')); % onset of first cycle
        STIM.photo_on{tr} = detected;
        newTP(tr,:)  = [detected nan];
    elseif STIM.temporal_freq(tr) <= 1
        detected  = refwin(find(dat > thresh,1,'first')); % onset of first cycle
        STIM.photo_on{tr} = detected;
        newTP(tr,:)  = [detected nan];
    else
        
        high         = refwin(dat > thresh);
        fewer        = high(diff([refwin(1) high]) > refresh);
        if isempty(fewer)
           thresh = 2;  
           high      = refwin(dat > thresh);
           fewer     = high(diff([refwin(1) high]) > refresh);
        end
        
        if length(fewer) > STIM.temporal_freq(tr)
            newTP(tr,:)       = [fewer(1) nan]; 
            STIM.photo_on{tr} = fewer(1:STIM.temporal_freq(tr));
        elseif isempty(fewer)
            newTP(tr,:)       = [nan nan]; 
            STIM.photo_on{tr} = nan; 
            warning('cant trigger photodiode\n'); 
        else
            newTP(tr,:)       = [fewer(1) nan]; 
            STIM.photo_on{tr} = fewer;
        end
        
        
        
        if any(strcmp(BRdatafile(1:8),'181207_B'))
            if isempty(fewer)
                newTP(tr,:)  = [TP(tr,1) nan];
            else
                newTP(tr,:)  = [fewer(1) nan];
            end
        end
        
%         if  mod(tr,10) == 0 & tr < 500
%         figure, plot(refwin,dat); v = vline(fewer(1)); set(v,'color','k','linewidth',2); 
%         end

    end
end


