function [STIM,spkTPs,all_cyc,realtr, veconsets] = sortTrialData(BRdatafile,brdrname,pEvC,pEvT,grating,stimfeatures,flag_RewardedTrialsOnly,badobs,extension,use_evcodes)

% BRdatafile   : filename
% brdrname     : where ns6/ns2 files are located
% pEvC, pEvT   : cell arrays with event codes and times for every trial
% grating      : variable from text file listing stimulus parameters every tr
% stimfeatures : list of relevant features for analysis to be saved in STIM
% flag_RewardedTrialsOnly : include only trs where animal received juice
% badobs       : vector , exclude these trials
% extension    : which file type (ns2 or ns6) to use for photodiode signal
% use_evcodes  : if 1, don't use photodiode signal, just use eventcodes


% spkTPs       : onset of stimulus for each trial (ex. badobs, )
% all_cyc      : onset of each cycle for each trial if drifting grating
%

% sessions before August 2016--> trigger to event codes
% can only look at spike rate etc. for these data, no time locked figures

n = datenum(BRdatafile(1:6),'yymmdd');
all_cyc = [];
veconsets = [];

if ~isempty(strfind(BRdatafile,'drft')) || ~isempty(strfind(BRdatafile,'flicker'))
    
    if n >= datenum('08/03/2016','mm/dd/yyyy') && use_evcodes ~= 1
        
        filename = fullfile(brdrname,BRdatafile);
        if strcmp(extension,'ns2')
            sampledown = 1;
        else
            sampledown = 0; % if you set this to 1, the stimulus onset times will be at a sampling rate of 1000
        end
        [BNC, BNCLabels,fs] = getBNCData({'ainp1'},filename,sampledown);
        
        realtr = [];
        obs = 0; clear spkTPs STIM
        for t = 1: length(pEvC)
            
            stim =  find(grating.trial == t); if any(diff(stim) ~= 1); error('check grating file'); end
            codes = pEvC{t};
            onsets = find(codes == 23);
            
            for p = 1:length(onsets)
                
                obs = obs +1;
                realtr = [realtr t];
                
                if t == badobs || flag_RewardedTrialsOnly && ~ismember(96,codes)
                    
                    init_onset(obs) = 0;
                    for f = 1:length(stimfeatures)
                        STIM.(stimfeatures{f})(obs,:) = grating.(stimfeatures{f})(stim(p));
                    end
                    all_cyc{obs} = [];
                    spkTPs(obs) = 0;
                    
                else
                    
                    refreshdelay = fs./85 +200;
                    
                    trstart      = floor((pEvT{t}(1)+ (500.*(fs/1000)))); % only consider samples > 50 ms in (there's an initial flash)
                    trend        = floor(pEvT{t}(end));
                    refwin  = trstart:trend;
                    sig     = BNC(refwin,1);
                    thr = mean(sig) -  2*std(sig,0,1);%mean(sig) - 1.5*std(sig,0,1);
                    spkTPs(obs) = refwin(find(sig<thr,1,'first')); % onset of first cycle
                    
                    if spkTPs(obs) == 0
                        fprintf('%u %u\n',t,spkTPs(obs))
                    end
                    if grating.temporal_freq(t) <= 1
                        subpts = spkTPs(obs);
                    else
                        allp =  find(sig<thr);
                        pts    = refwin(allp(diff(allp)>1));
                        subpts = [pts(1) pts((find(diff(pts)>(refreshdelay)))+1)];
                    end
                    all_cyc{obs} = double(subpts); % onset of all cycles in trial
                    
                    % %uncomment to plot a examples of triggering with photodiode signal
                    %                                 if mod(obs,4) == 0
                    %                                   figure, plot(refwin,BNC(refwin)); hold on; v = vline(subpts);
                    %                                   set(v,'LineWidth',2);
                    %
                    %                                 end
                    %
                    if length(grating.trial(stim)) == 2 || p == 1
                        for f = 1:length(stimfeatures)
                            STIM.(stimfeatures{f})(obs,:) = grating.(stimfeatures{f})(stim(p));
                        end
                    elseif length(grating.trial(stim)) > 2 && p > 1
                        for f = 1:length(stimfeatures)
                            STIM.(stimfeatures{f})(obs,:) = grating.(stimfeatures{f})(stim(p*2-1));
                        end
                    end
                    
                end
            end
        end
        
        
    else
        
        realtr = [];
        obs = 0; clear spkTPs STIM
        for t = 1: length(pEvC)
            
            stimon  =  pEvC{t} == 23 | pEvC{t} == 25 | pEvC{t} == 27 | pEvC{t} == 29 | pEvC{t} == 31;
            stimoff =  pEvC{t} == 24 | pEvC{t} == 26 | pEvC{t} == 28 | pEvC{t} == 30 | pEvC{t} == 32;
            
            st = pEvT{t}(stimon);
            en = pEvT{t}(stimoff);
            
            stim =  find(grating.trial == t); if any(diff(stim) ~= 1); error('check grating file'); end
            
            for p = 1:length(en)
                obs = obs + 1;
                realtr = [realtr t];
                if  any(t == badobs) || ...
                        (flag_RewardedTrialsOnly && ~any(pEvC{t} == 96))
                    
                    spkTPs(obs,:) = [0 0];
                    
                    for f = 1:length(stimfeatures)
                        STIM.(stimfeatures{f})(obs,:) = NaN;
                    end
                    
                else
                    
                    spkTPs(obs,:) = [st(p) en(p)];
                    
                    for f = 1:length(stimfeatures)
                        STIM.(stimfeatures{f})(obs,:) = grating.(stimfeatures{f})(stim(p));
                    end
                    
                end
            end
        end
        
    end
    
    
else
    
    % STATIC GRATINGS:
    
    if n >= datenum('08/03/2016','mm/dd/yyyy') && use_evcodes ~= 1
        
        filename = fullfile(brdrname,BRdatafile);
        if strcmp(extension,'ns2')
            sampledown = 1;
        else
            sampledown = 0; % if you set this to 1, the stimulus onset times will be at a sampling rate of 1000
        end
        [BNC, BNCLabels,fs] = getBNCData({'ainp1'},filename,sampledown);
        
        realtr = [];
        obs = 0; clear spkTPs STIM
        for f = 1:length(stimfeatures)
            STIM.(stimfeatures{f}) = [];
        end
        veconsets = [];
        for t = 1: length(pEvC)
            
            stim =  find(grating.trial == t); if any(diff(stim) ~= 1); error('check grating file'); end
            codes = pEvC{t};
            onsets = find( codes == 23 | codes == 25 | codes == 27 | codes == 29| codes == 31 );
            
            
            realtr = [realtr t];
            
            
            if t~= badobs && ismember(96,codes)
                %
                
                refreshdelay = fs./85 + 200;
                
                trstart      = floor((pEvT{t}(1)+ (50.*(fs/1000)))); % only consider samples > 50 ms in (there's an initial flash)
                trend        = floor(pEvT{t}(end));
                refwin  = trstart:trend;
                sig     = BNC(refwin,1);
                thr = mean(sig) -  1.5*std(sig,0,1);%mean(sig) - 1.5*std(sig,0,1);
                
                
                allp   = find(sig<thr);
                pts    = [refwin(allp(diff(allp)>1)+1) refwin(allp(end))];
                subpts = [pts(1) pts((find(diff(pts)>(refreshdelay)))+1)];
                
                if length(subpts) == 10
                    obs = obs +1;
                    spkTPs(obs) = refwin(find(sig<thr,1,'first')); % onset of first cycle
                    
                    all_cyc{obs} = double(subpts); % onset of all cycles in trial
                    
                    % % uncomment to plot a examples of triggering with photodiode signal
                    %                                 if mod(obs,8) == 0
                    %                                   figure, plot(refwin,BNC(refwin)); hold on;
                    %                                   v = vline(subpts(1:2:end));set(v,'Color','r');
                    %                                   set(v,'LineWidth',3,'LineStyle','-');
                    %                                      v = vline(subpts(2:2:end));set(v,'Color','b');
                    %                                   set(v,'LineWidth',3,'LineStyle','-');
                    %
                    %                                  end
                    
                    
                    for p = 1:length(onsets)
                        for f = 1:length(stimfeatures)
                            STIM.(stimfeatures{f}) =[STIM.(stimfeatures{f}); grating.(stimfeatures{f})(stim(p))'];
                        end
                    end
                end
                
            end
        end
        
        
        
    else
        
        realtr = [];
        obs = 0; clear spkTPs STIM
        for t = 1: length(pEvC)
            
            stimon  =  pEvC{t} == 23 | pEvC{t} == 25 | pEvC{t} == 27 | pEvC{t} == 29 | pEvC{t} == 31;
            stimoff =  pEvC{t} == 24 | pEvC{t} == 26 | pEvC{t} == 28 | pEvC{t} == 30 | pEvC{t} == 32;
            
            st = pEvT{t}(stimon);
            en = pEvT{t}(stimoff);
            
            stim =  find(grating.trial == t); if any(diff(stim) ~= 1); error('check grating file'); end
            
            for p = 1:length(en)
                obs = obs + 1;
                realtr = [realtr t];
                if  any(t == badobs) || ...
                        (flag_RewardedTrialsOnly && ~any(pEvC{t} == 96))
                    
                    spkTPs(obs,:) = [0 0];
                    
                    for f = 1:length(stimfeatures)
                        STIM.(stimfeatures{f})(obs,:) = NaN;
                    end
                    
                else
                    
                    spkTPs(obs,:) = [st(p) en(p)];
                    
                    for f = 1:length(stimfeatures)
                        STIM.(stimfeatures{f})(obs,:) = grating.(stimfeatures{f})(stim(p));
                    end
                    
                end
            end
        end
        
    end
    
end