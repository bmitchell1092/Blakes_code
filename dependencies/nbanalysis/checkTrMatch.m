function [allpass, message] =  checkTrMatch(grating,NEV)
% grating is the output of the function that reads the textfles written by the gen functions
% NEV is the matching nev file, nust be NEV because NSx timestamp is file stop time

message = cell(3,1);
pass    = ones(3,1);

% extract info from NEV
EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventTimes = floor(NEV.Data.SerialDigitalIO.TimeStampSec .* 1000); %ms, to match 1kHz
[pEvC, ~]  = parsEventCodesML(EventCodes,EventTimes);

% TEST 1 -- check that trial numbers as recorded in the text file are sequnetial
if any(unique(abs(diff(grating.trial))) > 1 )
    message{1} = 'FAIL: trial numbers as recorded in gTextFile are not sequential';
    pass(1)    = 0;
% elseif sum(grating.trial==1) ~= sum(grating.trial==grating.trial(end))
%     message{1} = 'FAIL: trial numbers as recorded in gTextFile may not be sequential';
%     pass(1)    = 0;
else
    message{1} = 'PASS: gTextFile trial numbers sequential';
end

% TEST 2 -- check that trial numbers as recorded in the event codes are sequnetial
if mode(cellfun(@(x) sum(x >= 116),pEvC)) == 3
    % this file sent block, condition number, and last diget of
    % trial number after the inital 999 trial start
    
    trlnums = double(cellfun(@(x) x(6) - 116, pEvC));
    if ~isequal(unique(diff(trlnums)),[-9 1]) 
        if length(trlnums) <= 9 && isequal(unique(diff(trlnums)),[1])
            message{2} = 'PASS: Event Codes trial numbers sequential';
        else
            message{2} = 'FAIL: trial numbers as recorded in Event Codes are not sequential';
            pass(2)    = 0;
        end
    else
        message{2} = 'PASS: Event Codes trial numbers sequential';
    end
else
    message{2} = 'NULL: could not test';
    pass(2)    = NaN;
end

% TEST 3 --  compare timestamps between ML and BR files
fields = fieldnames(grating);
if any(strcmp(fields,'timestamp'))
    mlstart = grating.timestamp(1);
    brstart = datenum(NEV.MetaTags.DateTime,'dd/mmm/yyyy dddd HH:MM:SS');
    if mlstart < brstart
        message{3} = 'FAIL: 1st TimeStamp in ML file is before BR record time';
        pass(3)    = 0;
    else
        message{3} = 'PASS: BR record time precedes ML timstamps';
    end
else
    message{3} = 'NULL: could not test';
    pass(3)    = NaN;
end

pass(isnan(pass)) = [];
if all(pass)
    allpass = true;
else
    allpass = false;
end
