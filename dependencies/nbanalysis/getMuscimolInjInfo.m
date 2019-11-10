function [prefname,postfname] = getMuscimolInjInfo(sdate)

switch sdate
    case '151218'
        prefname = {'005'};
        postfname = {'008'; '009'; '010'; '011'; '012'; '013'; '014'; '015'; '016';...
            '017'; '018'; '019'};
        
    case '160104'
        prefname = {'003'};
        postfname = {'004';'007';'008'; '009'; '010'; '011'; '012'; '013'; '014'; '015'; '016';...
            '017'};
        
    case '160202'
        prefname = {'005'};
        postfname = {'007'};
end

end