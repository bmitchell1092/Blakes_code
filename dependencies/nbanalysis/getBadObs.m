function badobs = getBadObs(BRdatafile)

switch BRdatafile
    
    case '151125_E_dotmapping001'
        badobs = [60 164 87 91];
    case '151203_E_rfori003'
        badobs = [156];
    case '151204_E_dotmapping002'
        badobs = [18 22 85 170 185 188 217];
    case '151204_E_dotmapping004'
        badobs = [198 169 137 20 16 24];
    case '151206_E_dotmapping001'
        badobs = 62;
    case '151208_E_kanizsa001'
        badobs = 369; 
    otherwise
        badobs = [];
        
end



