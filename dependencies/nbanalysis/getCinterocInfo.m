% session dates and good channels for cinteroc :

function [theseelectrodes, BRdatafiles, units,alldates,aligncodes,flipeye] = getCinterocInfo(date,ftype)
aligncodes = [];
flipeye = [];

cpatch_alldates = {'160406_I_cpatch013';'160407_I_cpatch005';'160412_I_cpatch006';'160412_I_cpatch012';...
    '160411_I_cpatch002';'160411_I_cpatch010';'160414_I_cpatch002';'160414_I_cpatch008';'160417_I_cpatch003';'160419_I_cpatch003';...
'160421_I_cpatch003';'160421_I_cpatch011';'160421_I_cpatch021'};


cinteroc_alldates = {'160427_I_cinterocdrft006';...
    '160429_I_cinterocdrft003';...
    '160502_I_cinterocdrft003';...
    '160502_I_cinterocdrft004';...
    '160502_I_cinterocdrft008';...
    '160503_I_cinterocdrft002';...
    '160505_I_cinterocdrft008';...
    '160505_I_cinterocdrft010';...
    '160506_I_cinterocdrft002';...
    '160506_I_cinterocdrft005';...
    '160509_I_cinterocdrft002';...
    '160509_I_cinterocdrft004';...
    '160510_I_cinterocdrft001';...
    '160510_I_cinterocdrft005';...
    '160512_I_cinterocdrft004';...
    '160512_I_cinterocdrft005';...
    '160513_I_cinterocdrft003';...
    '160513_I_cinterocdrft004';...
    '160513_I_cinterocdrft006';...
    '160602_I_cinterocdrft012';...
    '160609_I_cinterocdrft013';...
    '160609_I_cinterocdrft023';...
    '160611_I_cinterocdrft001';...
    '160611_I_cinterocdrft003';...
    '160613_I_cinterocdrft005';...
    '160613_I_cinterocdrft007';...
    '160615_I_cinterocdrft001';...
    '160615_I_cinterocdrft005';...
    '160616_I_cinterocdrft001';...
    '160616_I_cinterocdrft002';...
    '160616_I_cinterocdrft004';...
    '160616_I_cinterocdrft005';...
    '160616_I_cinterocdrft006';...
    '160623_I_cinterocdrft001';...
    '160623_I_cinterocdrft002';...
    '160623_I_cinterocdrft003';...
    '160623_I_cinterocdrft005';...
    '160623_I_cinterocdrft006';...
    '160624_I_cinterocdrft001';...
    '160624_I_cinterocdrft002';...
    '160624_I_cinterocdrft004';...
    '160625_I_cinterocdrft001';...
    '160627_I_cinterocdrft001';...
    '160627_I_cinterocdrft004';...
    '160627_I_cinterocdrft007';...
    '160628_I_cinterocdrft001';...
    '160628_I_cinterocdrft004';...
    '160628_I_cinterocdrft005';...
    '160628_I_cinterocdrft007';...
    '160628_I_cinterocdrft008';...
    '160628_I_cinterocdrft009'};
    
if strcmp(ftype,'cpatch')
    alldates = cpatch_alldates; 
else
    alldates = cinteroc_alldates; 
end
    BRdatafiles = {}; units = [];aligncodes = [];
    switch date
        case '151214'
            theseelectrodes = {'eC19'};
            
            
        case '151229'
            theseelectrodes = {'eC15'};
            BRdatafiles = {'151229_I_cinteroc007';'151229_I_cinteroc008'} ;
            
        case '160202'
            BRdatafiles = {'160202_I_cinteroc002'};
            
        case '160205'
            theseelectrodes = {'White'};
            BRdatafiles = {'160205_I_cinteroc001';'160205_I_cinteroc003';'160205_I_cinteroc004'};
        case '160331_I_cinteroc004'
            theseelectrodes = {'eD22'};
            BRdatafiles = {'160331_I_cinteroc004'};
            
        case '160406_I_cpatch013'
            theseelectrodes = {'eD18'; 'eD19'};
            units(1).num = 1;     % BR numbered units for channel eD18
            units(2).num = [1]; % BR numbered units for channel eD19
            BRdatafiles = {'160406_I_cpatch013'};
            
        case '160407_I_cpatch005'
            theseelectrodes = {'eD02';'eD06';'eD09'};
            flipeye  = [0 0 1];
            units(1).num = 1;     % BR numbered units for channel eD18
            units(2).num = 2;     % BR numbered units for channel eD18
            units(3).num = [1 2 3];     % BR numbered units for channel eD18
            
            BRdatafiles = {'160407_I_cpatch005'};
            
        case '160412_I_cpatch006'
            theseelectrodes = {'eD16';'eD15';'eD14'};
            units(1).num = [2];     % BR numbered units for channel eD18
            units(2).num = [1 2];     % BR numbered units for channel eD18
            units(3).num = [1 2];     % BR numbered units for channel eD18
            BRdatafiles = {'160412_I_cpatch006'};
            
        case '160412_I_cpatch012'
            theseelectrodes = {'eD16';'eD15';'eD14'};
            units(1).num = [1];     % BR numbered units for channel eD18
            units(2).num = [2];     % BR numbered units for channel eD18
            units(3).num = [1 2];     % BR numbered units for channel eD18
            BRdatafiles = {'160412_I_cpatch012'};
            
        case '160411_I_cpatch002'
            theseelectrodes = {'eD16'; 'eD18'};
            flipeye   = [0 1];
            units(1).num = [2];     % BR numbered units for channel eD18
            units(2).num = [1 2];     % BR numbered units for channel eD18
            BRdatafiles = {'160411_I_cpatch002'};
            
        case '160411_I_cpatch010'
            theseelectrodes = {'eD19'};
            flipeye  = [1]; 
            units(1).num = [1];     % BR numbered units for channel eD18
            
            BRdatafiles = {'160411_I_cpatch010'};
            
        case '160414_I_cpatch002'
            theseelectrodes = {'eD13'; 'eD12'; 'eD18'; 'eD19'};
            units(1).num = [1 2];
            units(2).num = [1 2];
            units(3).num = [1 2];
            units(4).num = [2];
            BRdatafiles = {'160414_I_cpatch002'};
            
        case '160414_I_cpatch008'
            theseelectrodes = {'eD15'; 'eD16'; 'eD18'};
            units(1).num = [1];
            units(2).num = [1];
            units(3).num = [1 2];
            
            BRdatafiles = {'160414_I_cpatch008'};
            
        case '160417_I_cpatch003'
            theseelectrodes = {'eD12'};
            units(1).num = [1 2 3];
            BRdatafiles = {'160417_I_cpatch003'};
            
        case '160419_I_cpatch003'
            theseelectrodes = {'eD18'};
            units(1).num = [1 2];
            BRdatafiles = {'160419_I_cpatch003'};
        case '160421_I_cpatch003'
            theseelectrodes = {'eD19'};
            units(1).num = [1];
            BRdatafiles = {'160421_I_cpatch003'};
        case '160421_I_cpatch011'
            theseelectrodes = {'eD12'};
            units(1).num = [];
            BRdatafiles = {'160421_I_cpatch011'};
        case '160421_I_cpatch021'
            theseelectrodes = {'eD22'; 'eD16'};
            units(1).num = [];
        BRdatafiles = {'160421_I_cpatch021'};
        %%%%%%%%%%%%%%%%%%  CINTEROCDRFT  v%%%%%%%%%%%%%%%%%%%%%
        %     case '160421_I_cinterocdrft006'
        %         theseelectrodes = {'eD11'; 'eD12'};
        %         units(1).num     = [1];
        %         units(1).rank    = [1];
        %         units(1).flipeye = [0];
        %         units(2).num = [1 2];
        %         units(2).flipeye = [0 0];
        %         units(2).rank = [4 4];
        %         BRdatafiles = {'160419_I_cinterocdrft006'};
        %
    case '160427_I_cinterocdrft006'
        theseelectrodes  = {'eD13'; 'eD17'};
        aligncodes       = 0;
        units(1).num     = [1 2];
        units(1).rank    = [1 1];
        units(1).flipeye = [0 0];
        BRdatafiles = {'160427_I_cinterocdrft006'};
        
    case '160429_I_cinterocdrft003'
        theseelectrodes  = {'eD12'; 'eD15'; 'eD18'; 'eD20'};
        flipeye          = [1 0 1 1];
        aligncodes       = 0;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [1];
        
        units(2).num     = [1];
        units(2).rank    = [1];
        units(2).flipeye = [0];
        
        units(3).num     = [2];
        units(3).rank    = [1];
        units(3).flipeye = [0];
        
        units(4).num     = [1];
        units(4).rank    = [1];
        units(4).flipeye = [0];
        
        BRdatafiles = {'160429_I_cinterocdrft003'};
        
    case '160502_I_cinterocdrft003'
        theseelectrodes  = {'eD12'; 'eD17'};
        flipeye = [1 1];
        aligncodes       = 0;
        units(1).num     = [1 2];
        units(1).rank    = [2 2];
        units(1).flipeye = [0 0];
        
        units(2).num     = [1 2];
        units(2).rank    = [2 2];
        units(2).flipeye = [0 0];
        
        BRdatafiles = {'160502_I_cinterocdrft003'};
        
        
    case '160502_I_cinterocdrft004'
        theseelectrodes  = {'eD03';'eD17'; 'eD18'};
        flipeye          = [0 1 1];
        aligncodes       = 0;
        units(1).num     = [1 2];
        units(1).rank    = [2 2];
        units(1).flipeye = [0 0];
        
        units(2).num     = [1];
        units(2).rank    = [2];
        units(2).flipeye = [1];
        
        units(1).num     = [1 2];
        units(1).rank    = [2 2];
        units(1).flipeye = [0 0];
        BRdatafiles = {'160502_I_cinterocdrft004'};
        
    case '160502_I_cinterocdrft008'
        theseelectrodes  = {'eD15';'eD16'; 'eD17'};
        flipeye          = [1 1 1]; 
        aligncodes       = 1;
        units(1).num     = [1];
        units(1).rank    = [2];
        units(1).flipeye = [0];
        
        units(2).num     = [1];
        units(2).rank    = [1];
        units(2).flipeye = [0];
        
        BRdatafiles = {'160502_I_cinterocdrft008'};
        
    case '160503_I_cinterocdrft002'
        theseelectrodes  = {'eD13';'eD14'; 'eD15'; 'eD17';'eD21'};
        flipeye          = [0 0 0 1 0]; 
        aligncodes       = 1;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [0];
        
        units(2).num     = [1];
        units(2).rank    = [1];
        units(2).flipeye = [1];
        
        units(3).num     = [1];
        units(3).rank    = [1];
        units(3).flipeye = [1];
        
        units(4).num     = [1];
        units(4).rank    = [1];
        units(4).flipeye = [1];
        BRdatafiles = {'160503_I_cinterocdrft002'};
        
    case '160505_I_cinterocdrft002'
        theseelectrodes  = {'eD15'; 'eD16'};
        
    case '160505_I_cinterocdrft008'
        theseelectrodes  = {'eD21'; 'eD24'};
        flipeye          = [0 1]; 
        aligncodes       = 0;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [1];
        
        units(2).num     = [1];
        units(2).rank    = [1];
        units(2).flipeye = [0];
        
        BRdatafiles = {'160505_I_cinterocdrft008'};
        
    case '160505_I_cinterocdrft010'
        theseelectrodes  = {'eD16'};
        aligncodes       = 0;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [1];
        
        
        BRdatafiles = {'160505_I_cinterocdrft010'};
        
    case '160506_I_cinterocdrft002'
        theseelectrodes  = {'eD18';'eD20'};
        flipeye          = [0 1];
        aligncodes       = 0;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [0];
        
        units(2).num     = [1];
        units(2).rank    = [1];
        units(2).flipeye = [0];
        
        BRdatafiles = {'160506_I_cinterocdrft002'};
        
    case '160506_I_cinterocdrft005'
        theseelectrodes  = {'eD18'};
        flipeye          = [1];
        aligncodes       = 0;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [1];
        
        
        BRdatafiles = {'160506_I_cinterocdrft005'};
        
        
    case '160509_I_cinterocdrft002'
        theseelectrodes  = {'eD09'; 'eD11'; 'eD08'};
        flipeye          = [1 0 1]; 
        aligncodes       = 1;
        units(1).num     = [1 2];
        units(1).rank    = [1 1];
        units(1).flipeye = [0 0];
        
        units(2).num     = [2];
        units(2).rank    = [1];
        units(2).flipeye = [0];
        
        units(3).num     = [1];
        units(3).rank    = [1];
        units(3).flipeye = [1];
        
        BRdatafiles = {'160509_I_cinterocdrft002'};
        
    case '160509_I_cinterocdrft004'
        theseelectrodes  = {'eD08';'eD10'};
        flipeye          = [1 0]; 
        aligncodes       = 0;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [0];
        
        units(2).num     = [1];
        units(2).rank    = [1];
        units(2).flipeye = [1];
        
        BRdatafiles = {'160509_I_cinterocdrft004'};
        
    case '160510_I_cinterocdrft001'
        theseelectrodes  = {'eD13';'eD12'};
        aligncodes       = 1;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [1];
        
        units(2).num     = [1];
        units(2).rank    = [1];
        units(2).flipeye = [1];
        
        BRdatafiles = {'160510_I_cinterocdrft001'};
        
        
    case '160510_I_cinterocdrft005'
        theseelectrodes  = {'eD07'; 'eD08'; 'eD11'; 'eD12'; 'eD13'};
        flipeye          = [1 1 0 0 0]; 
        aligncodes       = 1;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [1];
        
        units(2).num     = [1];
        units(2).rank    = [1];
        units(2).flipeye = [0];
        
        units(3).num     = [1];
        units(3).rank    = [1];
        units(3).flipeye = [1];
        
        units(4).num     = [1 2];
        units(4).rank    = [1 1];
        units(4).flipeye = [1 1];
        
        units(5).num     = [1];
        units(5).rank    = [1];
        units(5).flipeye = [1];
        BRdatafiles = {'160510_I_cinterocdrft005'};
        
    case '160512_I_cinterocdrft004'
        theseelectrodes  = {'eD14'};
        aligncodes       = 1;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [1];
        
        
        BRdatafiles = {'160512_I_cinterocdrft004'};
        
    case '160512_I_cinterocdrft005'
        theseelectrodes  = {'eD13';'eD14'};
        aligncodes       = 1;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [1];
        
        units(2).num     = [1];
        units(2).rank    = [1];
        units(2).flipeye = [0];
        
        BRdatafiles = {'160512_I_cinterocdrft005'};
        
    case '160513_I_cinterocdrft003'
        theseelectrodes  = {'eD16'};
        aligncodes       = 1;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [1];
        
        BRdatafiles = {'160513_I_cinterocdrft003'};
        
    case '160513_I_cinterocdrft004'
        theseelectrodes  = {'eD15';'eD16'};
        aligncodes       = 1;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [1];
        
        units(2).num     = [1];
        units(2).rank    = [1];
        units(2).flipeye = [1];
        BRdatafiles = {'160513_I_cinterocdrft004'};
        
    case '160513_I_cinterocdrft006'
        theseelectrodes  = {'eD15'};
        aligncodes       = 1;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [1];
        BRdatafiles = {'160513_I_cinterocdrft006'};
        
    case '160602_I_cinterocdrft012'
        theseelectrodes  = {'eD14'};
        aligncodes       = 1;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [0];
        BRdatafiles = {'160602_I_cinterocdrft012'};
        
    case '160609_I_cinterocdrft013'
        theseelectrodes  = {'eD11'};
        aligncodes       = 1;
        units(1).num     = [1];
        units(1).rank    = [1];
        units(1).flipeye = [0];
        BRdatafiles = {'160609_I_cinterocdrft013'};
        
    case '160609_I_cinterocdrft023';
        theseelectrodes  = {'eD17'};
        aligncodes       = 1;
        units(1).num     = [2];
        units(1).rank    = [1];
        units(1).flipeye = [0];
        BRdatafiles = {'160609_I_cinterocdrft023'};
        
    case '160611_I_cinterocdrft001';
         aligncodes       = 1;
        theseelectrodes  = {'eD14';'eD15'};
        BRdatafiles = {'160611_I_cinterocdrft001'};
        
        
    case '160611_I_cinterocdrft003';
         aligncodes       = 1;
        theseelectrodes  = {'eD10';'eD14';'eD15'};
        BRdatafiles = {'160611_I_cinterocdrft003'};
        
    case '160613_I_cinterocdrft003';
         aligncodes       = 1;
        theseelectrodes  = {'eD04';'eD07'};
        BRdatafiles = {'160613_I_cinterocdrft003'};
        
    case '160613_I_cinterocdrft005';
         aligncodes       = 1;
        theseelectrodes  = {'eD10';'eD14';'eD15'};
        BRdatafiles = {'160613_I_cinterocdrft005'};
        
    case '160613_I_cinterocdrft007';
         aligncodes       = 1;
        theseelectrodes  = {'eD09'};
        BRdatafiles = {'160613_I_cinterocdrft007'};
        
        
    case '160615_I_cinterocdrft001';
         aligncodes       = 1;
        theseelectrodes  = {'eD16'};
        BRdatafiles = {'160615_I_cinterocdrft001'};
        
    case '160615_I_cinterocdrft005';
         aligncodes       = 1;
        theseelectrodes  = {'eD21';'eD22'};
        BRdatafiles = {'160615_I_cinterocdrft005'};
        
    case '160616_I_cinterocdrft001';
         aligncodes       = 1;
        theseelectrodes  = {'eD15';'eD18';'eD20';'eD21';'eD24'};
        BRdatafiles = {'160616_I_cinterocdrft001'};
        
        
    case '160616_I_cinterocdrft002';
         aligncodes       = 1;
        theseelectrodes  = {'eD22'};
        BRdatafiles = {'160616_I_cinterocdrft002'};
        
    case '160616_I_cinterocdrft004';
         aligncodes       = 1;
        theseelectrodes  = {'eD10'};
        BRdatafiles = {'160616_I_cinterocdrft004'};
        
    case '160616_I_cinterocdrft005';
         aligncodes       = 1;
        theseelectrodes  = {'eD10'};
        BRdatafiles = {'160616_I_cinterocdrft005'};
        
    case '160616_I_cinterocdrft006';
         aligncodes       = 1;
        theseelectrodes  = {'eD15'};
        BRdatafiles = {'160616_I_cinterocdrft006'};
        
        
    case '160623_I_cinterocdrft001';
         aligncodes       = 1;
        theseelectrodes  = {'eD15'};
        BRdatafiles = {'160623_I_cinterocdrft001'};
        
    case '160623_I_cinterocdrft002';
         aligncodes       = 1;
        theseelectrodes  = {'eD18';'eD20'};
        BRdatafiles = {'160623_I_cinterocdrft002'};
        
    case '160623_I_cinterocdrft003';
         aligncodes       = 1;
        theseelectrodes  = {'eD20'};
        BRdatafiles = {'160623_I_cinterocdrft003'};
        
    case '160623_I_cinterocdrft005';
         aligncodes       = 1;
        theseelectrodes  = {'eD07'};
        BRdatafiles = {'160623_I_cinterocdrft005'};
        
    case '160623_I_cinterocdrft006';
         aligncodes       = 1;
        theseelectrodes  = {'eD23';'eD24'};
        BRdatafiles = {'160623_I_cinterocdrft006'};
        
    case '160624_I_cinterocdrft001';
         aligncodes       = 1;
        theseelectrodes  = {'eD23'};
        BRdatafiles = {'160624_I_cinterocdrft001'};
        
    case '160624_I_cinterocdrft002';
         aligncodes       = 1;
        theseelectrodes  = {'eD04'};
        BRdatafiles = {'160624_I_cinterocdrft002'};
        
    case '160624_I_cinterocdrft004';
         aligncodes       = 1;
        theseelectrodes  = {'eD24'};
        BRdatafiles = {'160624_I_cinterocdrft004'};
        
    case '160625_I_cinterocdrft001';
         aligncodes       = 1;
        theseelectrodes  = {'eD09'};
        BRdatafiles = {'160625_I_cinterocdrft001'};
        
    case '160627_I_cinterocdrft001';
         aligncodes       = 1;
        theseelectrodes  = {'eD22'};
        BRdatafiles = {'160627_I_cinterocdrft001'};
        
    case '160627_I_cinterocdrft004';
         aligncodes       = 1;
        theseelectrodes  = {'eD22';'eD23'};
        BRdatafiles = {'160627_I_cinterocdrft004'};
        
    case '160627_I_cinterocdrft007';
         aligncodes       = 1;
        theseelectrodes  = {'eD21';'eD23'};
        BRdatafiles = {'160627_I_cinterocdrft007'};
        
    case '160628_I_cinterocdrft001';
         aligncodes       = 1;
        theseelectrodes  = {'eD21';'eD22';'eD23'};
        BRdatafiles = {'160628_I_cinterocdrft001'};
        
    case '160628_I_cinterocdrft004';
         aligncodes       = 1;
        theseelectrodes  = {'eD21';'eD22'};
        BRdatafiles = {'160628_I_cinterocdrft004'};
        
    case '160628_I_cinterocdrft005';
         aligncodes       = 1;
        theseelectrodes  = {'eD20';'eD21'};
        BRdatafiles = {'160628_I_cinterocdrft005'};
        
    case '160628_I_cinterocdrft007';
         aligncodes       = 1;
        theseelectrodes  = {'eD11';'eD13'};
        BRdatafiles = {'160628_I_cinterocdrft007'};
        
    case '160628_I_cinterocdrft008';
         aligncodes       = 1;
        theseelectrodes  = {'eD24'};
        BRdatafiles = {'160628_I_cinterocdrft008'};
        
    case '160628_I_cinterocdrft009';
         aligncodes       = 1;
        theseelectrodes  = {'eD21'};
        BRdatafiles = {'160628_I_cinterocdrft009'};
end
