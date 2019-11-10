clear, close all force
warning('off','all')
addpath(genpath('/Users/kaciedougherty/Documents/GitHub/KiloSortUtils/offlineBRAutoSort/')); 
addpath(genpath('/users/kaciedougherty/documents/code/nbanalysis/')); 
addpath(genpath('/users/kaciedougherty/documents/code/lgn-dichoptic'));
addpath(genpath('/users/kaciedougherty/documents/code/nbanalysis/')); 
addpath(genpath('/users/kaciedougherty/documents/code/fNPMK')); 

bhvfilename   = '/users/kaciedougherty/documents/neurophysdata/190709_B/190709_B_calibration001.bhv'; 
nsfilename    = '/users/kaciedougherty/documents/neurophysdata/190709_B/190709_B_calibration001.ns2'; 

colors                = getSerialColors; 

[eye_dva,bhv,eye_ana] = getAnalogEyeinDVA(nsfilename,bhvfilename); 
[brdrname,BRdatafile] = fileparts(nsfilename); 
%[STIM]                = getEventTimeInfo([brdrname '/'],BRdatafile); 


filename        = [brdrname '/' BRdatafile];
NEV             = openNEV(strcat(filename,'.nev'),'noread','overwrite');
EventCodes      = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventSampels    = NEV.Data.SerialDigitalIO.TimeStamp;
[pEvC, pEvT]    = parsEventCodesML(EventCodes,EventSampels);


halfscreen       = (bhv.ScreenXresolution./2)./bhv.PixelsPerDegree; 

%%

for tr = 3%:10 %length(pEvC)
    
    if ~any(ismember(pEvC{tr},96)), continue,end
    
    clear fix_t refwin 
   
    fix_t  = floor(pEvT{tr}(ismember(pEvC{tr},35))./30); 
    
    refwin = fix_t - 96: fix_t + 400; 
    
    figure, set(gcf,'color','w','position',[50 50 800 500]); 
    subplot(1,2,1)
    yyaxis left
    h(1) = plot(eye_ana(refwin,1),'linewidth',2); 
    ylabel('10^-4 V')
    
    yyaxis right
    bhv_ana      = [bhv. AnalogData{tr}.General.Gen1  bhv. AnalogData{tr}.General.Gen2]; 
    h(2) = plot(bhv_ana(:,1),'linewidth',2); 
    legend(h,{'BR ana','BHV ana'})
    set(gca,'tickdir','out','box','off'); 
    ylabel('unknown units')
    
    subplot(1,2,2)
    yyaxis left
    bhv_dva   = bhv.AnalogData{tr}.EyeSignal; 
    h(1)      = plot(bhv_dva(:,1),'linewidth',2); 
    set(gca,'tickdir','out','box','off'); 
    ylabel('dva')
    
    yyaxis right
    h(2) = plot(eye_dva(refwin,1),'linewidth',2); 
    legend(h,{'BHV dva','BR dva'})
    ylabel('dva')

end

