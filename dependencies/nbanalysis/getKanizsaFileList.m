% kanizsa file list:
function [fnames] = getKanizsaFileList(sessdate);

switch sessdate
    
    case '140724'
        fnames.kanizsa = {'140724_B_qkanizsa001'};
        fnames.evp     = ['140724_B_evp002'];
        fnames.probe   = 'A';
        
    case '140725'
        fnames.kanizsa = {'140725_B_qkanizsa001'};
        fnames.evp     = ['140725_B_evp005'];
        fnames.probe   = 'A';
    case '140729'
        fnames.kanizsa = {'140729_B_qkanizsa001'};
        fnames.evp     = ['140729_B_evp003'];
        fnames.probe   = 'B';
    case '140731'
        fnames.kanizsa = {'140731_B_qkanizsa001'};
        fnames.evp     = ['140731_B_evp010'];
        fnames.probe   = 'B';
    case '140807'
        fnames.kanizsa = {'140807_B_qkanizsa_center001';'140807_B_qkanizsa_inctr001'};
        fnames.evp     = ['140807_B_evp003'];
        fnames.probe   = 'D';
    case '140818'
        fnames.kanizsa = {'140818_B_qkanizsa_center001';'140818_B_qkanizsa_inctr001'};
        fnames.evp     = ['140818_B_evp004'];
        fnames.probe   = 'D';
    case '140820'
        fnames.kanizsa = {'140820_B_qkanizsa_center001';'140818_B_qkanizsa_inctr001'};
        fnames.evp     = ['140820_B_evp002'];
        fnames.probe   = 'D';
        
        
end