function [unit,filelist] = getUnitInfo(unitid)

filelist = {'160609_I_cinterocdrft013_11';...
    '160609_I_cinterocdrft016_07';...
    '160609_I_cinterocdrft023_17';...
    '160611_I_cinterocdrft003_10';...
    
    '160613_I_cinterocdrft005_02';...
    '160613_I_cinterocdrft005_14';...
    '160613_I_cinterocdrft005_15';...
    '160613_I_cinterocdrft008_09';...
    '160613_I_cinterocdrft009_06';...
    
    '160616_I_cinterocdrft002_22';...
    '160616_I_cinterocdrft004_10';...
    '160616_I_cinterocdrft005_10';...
    '160616_I_cinterocdrft006_15';...
    
    '160623_I_cinterocdrft001_23';...
    '160623_I_cinterocdrft002_18';...
    '160623_I_cinterocdrft003_20';...
    '160623_I_cinterocdrft005_07';...
    
    '160624_I_cinterocdrft001_23';...
    '160624_I_cinterocdrft004_24';...
    
    '160625_I_cinterocdrft001_09';...
    '160627_I_cinterocdrft001_22';...
    '160627_I_cinterocdrft004_22';...
    '160627_I_cinterocdrft004_23';...
    
    '160628_I_cinterocdrft005_20';...
    '160628_I_cinterocdrft005_21';...
    '160628_I_cinterocdrft007_11';...
    '160628_I_cinterocdrft007_13';...
    '160628_I_cinterocdrft008_24';...
    '160628_I_cinterocdrft009_21';...
    
    '160629_I_cinterocdrft001_15';...
    '160629_I_cinterocdrft003_13';...
    
    '160506_I_cinterocdrft002_20';...
    '160506_I_cinterocdrft002_18';...
    
    '160509_I_cinterocdrft004_08';...
    '160509_I_cinterocdrft004_10';...
    
    '160510_I_cinterocdrft005_07';...
    '160510_I_cinterocdrft005_11';...
    '160510_I_cinterocdrft005_13';...
    
    '160512_I_cinterocdrft004_14';...
    '161119_I_cinterocdrft005_23'};



switch unitid
    
    case '160609_I_cinterocdrft013_11'
        
        unit.cinteroc = '160609_I_cinterocdrft013';
        unit.bw       = '160609_I_bwflicker001';
        unit.color    = '160609_I_colorflicker001';
        unit.channel  = 'eD11';
        unit.num      = [1];
        
    case '160609_I_cinterocdrft016_07'
        
        unit.cinteroc = '160609_I_cinterocdrft016';
        unit.bw       = '160609_I_bwflicker002';
        unit.channel  = 'eD07';
        unit.num      = [1];
        
    case '160609_I_cinterocdrft023_17'
        
        unit.cinteroc = '160609_I_cinterocdrft023';
        unit.bw       = '160609_I_bwflicker003';
        unit.channel  = 'eD17';
        unit.num      = [1];
        
    case '160611_I_cinterocdrft003_10'
        
        unit.cinteroc = '160611_I_cinterocdrft003';
        unit.bw       = '160611_I_bwflicker001';
        unit.channel  = 'eD10';
        unit.num      = [1];
        
    case '160613_I_cinterocdrft005_02'
        
        unit.cinteroc = '160613_I_cinterocdrft005';
        unit.bw       = '160613_I_bwflicker001';
        unit.channel  = 'eD02';
        unit.num      = [1];
        
    case '160613_I_cinterocdrft005_14'
        
        unit.cinteroc = '160613_I_cinterocdrft005';
        unit.bw       = '160613_I_bwflicker001';
        unit.channel  = 'eD14';
        unit.num      = [1];
    case '160613_I_cinterocdrft005_15'
        
        unit.cinteroc = '160613_I_cinterocdrft005';
        unit.bw       = '160613_I_bwflicker001';
        unit.channel  = 'eD15';
        unit.num      = [1];
    case '160613_I_cinterocdrft008_09'
        
        unit.cinteroc = '160613_I_cinterocdrft008';
        unit.bw       = '160613_I_bwflicker003';
        unit.channel  = 'eD09';
        unit.num      = [1];
        
    case '160613_I_cinterocdrft009_06'
        
        unit.cinteroc = '160613_I_cinterocdrft009';
        unit.bw       = '160613_I_bwflicker005';
        unit.channel  = 'eD06';
        unit.num      = [1];
        unit.refresh = [60];
        
    case '160616_I_cinterocdrft002_22'
        
        unit.cinteroc = '160616_I_cinterocdrft002';
        unit.bw       = '160616_I_bwflicker001';
        unit.channel  = 'eD22';
        unit.num      = [1];
    case '160616_I_cinterocdrft004_10'
        
        unit.cinteroc = '160616_I_cinterocdrft004';
        unit.bw       = '160616_I_bwflicker002';
        unit.channel  = 'eD10';
        unit.num      = [1];
        unit.refresh  = 60;
        
    case '160616_I_cinterocdrft005_10'
        
        unit.cinteroc = '160616_I_cinterocdrft005';
        unit.bw       = '160616_I_bwflicker003';
        unit.channel  = 'eD10';
        unit.num      = [1];
        unit.refresh  = 60;
        
    case '160616_I_cinterocdrft006_15'
        
        unit.cinteroc = '160616_I_cinterocdrft006';
        unit.bw       = '160616_I_bwflicker004';
        unit.channel  = 'eD15';
        unit.num      = [1];
    case '160623_I_cinterocdrft001_23'
        
        unit.cinteroc = '160623_I_cinterocdrft001';
        unit.bw       = '160623_I_bwflicker001';
        unit.channel  = 'eD23';
        unit.num      = [1];
    case '160623_I_cinterocdrft002_18'
        
        unit.cinteroc = '160623_I_cinterocdrft002';
        unit.bw       = '160623_I_bwflicker002';
        unit.channel  = 'eD18';
        unit.num      = [1];
    case '160623_I_cinterocdrft003_20'
        
        unit.cinteroc = '160623_I_cinterocdrft003';
        unit.bw       = '160623_I_bwflicker003';
        unit.channel  = 'eD20';
        unit.num      = [1];
    case '160623_I_cinterocdrft005_07'
        
        unit.cinteroc = '160623_I_cinterocdrft005';
        unit.channel  = 'eD07';
        unit.num      = [1];
    case '160623_I_cinterocdrft005_07'
        
        unit.cinteroc = '160623_I_cinterocdrft005';
        unit.channel  = 'eD07';
        unit.num      = [1];
    case '160624_I_cinterocdrft001_23'
        
        unit.cinteroc = '160624_I_cinterocdrft001';
        unit.bw       = '160624_I_bwflicker001';
        unit.channel  = 'eD23';
        unit.num      = [1];
    case '160623_I_cinterocdrft002_18'
        
        unit.cinteroc = '160623_I_cinterocdrft002';
        unit.bw       = '160623_I_bwflicker001';
        unit.channel  = 'eD18';
        unit.num      = [1];
    case '160623_I_cinterocdrft003_20'
        
        unit.cinteroc = '160623_I_cinterocdrft003';
        unit.bw       = '160623_I_bwflicker003';
        unit.channel  = 'eD20';
        unit.num      = [1];
    case '160623_I_cinterocdrft003_20'
        
        unit.cinteroc = '160623_I_cinterocdrft003';
        unit.bw       = '160623_I_bwflicker003';
        unit.channel  = 'eD20';
        unit.num      = [1];
    case '160623_I_cinterocdrft005_07'
        
        unit.cinteroc = '160623_I_cinterocdrft005';
        unit.bw       = '160623_I_bwflicker005';
        unit.channel  = 'eD07';
        unit.num      = [1];
    case '160624_I_cinterocdrft001_23'
        
        unit.cinteroc = '160624_I_cinterocdrft001';
        unit.bw       = '160624_I_bwflicker001';
        unit.channel  = 'eD23';
        unit.num      = [1];
    case '160624_I_cinterocdrft004_24'
        
        unit.cinteroc = '160624_I_cinterocdrft004';
        unit.bw       = '160624_I_bwflicker002';
        unit.channel  = 'eD24';
        unit.num      = [1];
    case '160625_I_cinterocdrft001_09'
        
        unit.cinteroc = '160625_I_cinterocdrft001';
        unit.channel  = 'eD09';
        unit.num      = [1];
        unit.refresh  = 60;
        
    case '160627_I_cinterocdrft001_22'
        
        unit.cinteroc = '160627_I_cinterocdrft001';
        unit.bw       = '160627_I_bwflicker001';
        unit.channel  = 'eD22';
        unit.num      = [1];
        unit.refresh  = 60;
    case '160627_I_cinterocdrft004_22'
        
        unit.cinteroc = '160627_I_cinterocdrft004';
        unit.bw       = '160627_I_bwflicker001';
        unit.channel  = 'eD22';
        unit.num      = [1];
        unit.refresh  = [60];
    case '160627_I_cinterocdrft004_23'
        
        unit.cinteroc = '160627_I_cinterocdrft004';
        unit.bw       = '160627_I_bwflicker001';
        unit.channel  = 'eD23';
        unit.num      = [1];
        unit.refresh  = [60];
    case '160628_I_cinterocdrft005_20'
        
        unit.cinteroc = '160628_I_cinterocdrft005';
        unit.bw       = '160628_I_bwflicker002';
        unit.channel  = 'eD20';
        unit.num      = [1];
    case '160628_I_cinterocdrft005_21'
        
        unit.cinteroc = '160628_I_cinterocdrft005';
        unit.bw       = '160628_I_bwflicker002';
        unit.channel  = 'eD21';
        unit.num      = [1];
    case '160628_I_cinterocdrft007_11'
        
        unit.cinteroc = '160628_I_cinterocdrft007';
        unit.bw       = '160628_I_bwflicker003';
        unit.channel  = 'eD11';
        unit.num      = [1];
    case '160628_I_cinterocdrft007_13'
        
        unit.cinteroc = '160628_I_cinterocdrft007';
        unit.bw       = '160628_I_bwflicker003';
        unit.channel  = 'eD13';
        unit.num      = [1];
    case '160628_I_cinterocdrft008_24'
        
        unit.cinteroc = '160628_I_cinterocdrft008';
        unit.bw       = '160628_I_bwflicker004';
        unit.channel  = 'eD24';
        unit.num      = [1];
    case '160628_I_cinterocdrft009_21'
        
        unit.cinteroc = '160628_I_cinterocdrft009';
        unit.bw       = '160628_I_bwflicker005';
        unit.channel  = 'eD21';
        unit.num      = [1];
    case '160629_I_cinterocdrft001_15'
        
        unit.cinteroc = '160629_I_cinterocdrft001';
        unit.bw       = '160629_I_bwflicker001';
        unit.channel  = 'eD15';
        unit.num      = [1];
    case '160629_I_cinterocdrft003_13'
        
        unit.cinteroc = '160629_I_cinterocdrft003';
        unit.bw       = '160629_I_bwflicker003';
        unit.channel  = 'eD13';
        unit.num      = [2];
    case '160506_I_cinterocdrft002_20'
        
        unit.cinteroc = '160506_I_cinterocdrft002';
        unit.channel  = 'eD20';
        unit.num      = [1];
        unit.flipeye  = [1];
        
    case   '160506_I_cinterocdrft002_18'
        
        unit.cinteroc = '160506_I_cinterocdrft002';
        unit.channel  = 'eD18';
        unit.num      = [1];
        unit.flipeye  = [0];
        
    case '160509_I_cinterocdrft004_08'
        
        unit.cinteroc = '160509_I_cinterocdrft004';
        unit.channel  = 'eD08';
        unit.num      = [1];
        unit.flipeye      = [1];
    case '160509_I_cinterocdrft004_10'
        
        unit.cinteroc = '160509_I_cinterocdrft004';
        unit.channel  = 'eD10';
        unit.num      = [1];
    case '160510_I_cinterocdrft005_07'
        
        unit.cinteroc = '160510_I_cinterocdrft005';
        unit.channel  = 'eD07';
        unit.num      = [1];
        unit.flipeye  = 1;
    case '160510_I_cinterocdrft005_11'
        
        unit.cinteroc = '160510_I_cinterocdrft005';
        unit.channel  = 'eD11';
        unit.num      = [1];
    case '160510_I_cinterocdrft005_13'
        
        unit.cinteroc = '160510_I_cinterocdrft005';
        unit.channel  = 'eD13';
        unit.num      = [1];
    case '160512_I_cinterocdrft004_14'
        
        unit.cinteroc = '160512_I_cinterocdrft004';
        unit.channel  = 'eD14';
        unit.num      = [1];
        
    case '161119_I_cinterocdrft005_23';
        
        unit.cinteroc = '161119_I_cinterocdrft005';
        unit.bw = '161119_I_bwflicker002';
        unit.channel  = 'eD23';
        unit.num      = [2];
        
        
        
        
        
end
        
        