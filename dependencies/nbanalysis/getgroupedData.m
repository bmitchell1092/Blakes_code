% datagroupfileslist: 
% datagroupfileslist: 

function [dates,files,stpath,kiloflick,rm] = getgroupedData
loc2 = '/volumes/drobo2/data/neurophys/rig021'; 
loc1 = '/volumes/drobo/data/neurophys/rig021'; 
rm = []; 

%kilo flick--visually responsive units detected by kilosort--run on
%concated bw/color flicker files
    
dates{1}  = '160609';
files{1}  = {'cinterocdrft009','bwflicker001','colorflicker003','cinterocdrft013'};
stpath{1} = {loc1}; 
kiloflick{1}.chan  = [11 11]; 
kiloflick{1}.clust = [56 49];

 
dates{length(dates)+1}  = '160609';
files{length(dates)}  = {'cinterocdrft016','bwflicker002','colorflicker004'};
stpath{length(dates)} = {loc1}; 
kiloflick{length(dates)}.chan  = [7]; 
kiloflick{length(dates)}.clust = [];
    
dates{length(dates)+1}  = '160609';
files{length(dates)}  = {'cinterocdrft023','bwflicker003','colorflicker005'};
stpath{length(dates)} = {loc1}; 
kiloflick{length(dates)}.chan  = [6 17]; 
kiloflick{length(dates)}.clust = [28 2];

dates{length(dates)+1}  = '160611';
files{length(dates)}  = {'cinterocdrft001','bwflicker001','colorflicker001'};
stpath{length(dates)} = {loc1}; 
kiloflick{length(dates)}.chan  = [13 14 15]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160611';
files{length(dates)}  = {'cinterocdrft003','bwflicker002','colorflicker002'};
stpath{length(dates)} = {loc1}; 
kiloflick{length(dates)}.chan  = [10]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160613';
files{length(dates)}  = {'cinterocdrft005','bwflicker001','colorflicker002'};
stpath{length(dates)} = {loc2}; 
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160613';
files{length(dates)}  = {'cinterocdrft008','bwflicker003','colorflicker003'};
stpath{length(dates)} = {loc2}; 
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160613';
files{length(dates)}  = {'cinterocdrft009','bwflicker005','colorflicker004'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [6]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160615';
files{length(dates)}  = {'cinterocdrft001','cinterocdrft002','bwflicker001','colorflicker001','tfsfdrft001'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160615';
files{length(dates)}  = {'cinterocdrft005','bwflicker002','colorflicker002'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [21 21 22]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160616';
files{length(dates)}  = {'cinterocdrft002','bwflicker001','colorflicker001'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [22]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160616';
files{length(dates)}  = {'cinterocdrft004','bwflicker002','colorflicker002'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [8 10]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160616';
files{length(dates)}  = {'cinterocdrft005','bwflicker003','colorflicker003'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [8 10]; 
kiloflick{length(dates)}.clust = [72 25];

dates{length(dates)+1}  = '160616';
files{length(dates)}  = {'cinterocdrft006','bwflicker002','colorflicker004'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = []; %have this for 16016_I_005
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160623';
files{length(dates)}  = {'cinterocdrft001','bwflicker001','colorflicker001'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160623';
files{length(dates)}  = {'cinterocdrft002','bwflicker002','colorflicker002'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [18 20]; 
kiloflick{length(dates)}.clust = [58 40];

dates{length(dates)+1}  = '160623';
files{length(dates)}  = {'cinterocdrft003','colorflicker003','bwflicker003'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [20]; 
kiloflick{length(dates)}.clust = [3];

dates{length(dates)+1}  = '160623';
files{length(dates)}  = {'cinterocdrft005','colorflicker004','bwflicker005'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [7]; 
kiloflick{length(dates)}.clust = [3];

dates{length(dates)+1}  = '160623';
files{length(dates)}  = {'cinterocdrft006','colorflicker005','bwflicker006'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [23]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160624';
files{length(dates)}  = {'cinterocdrft001','colorflicker001','bwflicker001'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160624';
files{length(dates)}  = {'cinterocdrft004','colorflicker002','bwflicker002'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160625';
files{length(dates)}  = {'cinterocdrft001','colorflicker001','bwflicker001'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160627';
files{length(dates)}  = {'cinterocdrft007','cinterocdrft008','cinterocdrft013','colorflicker004','bwflicker004'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [22]; 
kiloflick{length(dates)}.clust = [5];

dates{length(dates)+1}  = '160627';
files{length(dates)}  = {'cinterocdrft001','colorflicker001','bwflicker001'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160627';
files{length(dates)}  = {'cinterocdrft004','colorflicker002','bwflicker001'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = []; %[22 22]; 
kiloflick{length(dates)}.clust = []; %[25 2];

dates{length(dates)+1}  = '160628';
files{length(dates)}  = {'cinterocdrft004','colorflicker001','bwflicker001','cinterocdrft001'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [15 20 21]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160628';
files{length(dates)}  = {'cinterocdrft005','colorflicker002','bwflicker002'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [21]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160628';
files{length(dates)}  = {'cinterocdrft007','colorflicker003','bwflicker003'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [11 13 13]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160628';
files{length(dates)}  = {'cinterocdrft008','colorflicker004','bwflicker004'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160628';
files{length(dates)}  = {'cinterocdrft009','colorflicker005','bwflicker005'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [21]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160629';
files{length(dates)}  = {'cinterocdrft001','colorflicker001','bwflicker001'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [15]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160629';
files{length(dates)}  = {'cinterocdrft003','colorflicker003','bwflicker003'};
stpath{length(dates)} = {loc2};
kiloflick{length(dates)}.chan  = [13]; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160505';
files{length(dates)}  = {'cinterocdrft007','colorflicker001','cinterocdrft010'};
stpath{length(dates)} = {loc1};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160506';
files{length(dates)}  = {'cinterocdrft002','colorflicker001','cinterocdrft005'};
stpath{length(dates)} = {loc1};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160509';
files{length(dates)}  = {'cinterocdrft002','colorflicker002','colorflicker002'};
stpath{length(dates)} = {loc1};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160510';
files{length(dates)}  = {'cinterocdrft001','cinterocdrft005','colorflicker001'};
stpath{length(dates)} = {loc1};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160512';
files{length(dates)}  = {'cinterocdrft004','colorflicker001'};
stpath{length(dates)} = {loc1};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160517';
files{length(dates)}  = {'cinterocdrft003','colorflicker001'};
stpath{length(dates)} = {loc1};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];

dates{length(dates)+1}  = '160602';
files{length(dates)}  = {'cinterocdrft004','colorflicker001'};
stpath{length(dates)} = {loc1};
kiloflick{length(dates)}.chan  = []; 
kiloflick{length(dates)}.clust = [];
