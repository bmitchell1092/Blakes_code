% visualizing RFs in V1 chamber:
clear all;

grid   = nan(15,15);
center = median([1:size(grid,1)]);

notsampled = 0; 
% setup grid:
if notsampled == 0
grid([center-2:center+2],[1 size(grid,1)])     = 0;

grid([center-4:center+4],[2 size(grid,1)-1])   = 0;

grid([center-5:center+5],[3 (size(grid,1)-2)]) = 0;

grid([center-6:center+6],[4 (size(grid,1)-3)]) = 0;

grid([center-6:center+6],[5 (size(grid,1)-4)]) = 0;

grid(:,[6:10]) = 0;
else
grid([center-2:center+2],[1 size(grid,1)])     = nan;

grid([center-4:center+4],[2 size(grid,1)-1])   = nan;

grid([center-5:center+5],[3 (size(grid,1)-2)]) = nan;

grid([center-6:center+6],[4 (size(grid,1)-3)]) = nan;

grid([center-6:center+6],[5 (size(grid,1)-4)]) = nan;

grid(:,[6:10]) = nan;
end

nangrid = grid;
nangrid(nangrid == 0) = nan; 

% real data points (rfs from recordings)
pts(1).x = 9;
pts(1).y = 9;
pts(1).rx = -7.25;
pts(1).ry = -0.5;

pts(2).x = 10;
pts(2).y = 3;
pts(2).rx = -3.2;
pts(2).ry = -2.3;

pts(3).x = 3;
pts(3).y = 8;
pts(3).rx = -6.5;
pts(3).ry = -4.5;

pts(4).x = 5;
pts(4).y = 12;
pts(4).rx = -8.7;
pts(4).ry = -0.2;

pts(5).x = 13;
pts(5).y = 8;
pts(5).rx = -5.34;
pts(5).ry = -0.15;

pts(6).x = 13;
pts(6).y = 11;
pts(6).rx = -4.7;
pts(6).ry = -2;

pts(7).x = 5;
pts(7).y = 3;
pts(7).rx = -4.3;
pts(7).ry = -4;

pts(8).x = 7;
pts(8).y = 6;
pts(8).rx = -6.8;
pts(8).ry = -3.2;

pts(9).x = 6;
pts(9).y = 8;
pts(9).rx = -7.01;
pts(9).ry = -2.79;

% pts(10).x = 5;
% pts(10).y = 11;
% pts(10).rx = -4.0;
% pts(10).ry = -2.7;


pts(10).x = 4;
pts(10).y = 4;
pts(10).rx = -4.88;
pts(10).ry = -5.15;

pts(11).x = 12;
pts(11).y = 5;
pts(11).rx = -3.6;
pts(11).ry = -1;

pts(12).x = 7;
pts(12).y = 5;
pts(12).rx = -6;
pts(12).ry = -5;

pts(13).x = 9;
pts(13).y = 13;
pts(13).rx = -9;
pts(13).ry = -1.3;

pts(14).x = 11;
pts(14).y = 8;
pts(14).rx = -6;
pts(14).ry = -.1;

gridX = grid;
gridY = grid;
gridEcc = grid;
gridTh = grid;
   
for i = 1:length(pts)
    
    gridX(pts(i).x,pts(i).y) = pts(i).rx;
    
    gridY(pts(i).x,pts(i).y) = pts(i).ry;
    
    [gridth(i),gridecc(i)] = cart2pol(pts(i).rx,pts(i).ry); 
     gridEcc(pts(i).x,pts(i).y) = gridecc(i); 
     gridth(i) = rad2deg(gridth(i)); 
     if gridth(i) < 0
         gridth(i) = gridth(i) + 360; 
     end
     gridTh(pts(i).x,pts(i).y)  = gridth(i); 
     
     nangrid(pts(i).x,pts(i).y)  = gridth(i); 
end

figure, set(gcf,'Color','w'); 
Xmat = flipud(rot90(gridX,1)); 
h = imagesc(Xmat), view([45 90])
set(h,'alphadata',~isnan(gridX)); 
colorbar, set(gca,'Box','off','TickDir','out'); 
title(gca,'x positions'); 
cmap = colormap; cmap(end,:) = [0 0 0]; 
colormap(cmap); cb = colorbar, ylabel(cb,'deg'); 
ylabel('caudal-rostral'); xlabel('medial-lateral'); 

figure, set(gcf,'Color','w'); 
Ymat = flipud(rot90(gridY,1)); 
h = imagesc(Ymat), view([45 90])
set(h,'alphadata',~isnan(gridY));  
colorbar, set(gca,'Box','off','TickDir','out'); 
title(gca,'y positions'); 
cmap = colormap; cmap(end,:) = [0 0 0]; 
colormap(cmap); cb = colorbar, ylabel(cb,'deg');
ylabel('caudal-rostral'); xlabel('medial-lateral'); 
x = get(gca,'XLabel'); rotate(x,[0 0 1],deg2rad(90)); 

figure, set(gcf,'Color',[.5 .5 .5]); 
Emat = flipud(rot90(gridEcc,1)); 
h = imagesc(Emat), view([45 90])
set(h,'alphadata',~isnan(Emat));  
colorbar, set(gca,'Box','off','TickDir','out'); 
title(gca,'eccentricty'); 
cmap = colormap; cmap(1,:) = [0 0 0]; 
colormap(cmap); cb = colorbar, ylabel(cb,'deg');
ylabel('caudal-rostral'); xlabel('medial-lateral'); 
x = get(gca,'XLabel'); rotate(x,[0 0 1],deg2rad(90)); 


figure, set(gcf,'Color',[.5 .5 .5]); 
T2mat =flipud(rot90(nangrid,1)); 
Tmat = flipud(rot90(gridTh,1));  

h = imagesc(Tmat), view([45 90]); 
set(h,'alphadata',~isnan(Tmat)); 
hold on;
h2 = imagesc(T2mat), view([45 90]); 
colorbar, set(gca,'Box','off','TickDir','out'); 
title(gca,'angle'); caxis([179 270]); 
cmap = colormap; 
if notsampled == 0
    cmap(1,:) = [0 0 0];
end
colormap(cmap); cb = colorbar, ylabel(cb,'deg');
ylabel('caudal-rostral'); xlabel('medial-lateral'); 
x = get(gca,'XLabel'); rotate(x,[0 0 1],deg2rad(90)); 
 

% figure,
% [TH,R] = meshgrid(linspace(0,2*pi,100),linspace(1,10,100));
% [X,Y] = pol2cart(TH,log(R));
% surface(X,Y,zeros(size(X)),R);
% shading interp, colorbar

[X, Y] = ndgrid(1:15,1:15); 
V      = Emat; 
F      = griddedInterpolant(X,Y,V,'linear'); 
