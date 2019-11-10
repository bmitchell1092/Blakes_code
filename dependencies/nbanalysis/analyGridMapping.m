% analyze grid mapping data
clear all;

drname = 'C:\Users\MLab\Desktop\receptive-field-mapping-master';
cd(drname);

fkey = '*mods*.gGratingXY';
fname = dir(sprintf('%s',fkey));

ana = readFGXY(fname.name);

horzpos = round(ana.horzdva*10)./10; 
vertpos = round(ana.vertdva*10)./10; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot positions of presented stimuli:
figure,
set(gcf,'Color','w');
plot(ana.horzdva,ana.vertdva,'bo','MarkerFaceColor','b');
set(gca,'box','off','LineWidth',2,'TickDir','out');
title('Positions of presented stimuli');
xlabel('horizontal position (dva)');
ylabel('vertical position (dva)');

[xpos,ypos] = pol2cart(ana.theta,ana.eccentricity);
figure,
plot(xpos,ypos,'bo');
title('cartesian coordinates of presented stimuli');

figure, polar(ana.theta,ana.eccentricity,'ro');
title('polar coordinates of presented stimuli');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot positions of presented stimuli color-coded by other stimulus parameter:
param = 'tilt';

switch param
    
    case 'sf'
        
        G = grp2idx(ana.sf);
        figure,
        map = winter(max(G)-1);
        map = [0 0 0; map]; hold on;
        xlabel('x position');
        ylabel('y position');
        colormap(map); caxis([min(ana.sf) max(ana.sf)]); c = colorbar;
        ylabel(c,'spatial frequency');
        scatter(xpos,ypos,[],map(G,:),'o');
        
    case 'tilt'
        G = grp2idx(ana.tilt);
        figure,
        map = winter(max(G)-1);
        map = [0 0 0; map]; hold on;
        xlabel('x position');
        ylabel('y position');
        colormap(map); caxis([min(ana.tilt) max(ana.tilt)]); c = colorbar;
        ylabel(c,'tilt');
        scatter(xpos,ypos,[],map(G,:),'o');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot positions of presented stimuli color-coded by other stimulus parameter:
coord = cat(2,round(ana.horzdva*10)./10,round(ana.vertdva*10)/10);
unqvals  = unique(coord,'rows'); 
npos     = size(unique(coord)); 
possible = length([-90:10:0:10:30]) .* length([3:0.5:6]);

possf   = [0.2:0.2:0.6]; 
postilt = [0:30:359]; 

emp_sf = nan(length(unqvals),length(possf)); 
emp_tilt = nan(length(unqvals),length(postilt)); 

% sort trials by position
for i = 1:length(unqvals)
    % trials at these coords
    thesetrs = find(coord(:,1) == unqvals(i,1) & coord(:,2) == unqvals(i,2)); 
    ntrsatpos(i) = length(thesetrs); 
    
    %sf at this position on these trials:
    emp_sf(i,:) = hist(ana.sf(thesetrs),possf); 
    
    %tilt at this position on these trials: 
    emp_tilt(i,:) = hist(ana.tilt(thesetrs),postilt); 
    
end

figure, 
set(gcf,'Color','w'); 
bar([1:length(unqvals)],ntrsatpos); 
xlabel('position number'); 
ylabel('n total trials'); 
box off

figure, 
set(gcf,'Color','w','Position',[1 1 700 600],'PaperPositionMode','auto'); 
bar(emp_sf,'stacked');
xlabel('position number'); 
ylabel('n trials'); 
legend(gca,'sf 0.2','sf 0.4','sf 0.6'); 
box off
title('n trials per spatial freq by position'); 

figure, 
set(gcf,'Color','w','Position',[1 1 900 600],'PaperPositionMode','auto'); 
bar(emp_tilt,'stacked');
xlabel('position number'); 
ylabel('n trials'); 
title('n trials per orientation by position'); 
legend(gca,'0','30','60','90','120','150','180','210','240','270','300','330'); 
box off




