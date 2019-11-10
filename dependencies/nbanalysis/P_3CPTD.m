function P_3CPTD(data1, data2, tim, xlim, ylim, cmap, vl, hl)

data1 = H_2DSMOOTH(data1);
data2 = H_2DSMOOTH(data2);

limi1 = nanmax(nanmax(abs(data1(ylim(1):ylim(2),:))));
limi2 = nanmax(nanmax(abs(data2(ylim(1):ylim(2),:))));

limi = max([limi1 limi2]);

figure;
set(gcf, 'units', 'normalized', 'position', [.025 .025, .45, .9])
subplot(3,1,1)
imagesc(tim, 1:size(data1,1), data1);
box off
title('Target in RF')
xlabel('Time (ms)')
ylabel('Depth from L4/5 Boundary (mm)')
set(gca,'xlim', xlim, 'ylim', ylim, ...
    'linewidth', 2, 'fontsize', 12, 'fontweight', 'bold', ...
    'ytick', ylim(1):50:ylim(2), 'yticklabel',num2str(5-(([ylim(1):50:ylim(2)])./100)'))
caxis([-limi limi]);
v1 = vline(vl);
h1 = hline(hl);
set(v1, 'color', 'k', 'linewidth', 2, 'linestyle', '-.')
set(h1, 'color', 'k', 'linewidth', 2, 'linestyle', '-.')
c1 = colorbar;

subplot(3,1,2)
imagesc(tim, 1:size(data2,1), data2);
box off
title('Distractor in RF')
xlabel('Time (ms)')
ylabel('Depth from L4/5 Boundary (mm)')
set(gca,'xlim', xlim, 'ylim', ylim, ...
    'linewidth', 2, 'fontsize', 12, 'fontweight', 'bold', ...
    'ytick', ylim(1):50:ylim(2), 'yticklabel',num2str(5-(([ylim(1):50:ylim(2)])./100)'))
caxis([-limi limi]);
v1 = vline(vl);
h1 = hline(hl);
set(v1, 'color', 'k', 'linewidth', 2, 'linestyle', '-.')
set(h1, 'color', 'k', 'linewidth', 2, 'linestyle', '-.')
c1 = colorbar;

subplot(3,1,3)
imagesc(tim, 1:size(data1,1), data1-data2);
box off
title('Target Selection')
xlabel('Time (ms)')
ylabel('Depth from L4/5 Boundary (mm)')
set(gca,'xlim', xlim, 'ylim', ylim, ...
    'linewidth', 2, 'fontsize', 12, 'fontweight', 'bold', ...
    'ytick', ylim(1):50:ylim(2), 'yticklabel',num2str(5-(([ylim(1):50:ylim(2)])./100)'))
if strcmp('tej', cmap)
    tej = flipud(colormap('jet'));
    colormap(tej)
else
    colormap(cmap)
end
v1 = vline(vl);
h1 = hline(hl);
set(v1, 'color', 'k', 'linewidth', 2, 'linestyle', '-.')
set(h1, 'color', 'k', 'linewidth', 2, 'linestyle', '-.')
c1 = colorbar;

end