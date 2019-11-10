function rocplot(falsealarms,hits)
% rocplot(falsealarms,hits)
% plots ROC curves
% a) with 'plot' (smooth)
% b) with 'stairs' (sawtooth pattern)
% USAGE: rocplot(falsealarms,hits)
% INPUT: falsealarms, hits - must be normalized to 1!
% NOTES: 03.05.04 - AM
 
figure(99)
subplot(1,2,1)
set(gca,'FontSize',6)
plot(falsealarms,hits,'k.-')
hold on
xlim([0 1])
ylim([0 1])
plot([0 1],[0 1],'k:')
axis square
subplot(1,2,2)
set(gca,'FontSize',6)
stairs(falsealarms,hits,'k.-')
hold on
xlim([0 1])
ylim([0 1])
plot([0 1],[0 1],'k:')
axis square
 