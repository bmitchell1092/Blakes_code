function h = mybreakplot(x,y,y_break_start,y_break_end);

% SOME DFINITIONS
y_break_width = y_break_end - y_break_start;
y_break_mid   = y_break_width./2 + y_break_start;
y_range       = range(y(:,1));

% LOSE THE DATA IN THE BREAK, WE DON'T NEED IT ANYMORE
i   = y > y_break_start & y < y_break_end;
x(i)= [];
y(i)= [];

% MAP THE DATA
i    = y >= y_break_end;
y2   = y - i.*y_break_width;

% PLOT THE MAPPED DATA
h    = plot(x,y,'-'); hold on; 
ylim = get(gca,'ylim');
h    = plot(x,y2,'-');
set(gca,'ylim',ylim-[0 y_break_width]);

% CREATE THE "BREAK" EFFECT
xlim        = get(gca,'xlim');
xtick       = get(gca,'XTick');
ytick       = get(gca,'YTick');
yticklabel  = get(gca,'yticklabel');

y_gap_width = y_range;
y_half_gap  = y_gap_width./2;
y_gap_mid   = y_break_start + y_half_gap;

x_half_tick = abs(xtick(1));
xx          = [xlim(1) xlim(1)+x_half_tick];
set(gca,'xlim',xlim);

% MAP TICKS BACK
i_wrong_ticks = ytick > y_break_start;
ytick = ytick + i_wrong_ticks.*y_break_width;
integer_ticks = all(floor(ytick) == ytick);
label_width = size(yticklabel,2);
if integer_ticks
    format_string = sprintf('%%%dd\n',label_width);
else
    left_side = ceil(log10(max(ytick)));
    right_side = label_width-left_side-1;
    format_string = sprintf('%%%d.%df\n',label_width,right_side);
end;
set(gca, 'yticklabel', num2str(ytick'));

hl = plot(xlim(1):xlim(2),repmat(y_break_start,length(xlim(1):xlim(2)),1)); hold on; set(hl,'color','r'); 
hl2 = plot(xlim(1):xlim(2),repmat(y_break_end+.001,length(xlim(1):xlim(2)),1)); hold on; set(hl2,'color','b');
