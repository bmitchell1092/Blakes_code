function P_EVPLINE_BASIC_2(data_in, data_in2, t_in, varargin)

c1 = 'k';
c2 = 'k';
s1 = '-';
s2 = ':';

tst = [];

t1 = -50;
t2 = 250;

varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-c1'}
            c1 = varargin{varStrInd(iv)+1};
        case {'-c2'}
            c2 = varargin{varStrInd(iv)+1};
        case {'-s1'}
            s1 = varargin{varStrInd(iv)+1};
        case {'-s2'}
            s2 = varargin{varStrInd(iv)+1};
        case {'-tst'}
            tst = varargin{varStrInd(iv)+1};
        case {'-tocx'}
            tocx = varargin{varStrInd(iv)+1};
        case {'-bocx'}
            bocx = varargin{varStrInd(iv)+1};
        case {'-tol4'}
            tol4= varargin{varStrInd(iv)+1};
        case {'-bol4'}
            bol4 = varargin{varStrInd(iv)+1};
        case {'-t1'}
            t1 = varargin{varStrInd(iv)+1};
        case {'-t2'}
            t2 = varargin{varStrInd(iv)+1};
        case {'-ci1'}
            ci1 = varargin{varStrInd(iv)+1};
        case {'-ci2'}
            ci2 = varargin{varStrInd(iv)+1};
    end
end

t_min = nanmin(nanmin(data_in));
t_max = nanmax(nanmax(data_in));

t_min2 = nanmin(nanmin(data_in2));
t_max2 = nanmax(nanmax(data_in2));

sep = max([t_max t_max2]) - max([t_min t_min2])*.5;

hold on;

if exist('tol4', 'var') && exist('bol4', 'var')
    
    fill([min(t_in) min(t_in) max(t_in) max(t_in)], ...
        [ repmat(nanmean(data_in(bol4,t_in<0) - (sep*(bol4-1)),2), numel(t_in), 1) - sep * .5, ...
        repmat(nanmean(data_in(tol4,t_in<0) - (sep*(tol4-1)),2), numel(t_in), 1) + sep * .5, ...
        repmat(nanmean(data_in(tol4,t_in<0) - (sep*(tol4-1)),2), numel(t_in), 1) + sep * .5, ...
        repmat(nanmean(data_in(bol4,t_in<0) - (sep*(bol4-1)),2), numel(t_in), 1) - sep * .5], ...
        [ .9 .9 .9 ], 'linestyle', 'none');
    
end

if exist('tocx', 'var')
    plot(t_in,  repmat(nanmean(data_in(tocx,t_in<0) - (sep*(tocx-1)),2), numel(t_in), 1) + ...
        sep*.5, ...
        'linewidth', 2, 'linestyle', '--', 'color', 'k')
end

if exist('bocx', 'var')
    plot(t_in,  repmat(nanmean(data_in(bocx,t_in<0) - (sep*(tobx-1)),2), numel(t_in), 1) - ...
        sep*.5, ...,
        'linewidth', 2, 'linestyle', '--', 'color', 'k')
end

if exist('tol4', 'var')
    plot(t_in,  repmat(nanmean(data_in(tol4,t_in<0) - (sep*(tol4-1)),2), numel(t_in), 1) + ...
        sep*.5, ...,
        'linewidth', 2, 'linestyle', '--', 'color', 'k')
end

if exist('bol4', 'var')
    plot(t_in,  repmat(nanmean(data_in(bol4,t_in<0) - (sep*(bol4-1)),2), numel(t_in), 1) - ...
        sep*.5, ...,
        'linewidth', 2, 'linestyle', '--', 'color', 'k')
end

if size(data_in, 1) > 1
    for i_ch = 1 : size( data_in, 1)
        
        plot(t_in, data_in(i_ch,:) - (sep*(i_ch-1)), ...
            'linewidth', 2, 'color', c1, 'linestyle', s1)
        
        plot(t_in, data_in2(i_ch,:) - (sep*(i_ch-1)), ...
            'linewidth', 2, 'color', c2, 'linestyle', s2)
        
        if ~isempty(tst) && ~isnan(tst(i_ch))
            
            if i_ch == 1
                plot([tst(i_ch) tst(i_ch)], ...
                    [nanmean(data_in(i_ch, t_in < 0) - (sep*(i_ch-1)), 2), ...
                    0 + sep], 'color', 'r', 'linewidth', 2)
            else
                plot([tst(i_ch) tst(i_ch)], ...
                    [nanmean(data_in(i_ch, t_in < 0) - (sep*(i_ch-1)), 2), ...
                    nanmean(data_in(i_ch, t_in < 0) - (sep*(i_ch-2)), 2)], ...
                    'color', 'r', 'linewidth', 2)
            end
            
        end
        
    end
    
else
    
    plot(t_in, data_in, ...
        'linewidth', 2, 'color', c1, 'linestyle', s1)
    
    plot(t_in, data_in2, ...
        'linewidth', 2, 'color', c2, 'linestyle', s2)
    
    if exist('ci1', 'var')
        ciplot(ci1(1,:), ci1(2,:), t_in)
    end
    
    if exist('ci2', 'var')
        ciplot(ci2(1,:), ci2(2,:), t_in)
    end
    
    if ~isempty(tst) && ~isnan(tst)
        tv = vline(tst);
        set(tv, 'color', 'r', 'linewidth', 2);  
    end
    
end

vl = vline(0);
set(vl, 'color', 'k', 'linestyle', '-', 'linewidth', 2);

set(gca, 'xlim', [t1 t2], ...
    'linewidth', 2, ...
    'fontsize', 14, 'fontweight', 'bold', ...
    'layer', 'top')

if size(data_in, 1) > 1
    set(gca, 'ylim', [-sep*(i_ch), 0 + sep], 'ytick', []);
    ylabel('Deep <--- Depth ---> Shallow')
end

box off;

xlabel('Time (ms)');

end