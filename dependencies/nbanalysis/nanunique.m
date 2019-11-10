function u = nanunique(varargin)
  x = varargin{1};
    t = rand;
      while any(x(:)==t)
                t = rand;
      end
        x(isnan(x)) = t;
          u = unique(x,varargin{2:end});
            u(u==t)=[];
end