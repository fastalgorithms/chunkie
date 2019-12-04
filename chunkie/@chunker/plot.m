function plot(obj,varargin)
%PLOT plot the coordinates of the chunker object
% Must be a 2D chunker

ifhold = ishold();

for i = 1:length(obj)
    tmp = obj(i);
    assert(tmp.dim == 2,'for plot must be 2D chunker');
    [tmp,ifclosed] = sort(tmp);
    xs = tmp.r(1,:,:); xs = xs(:);
    ys = tmp.r(2,:,:); ys = ys(:);

    if ifclosed
        xs = [xs(:); xs(1)];
        ys = [ys(:); ys(1)];
    end

    plot(xs,ys,varargin{:})
    hold on
    
end

hold off

if ifhold
    hold on
end