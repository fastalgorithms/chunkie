function quiver(obj,varargin)
%QUIVER quiver plot of chunker normal vectors in 2 dimensions



ifhold = ishold();

for i = 1:length(obj)
    tmp = obj(i);
    assert(tmp.dim == 2,'for quiver plot must be 2D chunker');

    tmp = sort(tmp);
    xs = tmp.r(1,:,:); xs = xs(:);
    ys = tmp.r(2,:,:); ys = ys(:);

    rnorms = normals(tmp);

    u = rnorms(1,:,:); u = u(:);
    v = rnorms(2,:,:); v = v(:);

    quiver(xs,ys,u,v,varargin{:});
    hold on
    
end

hold off

if ifhold
    hold on
end
