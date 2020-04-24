function quiver(obj,varargin)
%QUIVER quiver plot of chunker normal vectors in 2 dimensions
% Uses standard quiver commands
%
% Syntax: quiver(chnkr,varargin)
%
% Input: 
%   chnkr - chunker object
%   varargin - any of the standard quiver commands
%
% Output:
%   none 
%
% Examples:
%   quiver(chnkr,'r','LineWidth',2) % quiver plot of chunker normals
%                                   % with thick red arrows
%
% see also PLOT, PLOT3

% author: Travis Askham (askhamwhat@gmail.com)

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
