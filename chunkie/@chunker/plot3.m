function plot3(obj,idata,varargin)
%PLOT3 plot values in the specified data row as the z-coordinate with 
% the xy-coordinates the nodes of the chunker object. Must be a 2D chunker.
% Uses standard plotting commands
%
% Syntax: plot3(chnkr,idata,varargin)
%
% Input: 
%   chnkr - chunker object
%   idata - chunker data row to plot
%   varargin - any of the standard plot commands
%
% Output:
%   none 
%
% Examples:
%   % add a data row and plot it
%   chnkr = chnkr.makedatarows(1);
%   fun1 = sin(chnkr.r(1,:));
%   chnkr.data(1,:) = fun1;
%   plot3(chnkr,1,'rx') % plots curve with red x's for points
%   plot3(chnkr,1,'b-','LineWidth',2) % plots curve with thick blue lines
%
% see also PLOT, QUIVER

% author: Travis Askham (askhamwhat@gmail.com)

ifhold = ishold();
if nargin < 2
    idata = 1;
end

for i = 1:length(obj)
    tmp = obj(i);
    assert(tmp.hasdata);
    assert(tmp.dim == 2,'for plot must be 2D chunker');
    [tmp,info] = sort(tmp);
    ifclosed = info.ifclosed;
    nchs = info.nchs;
    istart = 1;
    for ii = 1:length(nchs)
        iend = istart+nchs(ii)-1;
        xs = tmp.r(1,:,istart:iend); xs = xs(:);
        ys = tmp.r(2,:,istart:iend); ys = ys(:);
        zs = tmp.data(idata,:,istart:iend); zs = zs(:);
    
        if ifclosed
            xs = [xs(:); xs(1)];
            ys = [ys(:); ys(1)];
            zs = [zs(:); zs(1)];
        end

        plot3(xs,ys,zs,varargin{:})
        hold on
        istart = istart+nchs(ii);
    end
    
end

hold off

if ifhold
    hold on
end