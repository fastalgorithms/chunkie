function plot(obj,varargin)
%PLOT plot the xy-coordinates of the points of the chunker object
% Must be a 2D chunker. Uses standard plotting commands
%
% Syntax: plot(chnkr,varargin)
%
% Input: 
%   chnkr - chunker object
%   varargin - any of the standard plot commands
%
% Output:
%   none 
%
% Examples:
%   plot(chnkr,'rx') % plots chunker with red x's for points
%   plot(chnkr,'b-','LineWidth',2) % plots chunker with thick blue lines
%
% see also PLOT3, QUIVER

% author: Travis Askham (askhamwhat@gmail.com)

ifhold = ishold();

for i = 1:length(obj)
    tmp = obj(i);
    assert(tmp.dim == 2,'for plot must be 2D chunker');
    [tmp,info] = sort(tmp);
    ifclosed = info.ifclosed;
    nchs = info.nchs;
    istart = 1;
    for ii = 1:length(nchs)
        iend = istart+nchs(ii)-1;
        xs = tmp.r(1,:,istart:iend); xs = xs(:);
        ys = tmp.r(2,:,istart:iend); ys = ys(:);

        if ifclosed
            xs = [xs(:); xs(1)];
            ys = [ys(:); ys(1)];
        end

        if isreal(xs)
            plot(xs,ys,varargin{:})
        else
            plot(real(xs),ys,varargin{:},real(xs),imag(xs),'b--')
        end
        hold on
        istart = istart+nchs(ii);
    end
    
end

hold off

if ifhold
    hold on
end