function plot(obj,varargin)
%
%   '-x' for plots with points marked by 'x's with a different color
%        for each chunker object.
%
%

ifhold = ishold();

echnks =  obj.echnks;

hold on
for i=1:numel(echnks)
    plot(echnks(i),varargin{:});
end

hold off

if ifhold
    hold on
end