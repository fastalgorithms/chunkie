
function [r,d,d2] = starfish(t,varargin)
%STARFISH
% return position, first and second derivatives of parameterized starfish
% domain. 
%
% Inputs:
% t - points (in [0,2pi] to evaluate these quantities
%
% Optional inputs:
% narms - integer, number of arms on starfish (5)
% amp - float, amplitude of starfish arms relative to radius of length 1
%               (0.3)
% ctr - float(2), x,y coordinates of center of starfish ( [0,0] )
% phi - float, phase shift (0)
% scale - scaling factor (1.0)

narms = 5;
amp = 0.3;
x0 = 0.0;
y0 = 0.0;
phi = 0.0;
scale = 1.0;
if nargin > 1 && ~isempty(varargin{1})
    narms = varargin{1};
end
if nargin > 2 && ~isempty(varargin{2})
    amp = varargin{2};
end
if nargin > 3 && ~isempty(varargin{3})
    ctr = varargin{3};
    x0 = ctr(1); y0 = ctr(2);
end
if nargin > 4 && ~isempty(varargin{4})
    phi = varargin{4};
end
if nargin > 5 && ~isempty(varargin{5})
    scale = varargin{5};
end


fvals = zeros(length(t),6);
ct = cos(t);
st = sin(t);
cnt = cos(narms*(t+phi));
snt = sin(narms*(t+phi));

xs = x0+(1+amp*cnt).*ct*scale;
ys = y0+(1+amp*cnt).*st*scale;
dxs = -(1+amp*cnt).*st-narms*amp*snt.*ct;
dxs = dxs*scale;
dys = (1+amp*cnt).*ct-narms*amp*snt.*st;
dys = dys*scale;
d2xs = -dys-narms*amp*(narms*cnt.*ct-snt.*st);
d2xs = d2xs*scale;
d2ys = dxs-narms*amp*(narms*cnt.*st+snt.*ct);
d2ys = d2ys*scale;

r = [(xs(:)).'; (ys(:)).'];
d = [(dxs(:)).'; (dys(:)).'];
d2 = [(d2xs(:)).'; (d2ys(:)).'];

end

