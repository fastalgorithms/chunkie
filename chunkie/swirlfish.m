function [r,d,d2] = swirlfish(t,varargin)
%SWIRLFISH return position, first and second derivatives of a starfish 
% domain whose arms curl in a pinwheel pattern. The arms are no longer 
% mirror symmetric, but the domain keeps its narms-fold rotational 
% symmetry.
%
% The curve is the starfish pushed through the twist map 
% (rho,theta) -> (rho, theta + twist*(rho-1)), which rotates each circle 
% of radius rho by the angle twist*(rho-1):
%
% rho(t)   = 1 + amp*cos(narms*(t+phi))
% theta(t) = t + twist*(rho(t)-1)
% x(t) = x0 + scale*rho(t)*cos(theta(t))
% y(t) = y0 + scale*rho(t)*sin(theta(t))
%
% The twist map is an area-preserving homeomorphism of the plane, so for 
% any value of twist the curve is simple and encloses the same area as the 
% starfish, pi*(1+amp^2/2)*scale^2. When abs(twist)*amp*narms < 1, theta 
% is increasing in t and the domain is star-shaped about ctr. Larger 
% twists curl the arms further and narrow the gaps between them, so more 
% chunks are needed to resolve the boundary.
%
% Syntax: [r,d,d2] = swirlfish(t,narms,amp,twist,ctr,phi,scale)
%
% Input:
%   t - array of points (in [0,2pi])
%
% Optional input:
%   narms - integer, number of arms on swirlfish (5)
%   amp - float, amplitude of arms relative to radius of length 1 (0.3)
%   twist - float, rotation angle (in radians) per unit change in radius 
%           (1.0). each arm sweeps through an angle 2*amp*twist from base 
%           to tip. for twist > 0 the tips are rotated counterclockwise 
%           relative to the base; twist < 0 gives the mirror image and 
%           twist = 0 gives the starfish
%   ctr - float(2), x0,y0 coordinates of center of swirlfish ( [0,0] )
%   phi - float, phase shift (0)
%   scale - scaling factor (1.0)
%
% Output:
%   r - 2 x numel(t) array of positions, r(:,i) = [x(t(i)); y(t(i))]
%   d - 2 x numel(t) array of t derivative of r 
%   d2 - 2 x numel(t) array of second t derivative of r 
%
% Examples:
%   [r,d,d2] = swirlfish(t); % get default settings 
%   chnkr = chunkerfunc(@(t) swirlfish(t,5,0.3,2.0)); % stronger curl
%   [r,d,d2] = swirlfish(t,narms,[],[],ctr,[],scale); % change some settings
%
% see also STARFISH

narms = 5;
amp = 0.3;
twist = 1.0;
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
    twist = varargin{3};
end
if nargin > 4 && ~isempty(varargin{4})
    ctr = varargin{4};
    x0 = ctr(1); y0 = ctr(2);
end
if nargin > 5 && ~isempty(varargin{5})
    phi = varargin{5};
end
if nargin > 6 && ~isempty(varargin{6})
    scale = varargin{6};
end

cu = cos(narms*(t+phi));
su = sin(narms*(t+phi));

% radius and its derivatives
rho = 1+amp*cu;
drho = -narms*amp*su;
d2rho = -narms^2*amp*cu;

% angle and its derivatives
th = t+twist*amp*cu;
dth = 1+twist*drho;
d2th = twist*d2rho;

ct = cos(th);
st = sin(th);

xs = x0+scale*rho.*ct;
ys = y0+scale*rho.*st;
dxs = scale*(drho.*ct-rho.*dth.*st);
dys = scale*(drho.*st+rho.*dth.*ct);
d2xs = scale*((d2rho-rho.*dth.^2).*ct-(2*drho.*dth+rho.*d2th).*st);
d2ys = scale*((d2rho-rho.*dth.^2).*st+(2*drho.*dth+rho.*d2th).*ct);

r = [(xs(:)).'; (ys(:)).'];
d = [(dxs(:)).'; (dys(:)).'];
d2 = [(d2xs(:)).'; (d2ys(:)).'];

end
