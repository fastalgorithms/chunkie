function cg = geo_spiral(zk, ppw)
%GEO_SPIRAL chunkgraph for the Bruno-Lintner 2012 logarithmic spiral
% gamma(t) = exp(-1-5i + (2+10i) t), t in [0,1].
% Endpoints: gamma(0) = exp(-1-5i), gamma(1) = exp(1+5i).

if nargin < 2 || isempty(ppw), ppw = 10; end
pref = []; pref.k = 16;
maxchunklen = pref.k/ppw/abs(zk)*2*pi;

a = exp(-1 - 5i);
b = 2 + 10i;
fchnk = @(t) spiralcurve(t,a,b);

z0 = a;            % gamma(0)
z1 = a*exp(b);     % gamma(1)
verts = [real(z0), real(z1); imag(z0), imag(z1)];
edge2verts = [1;2];

cparams = struct('ta',0,'tb',1,'maxchunklen',maxchunklen,'nover',2);
cg = chunkgraph(verts, edge2verts, {fchnk}, cparams, pref);
cg = balance(cg);
end

function [r,d,d2] = spiralcurve(t,a,b)
% z(t)=a*exp(b*t),  z'=a*b*exp(b*t),  z''=a*b^2*exp(b*t)
ebt = exp(b*t(:).');
z   = a*ebt;       zp  = a*b*ebt;     zpp = a*b*b*ebt;
r   = [real(z);  imag(z)];
d   = [real(zp); imag(zp)];
d2  = [real(zpp);imag(zpp)];
end
