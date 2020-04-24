function trap = trapfunc(fcurve,cparams,pref)
%TRAPFUNC discretize a curve using uniform points in 
% parameter space. 
%
%	cparams = curve parameters (default value)
%
%	cparams.ta = left end of t interval (0)
% 	cparams.tb = right end of t interval (2*pi)
%   cparams.npt = number of points to use (100)
%
% See also TRAPPER

ta = 0.0; tb = 2*pi;
npt = 100;

if isfield(cparams,'ta')
    ta = cparams.ta;
end	 
if isfield(cparams,'tb')
    tb = cparams.tb;
end	 
if isfield(cparams,'npt')
    npt = cparams.npt;
end

dim = checkcurveparam(fcurve,ta);

if nargin < 2
    cparams = [];
end
if nargin < 3
    p = []; p.dim = dim;
    pref = trapperpref(p);
else
    pref = trapperpref(pref);
end

assert(pref.dim == dim);


nout = 3;
out = cell(nout,1);

trap = trapper(pref);
trap = trap.addpt(npt);

h = (tb-ta)/npt;
ts = ta + (0:(npt-1))*h;

[rs,ds,d2s] = fcurve(ts);
trap.r = rs;
trap.d = ds;
trap.d2 = d2s;
trap.h = h;

end


