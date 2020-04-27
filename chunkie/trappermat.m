function [sysmat] = trappermat(trap,kern,opts)
%TRAPPERMAT build matrix for given kernel and trapper description of 
% boundary. This is a wrapper for various quadrature routines
%
% Input:
%   trap - trapper object describing boundary
%   kern  - kernel function. By default, this should be a function handle
%           accepting input of the form kern(s,t,taus,taut), where s and t
%           are the source and target locations and taus and taut are the 
%           unit tangent at the source and target locations.
%   opts  - options structure. available options (default settings)
%           opts.quad = string ('balog'), specify quadrature routine to 
%                       use. 
%
%                       - 'balog' uses an Alpert-style rule for
%                          logarithmically singular kernels and 
%                       smooth kernels with removable singularities
%                       - 'smooth' selects standard smooth trapezoidal
%                       quadrature
%
%           opts.quadorder = integer (8), desired quadrature order.
%           opts.nonsmoothonly = boolean (false), if true, only compute the
%                         entries for which a special quadrature is used
%                         (e.g. self and neighbor interactoins) and return
%                         in a sparse array.
%
% Output:
%   sysmat - the system matrix for convolution of the kernel defined by
%            kern with a density on the domain defined by trap

if nargin < 3
    opts = [];
end

quadorder = 8;
quad = 'balog';
nonsmoothonly = false;

if trap.npt < 2
  warning('trapper description has 1 or fewer points, doing nothing');
  sysmat = zeros(trap.npt,trap.npt);
  return
end

% determine operator dimensions using first two points

srcinfo = []; srcinfo.r = trap.r(:,1); srcinfo.d = trap.d(:,1);
srcinfo.d2 = trap.d2(:,1);
targinfo = []; targinfo.r = trap.r(:,2); targinfo.d = trap.d(:,2);
ftemp = kern(srcinfo,targinfo);
opdims = size(ftemp);

% get opts from struct if available

if isfield(opts,'quadorder')
    quadorder = opts.quadorder;
end
if isfield(opts,'quad')
    quad = opts.quad;
end
if isfield(opts,'nonsmoothonly')
    nonsmoothonly = opts.nonsmoothonly;
end

% call requested routine

if or(strcmpi(quad,'balog'),strcmpi(quad,'smooth'))

    type = 'log';
    if strcmpi(quad,'smooth')
        type = 'smooth';
    end
    if nonsmoothonly
        sysmat = chnk.quadba.buildmattd(trap,kern,quadorder,opdims,type);
    else
        sysmat = chnk.quadba.buildmat(trap,kern,quadorder,opdims,type);
    end
    
else
    warning('specified quadrature method "%s" not available',quad);
    sysmat = [];
end
	 

end
