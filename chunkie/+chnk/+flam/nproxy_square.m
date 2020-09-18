function npxy = nproxy_square(kern,width,opts)
%CHNK.FLAM.NPROXY_SQUARE find the number of points to place on a
% side of a proxy surface for given kernel, width, and tolerance
%
% Syntax: npxy = chnk.flam.nproxy_square(kern,width,eps)
%
% Input:
%   kern - kernel function of the form kern(srcinfo,targinfo)
%   width - box/cube width (of center box, proxy surface at
%                1.5*width)
%   opts - options structure
%       opts.eps - tolerance (default 1e-13)
%       opts.nsrc - number of sources to use in test (default 200)  
%
% Output:
%   npxy - number of points on side of box to be passed
%       to proxy surface routine
%

  nsrc = 200;
  eps = 1e-13;

  if nargin < 3
    opts = [];
  end

  if isfield(opts,'eps')
    eps = opts.eps;
  end
  if isfield(opts,'nsrc')
    nsrc = opts.nsrc;
  end

		      % set up sources with randomly oriented d and d2

  srcinfo = [];
  srcinfo.r = [-0.5;-0.5]*width + rand(2,nsrc)*width;
  
  srcinfo.d = randn(2,nsrc);
  srcinfo.d2 = randn(2,nsrc);

  npxy = 16;

				% double until you get enough
  for i = 1:11
    [pr,ptau] = chnk.flam.proxy_square_pts(npxy);
    targinfo.r = pr*width;
    targinfo.d = ptau;
    targinfo.d2 = zeros(2,npxy);

    mat = kern(srcinfo,targinfo);

    [sk,~] = id(mat,eps);
    
    if length(sk) < min(nsrc,npxy)
      npxy = (floor((length(sk)-1)/4+0.1)+1)*4;
      break;
    end
    npxy = 2*npxy;
  end
  
end
