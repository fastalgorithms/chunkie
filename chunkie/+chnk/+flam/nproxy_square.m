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
%       opts.rank_or_tol - tolerance or maximum rank (default 1e-13)
%       opts.eps - alias for rank_or_tol for backward compat (default 1e-13)
%       opts.nsrc - number of sources to use in test (default 200)  
%
% Output:
%   npxy - number of points on perimeter of box to be sent to proxy routine
%

  nsrc = 200;
  rank_or_tol = 1e-13;

  if nargin < 3
    opts = [];
  end

  if isfield(opts,'eps')
    rank_or_tol = opts.eps;
  end
  if isfield(opts,'rank_or_tol')
    rank_or_tol = opts.rank_or_tol;
  end
  if isfield(opts,'nsrc')
    nsrc = opts.nsrc;
  end

		      % set up sources with randomly oriented d and d2

  srcinfo = [];
  srcinfo.r = [-0.5;-0.5]*width + rand(2,nsrc)*width;
  
  srcinfo.d = randn(2,nsrc);
  theta = randn(1,nsrc);
  srcinfo.n = [cos(theta); sin(theta)];
  srcinfo.d2 = randn(2,nsrc);

  npxy = 16;

				% double until you get enough
  for i = 1:11
    [pr,ptau] = chnk.flam.proxy_square_pts(npxy);
    targinfo.r = pr*width;
    targinfo.d = ptau;
    targinfo.d2 = zeros(2,npxy);

    mat = kern(srcinfo,targinfo);

    [sk,~] = id(mat,rank_or_tol);
    
    if length(sk) < min(nsrc,npxy)
      npxy = (floor((length(sk)-1)/4+0.1)+1)*4;
      break;
    end
    npxy = 2*npxy;
  end
  
end
