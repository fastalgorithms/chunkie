function npxy = nproxy_square(kern, width, opts)
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

    if isa(kern,'kernel')
        kern = kern.eval;
    end
  if nargin < 3
    opts = [];
  end

  nsrc = 200;           if isfield(opts,'nsrc'), nsrc = opts.nsrc; end
  rank_or_tol = 1e-13;  if isfield(opts,'rank_or_tol'), rank_or_tol = opts.rank_or_tol; end

  if isfield(opts,'eps') % kept for backward compat
    rank_or_tol = opts.eps;
  end
	
  % set up sources with randomly oriented d and d2

  srcinfo = [];
  srcinfo.r = [-0.5;-0.5]*width + rand(2,nsrc)*width;
  
  srcinfo.d = randn(2,nsrc);
  srcinfo.d2 = randn(2,nsrc);
  srcinfo.n = [-srcinfo.d(2,:);srcinfo.d(1,:)] ./ sqrt(sum(srcinfo.d.^2,1));

  stmp = [];
  stmp.r = randn(2,1); stmp.d = randn(2,1); stmp.d2 = randn(2,1);
  stmp.n = [-stmp.d(2,:);stmp.d(1,:)] ./ sqrt(sum(stmp.d.^2,1));

  ttmp = [];
  ttmp.r = randn(2,1); ttmp.d = randn(2,1); ttmp.d2 = randn(2,1);
  ttmp.n = [-ttmp.d(2,:);ttmp.d(1,:)] ./ sqrt(sum(ttmp.d.^2,1));
  
  f = kern(stmp, ttmp);
  [opdims(1), opdims(2)] = size(f);

  % strengths of random sources
  sig = randn(opdims(2)*nsrc,1);

  npxy = 64;

  % double until you get enough
  last_intgrl = NaN;
  one_more    = true;
  for i = 0:14
    [pr,ptau,pw] = chnk.flam.proxy_square_pts(npxy);
    pwuse =  repmat(pw.', [opdims(1),1]); pwuse = pwuse(:);
    targinfo.r = width*pr;
    targinfo.d = ptau;
    targinfo.d2 = zeros(2,npxy);
    targinfo.n = [-ptau(2,:);ptau(1,:)] ./ sqrt(sum(ptau.^2,1));

    % compute integral along proxy surface of potential from random sources
    % to see if we've resolved the kernel
    intgrl = pwuse' * kern(srcinfo,targinfo) * sig;
   
    % check self convergence
    err = abs(intgrl - last_intgrl) / abs(intgrl);
    if err < rank_or_tol || ~one_more
      return;
    end
    npxy = 2*npxy;
    last_intgrl = intgrl;

    % if using very low tolerance like 1e-14, just get to 1e-12 error
    % and double one more time to avoid mucking around at machine eps
    if rank_or_tol < 1e-12 && err < 1e-12
        one_more = false;
    end
  end

  % kernel was not resolved even with 2^14 points, report failure
  npxy = -1;

end