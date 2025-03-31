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
%       opts.npxy_method - 'wavenumber'
%                         use kernel class zk parameter to determine
%                         number of terms, overwritten if
%                         opts.npxy_wavenumber (zk) is provided, and set to 
%                         0 if none of the parameters in kernel class
%                         have zk as a parameter
%                     'id' use an id with increasing number of sources
%                          to determine npxy
%                     'integral' compute a integral of kernel with a 
%                                random function to determine 
%                                the number of proxy points
%       opts.npxy_wavenumber - wavenumber for determining number of proxy
%                              points
%
% Output:
%   npxy - number of points on perimeter of box to be sent to proxy routine
%

    
  if nargin < 3
    opts = [];
  end

  nsrc = 200;           if isfield(opts,'nsrc'), nsrc = opts.nsrc; end
  rank_or_tol = 1e-13;  if isfield(opts,'rank_or_tol'), rank_or_tol = opts.rank_or_tol; end
  method = 'wavenumber';  if isfield(opts,'npxy_method'), method = opts.npxy_method; end


  eps = rem(rank_or_tol, 1);

  zkuse = 0;
  msg = "CHNK.FLAM.NPROXY_SQUARE: unable to determine wavenumber from kernel, setting it to 0";
  if isfield(opts, 'npxy_wavenumber') 
    zkuse = opts.npxy_wavenumber;
  else
    if isa(kern, 'kernel') 
      ptmp = kern.params;
      if isa(ptmp, 'struct')
          bb = vertcat(ptmp(:));
          if isfield(bb, 'zk')
            zk = [bb.zk];
            zkuse = max(abs(real(zk(:))));
          end
          if isfield(bb, 'zks')
            zks = [bb.zks];
            zkuse = max([zkuse; abs(real(zks(:)))]);
          end
      elseif isa(ptmp, 'cell')
          [m, n] = size(ptmp);
          zkuse = 0;
          for j=1:n
              for i=1:m
                  a = ptmp{i,j};
                  if isa(a, 'struct')
                      if isfield(a, 'zk')
                          zkuse = max(zkuse, abs(real(a.zk)));
                      end
                      if isfield(a, 'zks')
                          zkuse = max([zkuse; abs(real(a.zks(:)))]);
                      end
                  end
              end
          end
      else
          warning(msg)
      end
    else
      warning(msg);
    end
  end

  if isa(kern, 'kernel'), kern = kern.eval; end

  if isfield(opts,'eps') % kept for backward compat
    eps = opts.eps;
  end
	
  % set up sources with randomly oriented d and d2

  srcinfo = [];
  srcinfo.r = [-0.5;-0.5]*width + rand(2,nsrc)*width;
  
  srcinfo.d = randn(2,nsrc);
  srcinfo.d2 = randn(2,nsrc);
  srcinfo.n = [-srcinfo.d(2,:);srcinfo.d(1,:)] ./ sqrt(sum(srcinfo.d.^2,1));

  % determine opdims of kernel
  stmp = [];
  stmp.r = randn(2,1); stmp.d = randn(2,1); stmp.d2 = randn(2,1);
  stmp.n = [-stmp.d(2,:);stmp.d(1,:)] ./ sqrt(sum(stmp.d.^2,1));

  ttmp = [];
  ttmp.r = randn(2,1); ttmp.d = randn(2,1); ttmp.d2 = randn(2,1);
  ttmp.n = [-ttmp.d(2,:);ttmp.d(1,:)] ./ sqrt(sum(ttmp.d.^2,1));
  
  f = proxy_square_mat(kern, stmp, ttmp);
  [opdims(1), opdims(2)] = size(f);

  % strengths of random sources
  sig = randn(opdims(2)*nsrc,1);

  npxy = 16;

  if strcmpi(method, 'wavenumber')
    nmax = 4*zkuse*width + 50;
    r1 = width/sqrt(2);
    r2 = width*1.5;
    if abs(zkuse*width) > 1e-6
      js = besselj(1:nmax, zkuse*r1);
      hs = besselh(1:nmax, zkuse*r2);
      
    
    else
      js = r1.^(1:nmax);
      hs = 1./r2.^(1:nmax);
      
    end
    rfacs = abs(hs).*abs(js);
    nn = find(rfacs < eps, 3, 'first');

    nmin = max(ceil(3*width*zkuse/2/pi*4),5);
    nn = max([nmin, nn, ceil(rank_or_tol)]);
    npxy = floor(nn/4+1.1)*4;
    return

  else

    % double until you get enough
    last_intgrl = NaN;
    one_more    = true;
    for i = 0:17
      [pr,ptau,pw] = chnk.flam.proxy_square_pts(npxy);
      pwuse =  repmat(pw.', [opdims(1),1]); pwuse = pwuse(:);
      targinfo.r = width*pr;
      targinfo.d = ptau;
      targinfo.d2 = zeros(2,npxy);
      targinfo.n = [-ptau(2,:);ptau(1,:)] ./ sqrt(sum(ptau.^2,1));

      if strcmpi(method, 'integral')

        % compute integral along proxy surface of potential from random sources
        % to see if we've resolved the kernel
        intgrl = pwuse' * proxy_square_mat(kern, srcinfo, targinfo) * sig;
   
        % check self convergence
        err = abs(intgrl - last_intgrl) / sum(abs(sig(:)));
        if err < eps || ~one_more
          npxy = max(npxy/2, floor(rank_or_tol/4+1.1)*4);
          return;
        end
        npxy = 2*npxy;
        last_intgrl = intgrl;

        % if using very low tolerance like 1e-14, just get to 1e-12 error
        % and double one more time to avoid mucking around at machine eps
        if eps < 1e-12 && err < 1e-12
          one_more = false;
        end
      elseif strcmpi(method, 'id')
         
        if npxy > nsrc
          nsrc = ceil(1.2*npxy);
  
          srcinfo = [];
          srcinfo.r = [-0.5;-0.5]*width + rand(2,nsrc)*width;
  
          srcinfo.d = randn(2,nsrc);
          srcinfo.d2 = randn(2,nsrc);
          srcinfo.n = [-srcinfo.d(2,:);srcinfo.d(1,:)] ./ sqrt(sum(srcinfo.d.^2,1));
        end
      
        mat = proxy_square_mat(kern, srcinfo, targinfo);
        [sk, ~] = id(mat, rank_or_tol);

        if length(sk) < npxy
          npxy = floor((length(sk)/opdims(2)-1)/4+1.1)*4;
          return
        end
        npxy = 2*npxy;
      else
        msg = "NPROXY_SQUARE: Invalid method type for determining number of proxy points";
        error(msg);
      end
    end
  end

  % kernel was not resolved even with 2^14 points, report failure
  npxy = -1;

end


function mat = proxy_square_mat(kern,srcinfo,targinfo)

  if numel(kern) == 1
    mat = kern(srcinfo,targinfo);
    return;
  end

  % kern is a cellmat

  [m,n] = size(kern);

  mat = cell(m,n);

  for i=1:m
    for j=1:n
      mat{i,j} = kern{i,j}(srcinfo,targinfo);
    end
  end

  mat = cell2mat(mat);
end
