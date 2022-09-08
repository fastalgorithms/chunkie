function [targdomain,varargout] = finddomain_gui(chnkr,clmparams,targs)
% use Cauchy's integral formula to determine which domain each target lies
% in.
ncurve = clmparams.ncurve;
ndomain = clmparams.ndomain;
k = clmparams.ngl;
[~,w] = lege.exps(k);

wts = cell(1,ncurve);
for i=1:ncurve
  zp = chnkr(i).d(1,:)+1i*chnkr(i).d(2,:);
  hs = chnkr(i).h;
  ws = kron(hs(:),w(:));
  wts{i} = zp(:).*ws(:);
end

zsrc = cell(1,ncurve);
for i=1:ncurve
  zs = chnkr(i).r(1,:) + 1i*chnkr(i).r(2,:);
  zsrc{i} = zs(:);
end

ntarg = size(targs,2);
targdomain = zeros(ntarg,1);
targnear = zeros(ntarg,1);
    
idomup = find(clmparams.is_inf == 1);
idomdown = find(clmparams.is_inf == -1);
eps0 = 1e-3;
for j=1:ndomain
    if(clmparams.is_inf(j) == 0)
      srcinfo = [];
      ns = 0;
      for k=1:length(clmparams.clist{j})
          kcurve = abs(clmparams.clist{j}(k));
          ns = ns + length(zsrc{kcurve});
      end
      srcinfo.sources = zeros(2,ns);
      srcinfo.dipstr = complex(zeros(1,ns));
      istart = 1;
      for k=1:length(clmparams.clist{j})


        kcurve = abs(clmparams.clist{j}(k));
        nc = length(zsrc{kcurve});
        iend = istart + nc-1;

        srcinfo.sources(1,istart:iend) = real(zsrc{kcurve});
        srcinfo.sources(2,istart:iend) = imag(zsrc{kcurve});
        if clmparams.clist{j}(k) > 0
          srcinfo.dipstr(istart:iend) = -wts{kcurve};
        else
          srcinfo.dipstr(istart:iend) = wts{kcurve};
        end
        istart = istart + nc;
      end
      pg = 0;
      pgt = 1;
      U = cfmm2d(eps0,srcinfo,pg,targs,pgt);
      pot = round(real(U.pottarg/(2*pi*1i)));
      targdomain(pot==1) = j;
    end
end
targdomain(targdomain(:) == 0 & (targs(2,:)>0)') = idomup;
targdomain(targdomain == 0 & (targs(2,:)<0)') = idomdown;

chnkrtotal = merge(chnkr);
flag = flagnear(chnkrtotal,targs);

[tind,sind,~] = find(flag);
targs_test = targs(:,tind);
src_test = chnkrtotal.r(:,:,sind);
[~,m,n] = size(src_test);
xt = targs_test(1,:);
yt = targs_test(2,:);
xt = repmat(xt,[chnkrtotal.k,1]);
xt = xt(:).';

yt = repmat(yt,[chnkrtotal.k,1]);
yt = yt(:).';
src_test2 = reshape(src_test,[2,m*n]);
dist = (xt - src_test2(1,:)).^2 + (yt-src_test2(2,:)).^2;
dist = reshape(dist,[m,n]);
dmin = min(dist,[],1);
iind = dmin < 0.02^2;
targdomain(tind(iind))=0;
varargout{1} = tind;
varargout{2} = flag;
