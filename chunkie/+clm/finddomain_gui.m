function targdomain = finddomain_gui(chnkr,clmparams,targs)
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
      srcinfo.dipstr = complex(zeros(ns,1));
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

      pot = cfmm2d(eps0,srcinfo,targs);
      pot = round(pot/(2*pi*1i));
      targdomain(pot==1) = j;
    end
end
targdomain(targdomain(:) == 0 & (targs(2,:)>0)') = 1;
targdomain(targdomain == 0 & (targs(2,:)<0)') = 2;


for i=1:ntarg
  zt = targs(1,i)+1i*targs(2,i);
  d0 = 1000;
  for j=1:ncurve
    d = min(abs(zt-zsrc{j}));
    if d<d0
      d0=d;
    end
  end
  if d0<0.02
    targdomain(i)=0;
  end
end