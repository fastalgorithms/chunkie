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
for i=1:ntarg
  zt = targs(1,i)+1i*targs(2,i);
  for j=1:ndomain
    if(clmparams.is_inf(j) == 0)
      n = 0;
      for k=1:length(clmparams.clist{j})
        kcurve = abs(clmparams.clist{j}(k));
        tmp = sum(wts{kcurve}./(zsrc{kcurve}-zt));
        if clmparams.clist{j}(k) > 0
          n = n + tmp;
        else
          n = n - tmp;
        end
      end
      n = round(n/(2*pi*1i));
      if n==1
        targdomain(i) = j;
        break;
      end
    end
  end
  
  if targs(2,i)>0 && targdomain(i)==0
    targdomain(i) = idomup;
  elseif targs(2,i)<0 && targdomain(i)==0
    targdomain(i) = idomdown;
  end
end

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