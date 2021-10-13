function targdomain = finddomain(chnkr,clist,targs)
% use Cauchy's integral formula to determine which domain each target lies
% in.
ncurve = length(chnkr);
k = chnkr(1).k;
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
    
for i=1:ntarg
  zt = targs(1,i)+1i*targs(2,i);
  for j=3:5
    n = 0;
    for k=1:length(clist{j})
      kcurve = abs(clist{j}(k));
      tmp = sum(wts{kcurve}./(zsrc{kcurve}-zt));
      if clist{j}(k) > 0
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
  if targs(2,i)>0 && targdomain(i)==0
    targdomain(i) = 1;
  elseif targs(2,i)<0 && targdomain(i)==0
    targdomain(i) = 2;
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