function submat = diagbuildmat(r,d,n,d2,data,i,fkern,opdims,...
			 xs0,whts0,ainterps0kron,ainterps0,corrections,wtss,indd)
%CHNK.QUADGGQ.DIAGBUILDMAT                  

% grab specific boundary data
                
rs = r(:,:,i); ds = d(:,:,i); d2s = d2(:,:,i);
ns = n(:,:,i);
if(isempty(data))
    dd = [];
else
    dd = data(:,:,i);
end

if nargin < 13
    corrections = false;
    wtss = [];
    indd = [];
end


% interpolate boundary info

% get relevant coefficients

[~,k] = size(rs);            
            
rfine = cell(k,1);
dfine = cell(k,1);
nfine = cell(k,1);
d2fine = cell(k,1);
dsdt = cell(k,1);

for j = 1:k
    rfine{j} = (ainterps0{j}*(rs.')).';
    dfine{j} = (ainterps0{j}*(ds.')).';
    d2fine{j} = (ainterps0{j}*(d2s.')).';
    dfinenrm = sqrt(sum(dfine{j}.^2,1));
    nfine{j} = [dfine{j}(2,:); -dfine{j}(1,:)]./dfinenrm;
    dsdt{j} = (dfinenrm(:)).*whts0{j};
    if(not (isempty(data)))
	ddfine{j} = ((ainterps0{j}*(dd.'))).';
    end
end

srcinfo = [];
targinfo = [];

for j = 1:k
  srcinfo.r = rfine{j}; srcinfo.d = dfine{j};
  srcinfo.d2 = d2fine{j}; srcinfo.n = nfine{j};
  targinfo.r = rs(:,j); targinfo.d = ds(:,j);
  targinfo.d2 = d2s(:,j); targinfo.n = ns(:,j);
  if(isempty(dd))
      targinfo.data = [];
      srcinfo.data = [];
  else
    srcinfo.data = ddfine{j};
    targinfo.data = dd(:,j);
  end
    
  smatbigi = fkern(srcinfo,targinfo);
  dsdtndim2 = repmat(dsdt{j}.',opdims(2),1);
  dsdtndim2 = dsdtndim2(:);
  smatbigi = bsxfun(@times,smatbigi,dsdtndim2.');
  submat(opdims(1)*(j-1)+1:opdims(1)*j,:) = ...
    smatbigi*ainterps0kron{j};
end

if corrections
    srcinfo.r = rs;
    srcinfo.d = ds;
    srcinfo.d2 = d2s;
    srcinfo.n = ns;
    srcinfo.data = dd;

    targinfo.r = rs; targinfo.d = ds; targinfo.d2 = d2s; targinfo.n = ns;
    targinfo.data = dd;

    sc = fkern(srcinfo,targinfo);
    sc(indd) = 0;
    
    wtsi = wtss(:,i);
    submat = submat - sc.*(wtsi(:).');

end
   

