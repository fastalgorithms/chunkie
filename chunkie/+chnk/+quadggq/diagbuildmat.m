function submat = diagbuildmat(r,d,n,d2,h,data,i,fkern,opdims,...
			 xs0,whts0,ainterps0kron,ainterps0)
%CHNK.QUADGGQ.DIAGBUILDMAT                  

% grab specific boundary data
                
rs = r(:,:,i); ds = d(:,:,i); d2s = d2(:,:,i); hs = h(i); 
ns = n(:,:,i);
if(isempty(data))
    dd = [];
else
    dd = data(:,:,i);
end
% interpolate boundary info

% get relevant coefficients

[~,k] = size(rs);            
            
rfine = cell(k,1);
dfine = cell(k,1);
nfine = cell(k,1);
d2fine = cell(k,1);
dsdt = cell(k,1);

for i = 1:k
    rfine{i} = (ainterps0{i}*(rs.')).';
    dfine{i} = (ainterps0{i}*(ds.')).';
    d2fine{i} = (ainterps0{i}*(d2s.')).';
    dfinenrm = sqrt(sum(dfine{i}.^2,1));
    nfine{i} = [dfine{i}(2,:); -dfine{i}(1,:)]./dfinenrm;
    dsdt{i} = (dfinenrm(:)).*whts0{i}*hs;
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
  else
    targinfo.data = dd(:,j);
  end
    
  smatbigi = fkern(srcinfo,targinfo);
  dsdtndim2 = repmat(dsdt{j}.',opdims(2),1);
  dsdtndim2 = dsdtndim2(:);
  submat(opdims(1)*(j-1)+1:opdims(1)*j,:) = ...
    smatbigi*diag(dsdtndim2)*ainterps0kron{j};
end

end
   

