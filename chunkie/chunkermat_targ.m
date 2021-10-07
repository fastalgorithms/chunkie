function mat = chunkermat_targ(chnkrsrc,chnkrtarg,srcid,targid,fkern,opdims,type)

if strcmpi(type,'smooth')
  targinfo = []; targinfo.r = chnkrtarg.r(:,:,targid); 
  targinfo.d = chnkrtarg.d(:,:,targid); targinfo.d2 = chnkrtarg.d2(:,:,targid);
  nt = chnkrtarg.k*length(targid);

  srcinfo = []; srcinfo.r = chnkrsrc.r(:,:,srcid); 
  srcinfo.d = chnkrsrc.d(:,:,srcid); srcinfo.d2 = chnkrsrc.d2(:,:,srcid);

  mat = fkern(srcinfo,targinfo);

  wts = weights(chnkrsrc);
  wts = wts(:,srcid);
  wts2 = repmat( (wts(:)).', nt, 1);

  mat = mat.*wts2;
end

if strcmpi(type,'nearlog')
  k = chnkrsrc.k;
  
  qavail = chnk.quadggq.logavail();
  [~,i] = min(abs(qavail-k));
  assert(qavail(i) == k,'order %d not found, consider using order %d chunks', ...
      k,qavail(i));
  [xs1,wts1,xs0,~] = chnk.quadggq.getlogquad(k);

  ainterp1 = lege.matrin(k,xs1);
  temp = eye(opdims(2));
  ainterp1kron = kron(ainterp1,temp);

  nquad0 = size(xs0,1);

  ainterps0kron = zeros(opdims(2)*nquad0,opdims(2)*k,k);
  ainterps0 = zeros(nquad0,k,k);

  for j = 1:k
    xs0j = xs0(:,j);
    ainterp0_sm = lege.matrin(k,xs0j);
    ainterps0(:,:,j) = ainterp0_sm;
    ainterps0kron(:,:,j) = kron(ainterp0_sm,temp);
  end
  
  rs = chnkrsrc.r(:,:,srcid); ds = chnkrsrc.d(:,:,srcid); d2s = chnkrsrc.d2(:,:,srcid); 
  rt = chnkrtarg.r(:,:,targid); dt = chnkrtarg.d(:,:,targid); d2t = chnkrtarg.d2(:,:,targid); 
  hs = chnkrsrc.h(srcid);
  
  rfine = (ainterp1*(rs.')).'; dfine = (ainterp1*(ds.')).'; 
  d2fine = (ainterp1*(d2s.')).';

  srcinfo = []; srcinfo.r = rfine; srcinfo.d = dfine; 
  srcinfo.d2 = d2fine;

  targinfo = []; targinfo.r = rt; targinfo.d = dt; 
  targinfo.d2 = d2t;

  dfinenrm = sqrt(sum(dfine.^2,1));
  %dfinenrm = dfine(1,:,:); % for complex contour, by SJ 09/30/21
  dsdt = dfinenrm(:).*wts1(:)*hs;

  dsdtndim2 = repmat(dsdt(:).',opdims(2),1);
  dsdtndim2 = dsdtndim2(:);

  % get kernel values and then premultiply by interpolating matrix

  smatbig = fkern(srcinfo,targinfo);
  mat = smatbig*diag(dsdtndim2)*ainterp1kron;
end

end