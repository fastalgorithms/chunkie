function submat = diagbuildmat(r,d,d2,h,i,fkern,opdims,...
				      xs0,whts0,ainterps0kron,ainterps0)
%CHNK.QUADGGQ.DIAGBUILDMAT                  

% grab specific boundary data
                
rs = r(:,:,i); ds = d(:,:,i); d2s = d2(:,:,i); hs = h(i); 

% interpolate boundary info

% get relevant coefficients

[dim,k] = size(rs);            
            
[nq,~,nt] = size(ainterps0);
rfine = zeros(dim,nq,nt);
dfine = zeros(dim,nq,nt);
d2fine = zeros(dim,nq,nt);
for i = 1:nt
    rfine(:,:,i) = (ainterps0(:,:,i)*(rs.')).';
    dfine(:,:,i) = (ainterps0(:,:,i)*(ds.')).';
    d2fine(:,:,i) = (ainterps0(:,:,i)*(d2s.')).';    
end

%  rsc = u*(rs.');
%  dsc = u*(ds.');
% d2sc = u*(d2s.'); 
%  
% % % then interpolate 
%  
%  
%  rfine = lege.exev(xs0,rsc);
%  dfine = lege.exev(xs0,dsc);
%  d2fine = lege.exev(xs0,d2sc);
% % 
% 
% %  rfine = zeros(size(xs0,1),k,dim);
% %  dfine = zeros(size(xs0,1),k,dim);
% % d2fine = zeros(size(xs0,1),k,dim); 
% % %    
% %  for j = 1:dim
% %      rfine(:,:,j) = lege.exev(xs0,rsc(:,j));
% %      dfine(:,:,j) = lege.exev(xs0,dsc(:,j));
% %      d2fine(:,:,j) = lege.exev(xs0,d2sc(:,j));
% %  end
% 
% %dspecnrms = sqrt(sum(abs(dfine).^2,3));
% %tauspec = bsxfun(@rdivide,dspec,dspecnrms);
% 
% %dsnrms = sqrt(sum(abs(ds).^2,1));
% %taus = bsxfun(@rdivide,ds,dsnrms);
% 
% %dsdt = dspecnrms.*whts0*hs;
% 
% % get kernel values and then premultiply by interpolating matrix
% % special nodes are the sources and the targs are the regular points
% 
% submat = zeros(opdims(1)*k,opdims(2)*k);
% 
% %rspec = permute(rspec,[3 1 2]);
% %tauspec = permute(tauspec,[3,1,2]);
% rfine = permute(rfine,[3,1,2]);
% dfine = permute(dfine,[3,1,2]);
% d2fine = permute(d2fine,[3,1,2]);

dfinenrm = squeeze(sqrt(sum(dfine.^2,1)));
%dfinenrm = squeeze(dfine(1,:,:)); % for complex contour, added by SJ 9/30/21
dsdt = dfinenrm.*whts0*hs;

srcinfo = [];
targinfo = [];

for j = 1:k
    srcinfo.r = rfine(:,:,j); srcinfo.d = dfine(:,:,j); 
    srcinfo.d2 = d2fine(:,:,j);
    targinfo.r = rs(:,j); targinfo.d = ds(:,j); targinfo.d2 = d2s(:,j);
    smatbigi = fkern(srcinfo,targinfo);
    dsdtndim2 = repmat(dsdt(:,j).',opdims(2),1);
    dsdtndim2 = dsdtndim2(:);
    submat(opdims(1)*(j-1)+1:opdims(1)*j,:) = ...
        smatbigi*diag(dsdtndim2)*ainterps0kron(:,:,j);
end

end
   

