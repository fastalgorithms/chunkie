function submat = chunksdiagbuildmat(r,d,h,i,fkern,ndims,...
				      u,xs0,whts0,ainterps0)
%CHUNKSDIAGBUILDMAT                  

% grab specific boundary data
                
rs = r(:,:,i); ds = d(:,:,i); hs = h(i); 

% interpolate boundary info

% get relevant coefficients

[dim,k] = size(rs);            
            
rsc = u*(rs.');
dsc = u*(ds.');

% then interpolate 

rspec = zeros(size(xs0,1),k,dim);
dspec = zeros(size(xs0,1),k,dim);
   
for j = 1:dim
    rspec(:,:,j) = lege.exev(xs0,rsc(:,j));
    dspec(:,:,j) = lege.exev(xs0,dsc(:,j));
end

dspecnrms = sqrt(sum(abs(dspec).^2,3));
tauspec = bsxfun(@rdivide,dspec,dspecnrms);

dsnrms = sqrt(sum(abs(ds).^2,1));
taus = bsxfun(@rdivide,ds,dsnrms);

dsdt = dspecnrms.*whts0*hs;

% get kernel values and then premultiply by interpolating matrix
% special nodes are the sources and the targs are the regular points

submat = zeros(ndims(1)*k,ndims(2)*k);
for j = 1:k
    rspecj = (squeeze(rspec(:,j,:))).';
    tauspecj = (squeeze(tauspec(:,j,:))).';
    smatbigi = fkern(rspecj,rs(:,j),tauspecj,taus(:,j));
    dsdtndim2 = repmat(dsdt(:,j).',ndims(2),1);
    dsdtndim2 = dsdtndim2(:);
    submat(ndims(1)*(j-1)+1:ndims(1)*j,:) = ...
        smatbigi*diag(dsdtndim2)*ainterps0(:,:,j);
end

 end
   

