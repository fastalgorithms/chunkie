function submat = chunksdiagbuildmat(r,d,h,i,fkern,opdims,...
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


rspec = lege.exev(xs0,rsc);
dspec = lege.exev(xs0,dsc);
% 
% rspec = zeros(size(xs0,1),k,dim);
% dspec = zeros(size(xs0,1),k,dim);
%    
% for j = 1:dim
%     rspec(:,:,j) = lege.exev(xs0,rsc(:,j));
%     dspec(:,:,j) = lege.exev(xs0,dsc(:,j));
% end

dspecnrms = sqrt(sum(abs(dspec).^2,3));
tauspec = bsxfun(@rdivide,dspec,dspecnrms);

dsnrms = sqrt(sum(abs(ds).^2,1));
taus = bsxfun(@rdivide,ds,dsnrms);

dsdt = dspecnrms.*whts0*hs;

% get kernel values and then premultiply by interpolating matrix
% special nodes are the sources and the targs are the regular points

submat = zeros(opdims(1)*k,opdims(2)*k);

rspec = permute(rspec,[3 1 2]);
tauspec = permute(tauspec,[3,1,2]);

for j = 1:k
    smatbigi = fkern(rspec(:,:,j),rs(:,j),tauspec(:,:,j),taus(:,j));
    dsdtndim2 = repmat(dsdt(:,j).',opdims(2),1);
    dsdtndim2 = dsdtndim2(:);
    submat(opdims(1)*(j-1)+1:opdims(1)*j,:) = ...
        smatbigi*diag(dsdtndim2)*ainterps0(:,:,j);
end

 end
   

