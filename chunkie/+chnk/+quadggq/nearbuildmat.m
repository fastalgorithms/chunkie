function submat = chunksnearbuildmat(r,d,h,i,j,fkern,opdims,...
				      u,xs1,whts1,ainterp1)
%CHUNKSNEARBUILDMAT

% grab specific boundary data
                
rs = r(:,:,j); rt = r(:,:,i); ds = d(:,:,j); dt = d(:,:,i);
hs = h(j); ht = h(i);

% interpolate boundary info

% get relevant coefficients

[dim,~] = size(rs);            



            
rsc = u*rs.';
dsc = u*ds.';

% then interpolate 

% rfine = zeros(length(xs1),dim);
% dfine = zeros(length(xs1),dim);\



xs1 = xs1(:);

rfine = lege.exev(xs1,rsc);
dfine = lege.exev(xs1,dsc);

% for i = 1:dim
%     rfine(:,i) = lege.exev(xs1,rsc(:,i));
%     dfine(:,i) = lege.exev(xs1,dsc(:,i));
% end



rfine = rfine.';

dfinenrms = sqrt(sum(abs(dfine).^2,2));
taufine = (bsxfun(@rdivide,dfine,dfinenrms)).';

dtnrms = sqrt(sum(abs(dt).^2,2));
taut = (bsxfun(@rdivide,dt,dtnrms)).';
    
dsdt = dfinenrms.*whts1*hs;

dsdtndim2 = repmat(dsdt(:).',opdims(2),1);
dsdtndim2 = dsdtndim2(:);

% get kernel values and then premultiply by interpolating matrix

smatbig = fkern(rfine,rt,taufine,taut);
submat = smatbig*diag(dsdtndim2)*ainterp1;

 end
   
