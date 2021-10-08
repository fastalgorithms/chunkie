function submat = nearbuildmat(r,d,d2,h,i,j,fkern,opdims,...
				      xs1,whts1,ainterp1kron,ainterp1)
%CHNKR.QUADGGQ.NEARBUILDMAT

% grab specific boundary data
                
rs = r(:,:,j); ds = d(:,:,j); d2s = d2(:,:,j); 
rt = r(:,:,i); dt = d(:,:,i); d2t = d2(:,:,i); 
hs = h(j);

% interpolate boundary info

% get relevant coefficients

% [dim,~] = size(rs);            
% 
% 
% 
%             
% rsc = u*rs.';
% dsc = u*ds.';
% 
% % then interpolate 
% 
% % rfine = zeros(length(xs1),dim);
% % dfine = zeros(length(xs1),dim);\
% 
% 
% 
% xs1 = xs1(:);
% 
% rfine = lege.exev(xs1,rsc);
% dfine = lege.exev(xs1,dsc);
% 
% % for i = 1:dim
% %     rfine(:,i) = lege.exev(xs1,rsc(:,i));
% %     dfine(:,i) = lege.exev(xs1,dsc(:,i));
% % end
% 
% 
% 
% rfine = rfine.';
% 
% dfinenrms = sqrt(sum(abs(dfine).^2,2));
% taufine = (bsxfun(@rdivide,dfine,dfinenrms)).';
% 
% dtnrms = sqrt(sum(abs(dt).^2,1));
% taut = (bsxfun(@rdivide,dt,dtnrms));
%     
% dsdt = dfinenrms.*whts1*hs;

rfine = (ainterp1*(rs.')).'; dfine = (ainterp1*(ds.')).'; 
d2fine = (ainterp1*(d2s.')).';

srcinfo = []; srcinfo.r = rfine; srcinfo.d = dfine; 
srcinfo.d2 = d2fine;

targinfo = []; targinfo.r = rt; targinfo.d = dt; 
targinfo.d2 = d2t;

dfinenrm = sqrt(sum(dfine.^2,1));
%dfinenrm = dfine(1,:,:); % for complex contour, by SJ 09/30/21
dsdt = dfinenrm(:).*whts1(:)*hs;

dsdtndim2 = repmat(dsdt(:).',opdims(2),1);
dsdtndim2 = dsdtndim2(:);

% get kernel values and then premultiply by interpolating matrix

smatbig = fkern(srcinfo,targinfo);
submat = smatbig*diag(dsdtndim2)*ainterp1kron;

 end
   
