function [val, grad, hess] = green_0th_mode(src, targ, origin)
%
% CHNK.AXISSYMHELM2D.GREEN evaluate the Laplace green's function
% for the given sources and targets. 
%
% Note: that the first coordinate is r, and the second z.
% The code relies on precomputed tables and hence loops are required for 
% computing various pairwise interactions.
% Finally, the code is not efficient in the sense that val, grad, hess 
% are always internally computed independent of nargout
%
% Returns for gradient are:
% grad = d_{r}, d_{r'}, d_{z}, d_{z'}
%
% Returns for hess are:
% hess = d_{rr'}, d_{zz'}, d_{rz'}, d_{r'z}

[~, ns] = size(src);
[~, nt] = size(targ);

rt = repmat(targ(1,:).',1,ns); % r
rs = repmat(src(1,:),nt,1); % r'
r  = (rt + origin(1));
rp = (rs + origin(1));
dr = rt-rs; % r - r'
z  = repmat(targ(2,:).',1,ns);
zp = repmat(src(2,:),nt,1);
dz = z-zp; % z - z'

[gs,gdzs,gdrs,gdrps,gdrpr,gdzz,gdrz,gdrpz] = chnk.axissymlap2d.gfunc(r,rp,dr,z,zp,dz);

const = 1/(4*pi^2);

val = gs*const;
grad = zeros(nt, ns, 4);
grad(:,:,1) = gdrs;
grad(:,:,2) = gdrps;
grad(:,:,3) = gdzs;
grad(:,:,4) = -gdzs;
grad = grad*const;

hess = zeros(nt, ns, 4);
hess(:,:,1) = gdrpr;
hess(:,:,2) = -gdzz;
hess(:,:,3) = -gdrz;
hess(:,:,4) = gdrpz;
hess = hess*const;

% if anynan(hess)
%   error('Execution stopped: NaN value detected.');
% end

end
