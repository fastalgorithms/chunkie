function [sysmat] = buildmat(trap,kern,quadorder,opdims,type)
%CHNK.QUADBA.BUILDMAT build matrix for given kernel and trapper
% description of boundary using an Alpert-style rule
%
%  Input:
%   trap - trapper description of geometry
%   kern - kernel function of the form kern(src,targ,srctau,targtau)
%   quadorder - requested order for quadrature rule (0,2,4,8,16 available
%          for log type singularities. does nothing for smooth)
%   opdims - the dimensions of the kern operator (standard would be [1 1]
%         but matrix valued kernels, like 2d stokes, would be [2 2])
%   type - singularity type. only 'smooth' and 'log' are implemented
%   
%  Output:
%   sysmat - system matrix


if strcmpi(type,'smooth')
    % apply smooth rule immediately and exit
    r = trap.r;
    d = trap.d;
    d2 = trap.d2;
    n = trap.n;

    wts = trap.wts;
    wts2 = repmat(wts,opdims(2),1);
    wts2 = wts2(:);
    
    srcinfo = []; srcinfo.r = r; srcinfo.d = d; srcinfo.d2 = d2;
    srcinfo.n = n;
    sysmat = kern(srcinfo,srcinfo);
    sysmat = bsxfun(@times,sysmat,(wts2(:)).');
    return
elseif strcmpi(type,'log')
    % get log singularity quadrature nodes and weights
    qavail = chnk.quadba.logavail();
    [~,i] = min(abs(qavail-quadorder));
    if (qavail(i) ~= quadorder)
        warning('order %d not found, using order %d', ...
            quadorder,qavail(i));
        quadorder = qavail(i);
    end    
    [xs,ws,nskip] = chnk.quadba.getlogquad(quadorder);
else
    error('quadrature type "%s" not available',type)
end

left = min(floor(xs)); right = max(ceil(xs));
ninterphalf = floor(quadorder/2)+1;
ninterp = ninterphalf*2;
indrel = (left-ninterphalf):(right+ninterphalf);
assert(trap.npt >= 2*nskip-1,'too few points on curve to use this rule');

% find interpolation indices (relative for now) for each
% support node (use the nuse closest points)

nuse = ninterp;
induse = zeros(nuse,length(xs));

for i = 1:length(xs)
	 [~,indtemp] = sort(abs(indrel-xs(i)));
     induse(:,i) = sort(indtemp(1:nuse));
end

% get interpolation coefficients for each support node

interp_coeffs = zeros(length(indrel),length(xs));
for i = 1:length(xs)
    indusei = induse(:,i);
    cfs = barycoefs(indrel(indusei),xs(i));
    interp_coeffs(indusei,i) = cfs;
end

% grab curve data

r = trap.r;
d = trap.d;
d2 = trap.d2;
ds = sqrt(sum(d.^2,1));
n = trap.n;

%% use conv routine to get interpolated values of r, d, ds

% first pad r, d, ds

npt = trap.npt;
endpad = (npt-(length(indrel)-1)/2+1):npt;
startpad = 1:(length(indrel)-1)/2;
rpad = [r(:,endpad), r, r(:,startpad)];
npad = [n(:,endpad), n, n(:,startpad)];
dpad = [d(:,endpad), d, d(:,startpad)];
d2pad = [d2(:,endpad), d2, d2(:,startpad)];
dspad = [ds(endpad), ds, ds(startpad)];

rpad = rpad.';
npad = npad.';
dpad = dpad.';
d2pad = d2pad.';
dspad = dspad(:);

% interpolated values would be from a sort of transposed convolution,
% so flip coeffs here

tmp = fliplr(interp_coeffs);
dim = trap.dim;
rinterp = zeros(npt,dim,length(xs));
ninterp = zeros(npt,dim,length(xs));
dinterp = zeros(npt,dim,length(xs));
d2interp = zeros(npt,dim,length(xs));
dsinterp = zeros(npt,length(xs));

for i = 1:length(xs)
    interpi = tmp(:,i);
    rinterp(:,:,i) = conv2(rpad,interpi,'valid');
    ninterp(:,:,i) = conv2(npad,interpi,'valid');
    dinterp(:,:,i) = conv2(dpad,interpi,'valid');
    d2interp(:,:,i) = conv2(d2pad,interpi,'valid');    
    dsinterp(:,i) = conv2(dspad,interpi,'valid');
end

rinterp = permute(rinterp,[2 3 1]);
ninterp = permute(ninterp,[2 3 1]);
dinterp = permute(dinterp,[2 3 1]);
d2interp = permute(d2interp,[2 3 1]);
dsinterp = dsinterp.';

%% pretend smooth works

r = trap.r;
d = trap.d;
n = trap.n;
d2 = trap.d2;

wts = trap.wts;
wts2 = repmat((wts(:)).',opdims(2),1);
wts2 = wts2(:);

srcinfo = []; srcinfo.r = r; srcinfo.d = d; srcinfo.d2 = d2;
srcinfo.n = n;
sysmat = bsxfun(@times,kern(srcinfo,srcinfo),(wts2(:)).');

%% kill elements too close to diagonal (cyclically speaking)

ind2kill = bsxfun(@plus,(1:npt).',-(nskip-1):(nskip-1));
ind2kill = indwrap(ind2kill,npt);
linind2kill = (ind2kill-1)*npt+repmat( (1:npt).',1,2*nskip-1);
linind2kill = sort(linind2kill(:));
sysmat(linind2kill) = 0;

%% compute corrected values in a loop

ifix = zeros(size(interp_coeffs,1)*npt,1);
jfix = ifix;
vfix = ifix;

ihave = 0;

interp_coeffs = interp_coeffs.';
wh = (ws(:)).'*trap.h;

if quadorder == 0
    fixmat = sparse(npt*opdims(1),npt*opdims(2));
else
    srcinfo = [];
    targinfo = [];    
    for i = 1:npt

        % spatial locations of support nodes and scaled weights
        ri = rinterp(:,:,i);
        ni = ninterp(:,:,i);
        di = dinterp(:,:,i);
        d2i = d2interp(:,:,i);
        ds_i = dsinterp(:,i);
        
        
        wh2 = (wh(:)).*(ds_i(:));
        wh2 = repmat((wh2(:)).',opdims(2),1);
        wh2 = (wh2(:)).';

        % eval kernel and combine
        srcinfo.r = ri; srcinfo.d = di; srcinfo.d2 = d2i;
        srcinfo.n = ni;
        targinfo.r = r(:,i); targinfo.d = d(:,i); 
        targinfo.d2 = d2(:,i); targinfo.n = n(:,i);
        kmat = kern(srcinfo,targinfo).*wh2;
        mattemp = sparse(kmat*interp_coeffs);
        [itemp,jtemp,vtemp] = find(mattemp);
        ntemp = nnz(mattemp);
        
        % insert indices in the right place
        istart = ihave+1;
        iend = ihave+ntemp;
        ifix(istart:iend) = (i-1)*opdims(1)+itemp;
        jpts = double(idivide(int64(jtemp-1),int64(opdims(2)))+1);
        jdiff = jtemp-(jpts-1)*opdims(2);
        jpts = indwrap(indrel(jpts)+i,npt);
        jfix(istart:iend) = (jpts(:)-1)*opdims(2) + jdiff(:);
        vfix(istart:iend) = vtemp;
        ihave = ihave+ntemp;
    end
    ifix = ifix(1:ihave);
    jfix = jfix(1:ihave);
    vfix = vfix(1:ihave);

    fixmat = sparse(ifix,jfix,vfix);
end

sysmat = sysmat + fixmat;


end

function cfs = barycoefs(pts,targ)

% if targ is within 10 eps of a point, just use that value
icheck = (abs(pts-targ) < 10*eps);
if any(icheck)
    ii = find(icheck);
    ii = ii(1);
    cfs = zeros(length(pts),1);
    cfs(ii) = 1.0;
    return
end

% weights
wbig = 1./(pts(:)- (pts(:)).');
wbig(eye(length(pts),'logical')) = 1;
w = prod(wbig,1);

% form coefs
cfs = w(:)./(targ-pts(:));
cfs = cfs/(sum(cfs));

end



function indout = indwrap(indin,n)

indout = indin;
indout(indout < 1) = n + indout(indout < 1);
indout(indout > n) = indout(indout > n) - n;

end
