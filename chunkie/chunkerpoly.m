function [chnkr, varargout] = chunkerpoly(verts,cparams,pref,edgevals)
%CHUNKPOLY return a chunker object corresponding to the polygon with the
% given vertices. By default, a polygon with rounded corners is returned.
% Open and higher dimensional "polygons" are allowed. Optional dyadic
% refinement near corners of true polygons (no rounding).
%
% Syntax: chnkr = chunkpoly(verts,cparams,pref,edgevals)
%
% Input:
%    verts - (dimv,nverts) array of "polygon" vertices
%            in order
%
% Optional input:
%    cparams - options structure
%       cparams.rounded = true if corner rounding is to be used.
%                         false if no rounding is used (true)
%       cparams.widths = radius around each corner point where either the
%                        rounding or dyadic refinement is applied 
%                        (defaults to 1/10th of minimum
%                         of length of adjoining edges)
%       cparams.autowidths = automatically compute widths (false)
%       cparams.autowidthsfac = if using autowidths, set widths
%                             to autowidthsfac*minimum of adjoining
%                             edges (0.1)
%       cparams.smoothwidths = use a smooth version of autowidth
%                              detection, to allow for shape
%                              derivative computation (false)
%    	  cparams.ifclosed = true, if it's a closed polygon
%                          false, if it's an open segment (true)
%
%           ~ Rounding parameters ~
% 	    cparams.eps - resolve curve to tolerance eps
%                    resolve coordinates, arclength,
%          	     and first and second derivs of coordinates
%		         to this tolerance (1.0e-6)
%
%           ~ Parameters without rounding ~
%       cparams.dyadic = do dyadic refinement into corners (true)
%       cparams.depth = depth of dyadic refinement (30)
%
%   pref - chunkerpref object/ preference structure 
%       pref.k determines order of underlying Legendre nodes
%   edgevals - very optional input. specifies constant values along each
%              edge. The routine then stores a smoothed interpolant 
%              of these edge values on the rounded structure in the 
%              output chunker's data field
%
% Output:
%   chnkr - chunker object corresponding to (rounded) polygon
%                  
% Examples:
%   barbell_verts = chnk.demo.barbell();
%   chnkr = chunkpoly(barbell_verts); % rounded "barbell" domain
%                                     % with standard options
%   cparams = []; cparams.rounded = false;
%   pref = []; pref.k = 30;
%   chnkr = chunkpoly(barbell_verts,cparams,pref); % not rounded
%
% See also CHUNKERFUNC, CHUNKER, CHUNKERPREF

[dimv,nv] = size(verts);
ilabel = 1:nv;
nvin = nv;

assert(dimv > 1)

if nargin < 2
    cparams = [];
end
if nargin < 3
    p = [];
    p.dim = dimv;
    pref = chunkerpref(p);
else
    pref = chunkerpref(pref);
    if pref.dim ~= dimv
        warning('dimensions dont match, overwriting with vertex dim');
        pref.dim = dimv;
    end
end
if nargin < 4
    edgevals = [];
end

dim = pref.dim;

rounded = true;
dyadic = true;
depth = 30;

autowidths = false;
autowidthsfac = 0.1;
smoothwidths = false;
ifclosed = true;
eps = 1e-6;
nover = 0;
widths = zeros(nv, 1);

if isfield(cparams,'ifclosed')
   ifclosed = cparams.ifclosed;
end
if isfield(cparams,'rounded')
   rounded = cparams.rounded;
end
if isfield(cparams,'eps')
   eps = cparams.eps;
end
if isfield(cparams,'nover')
   nover = cparams.nover;
end
   
if (ifclosed)
   verts2 = [verts(:,end), verts, verts(:,1)];
   edges2 = sqrt(sum(diff(verts2,1,2).^2,1));
else
   edges2 = [0, sqrt(sum(diff(verts,1,2).^2,1)), 0];
end

widths_not_set = true;

if isfield(cparams,'widths')
   widths = cparams.widths;
   widths_not_set = false;
end
if isfield(cparams,'autowidths')
   autowidths = cparams.autowidths;
end
if isfield(cparams,'autowidthsfac')
   autowidthsfac = cparams.autowidthsfac;
end

if isfield(cparams, 'smoothwidths')
   smoothwidths = true;
end

if isfield(cparams,'depth')
    depth = cparams.depth;
    assert(depth >= 0, 'depth must be a nonnegative integer');
end

% fields for using smoothwidths, unused if smoothwidths is false 
widthspre = zeros(dim, nv);
widthsdc = zeros(dim, nv);
widthsnext = zeros(dim, nv);


if (autowidths || widths_not_set)
   if ~smoothwidths
     widths = autowidthsfac*...
     	  min(edges2(1:end-1),edges2(2:end));
   else
      if (ifclosed)
         rl = verts(:,end); rc = verts(:,1); rr = verts(:,2);
         [w,dwdl,dwdc,dwdr] = smoothminwidth(rl,rc,rr,autowidthsfac);
         widths(1) = w;
         widthspre(:,1) = dwdl;
         widthsdc(:,1) = dwdc;
         widthsnext(:,1) = dwdr;
         rl = verts(:,end-1); rc = verts(:,end); rr = verts(:,1);
         [w,dwdl,dwdc,dwdr] = smoothminwidth(rl,rc,rr,autowidthsfac);
         widths(end) = w;
         widthspre(:,end) = dwdl;
         widthsdc(:,end) = dwdc;
         widthsnext(:,end) = dwdr;
      else
         rc = verts(:,1); rr = verts(:,2);
         w = sqrt(sum( (rc-rr).^2 ));
         widths(1) = autowidthsfac*w;
         widthspre(:,1) = zeros(dim,1);
         widthsdc(:,1) = (rc-rr)/w;   
         widthsnext(:,1) = (rr-rc)/w;
         rl = verts(:,end-1); rc = verts(:,end);
         w = sqrt(sum( (rc-rl).^2 ));
         widths(end) = autowidthsfac*w;
         widthspre(:,end) = (rl-rc)/w;
         widthsdc(:,end) = (rc-rl)/w;
         widthsnext(:,end) = zeros(dim,1);
      end
      for i = 2:(nv-1)
         [widths(i),widthspre(:,i),widthsdc(:,i),widthsnext(:,i)] = ...
            smoothminwidth(verts(:,i-1),verts(:,i),verts(:,i+1), ...
            autowidthsfac);
      end
   end
end

if ifclosed
    widths = [widths(:); widths(1)];
    widthspre = [widthspre, widthspre(:,1)];
    widthsdc = [widthsdc, widthsdc(:,1)];
    widthsnext = [widthsnext, widthsnext(:,1)];
    verts = [verts, verts(:,1)];
    ilabel = [ilabel,1];
    nv = size(verts,2);
end


nedge = nv - 1;
nvals = 0;
if numel(edgevals) ~= 0
    assert(rem(numel(edgevals),nedge) == 0, ...
        'number of edge values should be multiple of number of edges');
    nvals = numel(edgevals)/nedge;
    edgevals = reshape(edgevals,nvals,nedge);
end

chnkr = chunker(pref);

nstorgrad = 0;
k = chnkr.k; dim = chnkr.dim;
if (smoothwidths), nstorgrad = nvin*(dim^2 + dim); end
ndata = nvals + nstorgrad;
igraddatarows = (nvals + 1):(nvals + nstorgrad);
chnkr = chnkr.makedatarows(ndata);
[t,w] = lege.exps(k);
onek = ones(k,1);

if smoothwidths
  itmp1 = 0:(nvin-1);
  itmp2 = 1:nvin;
  igrad = [nvals+1; nvals] + dim^2*[itmp1; itmp2];
  isgrad = [nvals + dim^2*nvin + 1; nvals + dim^2*nvin] + dim*[itmp1; itmp2];
end
dim2 = dim^2;
onedim = eye(dim); onedim = onedim(:);

% pre-define graded mesh for adaptivity in corners

if ~rounded
    tadap = zeros(k*(depth+1),1);
    hadap = zeros(depth+1,1);
    istart = 1;
    iend = k;
    t2 = 0.5*(t + 1);
    scale = 2^(-depth);
    shift = 0;
    tadap(istart:iend) = shift + t2*scale;
    hadap(1) = scale;
    shift = scale;
    for i = 1:depth
        istart = istart+k;
        iend = iend+k;
        tadap(istart:iend) = shift + t2*scale;
        hadap(i+1) = scale;
        shift = shift + scale;
        scale = scale*2;
    end
    tadap = reshape(tadap,1,k,depth+1);
end

if nvals > 0
    aint = lege.intmat(chnkr.k);
end

for i = 1:nv-1
    % grab vertices
    r1 = verts(:,i); r2 = verts(:,i+1);
    w1 = widths(i); w2 = widths(i+1);
    l = sqrt(sum((r1-r2).^2));
    assert(l > w1+w2+2*eps(1)*l,'widths too large for side');

    if smoothwidths
      dw2dr1 = widthspre(:,i+1); dw2dr2 = widthsdc(:,i+1);
      dw2dr3 = widthsnext(:,i+1);
    end
    
    % make chunk in middle
    v = (r2-r1);
    ts = (w1+(l-w2-w1)*(t+1)/2.0)/l;
    chnkr = chnkr.addchunk();
    nch = chnkr.nch;
    
    if nvals > 0
        val1 = edgevals(:,i);
        chnkr.data(:,:,nch) = bsxfun(@times,val1,onek(:).');
    end
    h = (l-w2-w1)/2.0/l;
    chnkr.r(:,:,nch) = r1 + bsxfun(@times,v(:),(ts(:)).');
    chnkr.d(:,:,nch) = repmat(v(:),1,k)*h;
    chnkr.d2(:,:,nch) = zeros(dim,k)*h*h;
    chnkr.adj(1,nch) = 0; chnkr.adj(2,nch) = 0;
    if nch > 1
        chnkr.adj(1,nch) = nch-1;
        chnkr.adj(2,nch-1) = nch;
    end
   

    if smoothwidths
      % derivative with respect to left vertex/right vertex
      ilgrad = igrad(1,ilabel(i)):igrad(2,ilabel(i)); 
      irgrad = igrad(1,ilabel(i+1)):igrad(2,ilabel(i+1));
      chnkr.data(ilgrad,:,nch) = bsxfun(@times,onedim(:),1-ts(:).');
      chnkr.data(irgrad,:,nch) = bsxfun(@times,onedim(:),ts(:).');       
    
      dsdt = sqrt(sum( (chnkr.d(:,:,nch)).^2 ,1));
      tau = bsxfun(@rdivide,chnkr.d(:,:,nch),dsdt);
      ilsgrad = isgrad(1,ilabel(i)):isgrad(2,ilabel(i)); 
      irsgrad = isgrad(1,ilabel(i+1)):isgrad(2,ilabel(i+1));
      chnkr.data(ilsgrad,:,nch) = -tau/norm(v(:));
      chnkr.data(irsgrad,:,nch) = tau/norm(v(:));    
    end
    

    if or(i < nv-1, ifclosed)
        if i == nv-1
            if nvals > 0
                val2 = edgevals(:,1);
            end
            r3 = verts(:,2);
            i1 = ilabel(i); i2 = ilabel(i+1); i3 = ilabel(2);
        else
            if nvals > 0
                val2 = edgevals(:,i+1);
            end
            r3 = verts(:,i+2);
            i1 = ilabel(i); i2 = ilabel(i+1); i3 = ilabel(i+2);
        end
        if rounded

            if ~smoothwidths
                
               % chunk up smoothed corner made by three verts
               l2 = sqrt(sum((r2-r3).^2));
               v = -v/l;
               v2 = (r3-r2)/l2;
               cosphi = dot(v2,v); %cosine of angle between edges
               sinphi = sqrt(1-cosphi*cosphi);
               trange = w2*sqrt((1-cosphi)/2); %parameter space range of gaussian
               hbell = abs(trange)/8.0; % width parameter of gaussian
               m = sinphi/(1-cosphi); % slope of abs approx
               cpt.ta = -trange;
               cpt.tb = trange;
               cpt.eps = eps; cpt.levrestr = 0; cpt.ifclosed = 0;
               chnkrt = sort(chunkerfunc(@(t) fround(t,m,hbell,dim), cpt, pref));

               % do optimal procrustes match of left and right ends
               rl = fround(-trange,m,hbell,dim); rr = fround(trange,m,hbell,dim);
               [um,~,vm] = svd([w2*v w2*v2]*([rl rr].'));
               rotmat = um*vm.';
               

              % copy in rotated and translated chunks
               ncht = chnkrt.nch;
               chnkr = chnkr.addchunk(ncht);
               chnkr.r(:,:,nch+1:nch+ncht) = reshape(rotmat*chnkrt.r(:,:),dim,k,ncht)+r2;
               chnkr.d(:,:,nch+1:nch+ncht) = reshape(rotmat*chnkrt.d(:,:),dim,k,ncht);
               chnkr.d2(:,:,nch+1:nch+ncht) = reshape(rotmat*chnkrt.d2(:,:),dim,k,ncht);
               chnkr.adj(:,(nch+1):(nch+ncht)) = chnkrt.adj+nch;
               chnkr.adj(2,nch) = nch+1;
               chnkr.adj(1,nch+1) = nch;    
            else
               l2 = sqrt(sum((r2-r3).^2));
               l = sqrt(sum((r2-r1).^2));
            
               sig = 1/8;

               cpt.ta=-1; cpt.tb = 0; cpt.ifclosed = 0;
               fun = @(t) froundnew(t, r1, r2, r3, w2, sig);
               [chnkrt1, tab1] = chunkerfunc(fun, cpt, pref);
               cpt.ta=0; cpt.tb = 1; cpt.ifclosed = 0;
               [chnkrt2, tab2] = chunkerfunc(fun, cpt, pref);
            
               % copy in rounded corner chunks
               scal1=l/w2;
               dscal1r1 = (r1-r2)/(l*w2) - dw2dr1*l/(w2^2);
               dscal1r2 = (r2-r1)/(l*w2) - dw2dr2*l/(w2^2);
               dscal1r3 = -dw2dr3*l/(w2^2);
            
               ncht1 = chnkrt1.nch;
               chnkr = chnkr.addchunk(ncht1);

               iind1 = nch+1:(nch+ncht1);
               chnkr.r(:,:,iind1) = reshape(chnkrt1.r(:,:), dim, k, ncht1);
               chnkr.d(:,:,iind1) = reshape(chnkrt1.d(:,:), dim, k, ncht1);
               chnkr.d2(:,:,iind1) = reshape(chnkrt1.d2(:,:), dim, k, ncht1);
               chnkr.adj(:,iind1) = chnkrt1.adj+nch;
               chnkr.adj(2,nch) = nch+1;
               chnkr.adj(1,nch+1) = nch;

            % derivative with respect to left/center/right vertex
               tused1 = tab1(1,:) + (t(:)+1)/2.*(tab1(2,:) - tab1(1,:)); 
               tused1 = tused1(:);

               [~, dtmp, d2tmp] = froundnew(tused1, r1, r2, r3, w2, sig);
            
               [drr1, drr2, drr3, dsr1, dsr2, dsr3] = ...
                    froundnewgrad(tused1, r1, r2, r3, w2, ...
                                  dw2dr1, dw2dr2, dw2dr3, sig);            
            
               treal1 = tused1/scal1;
               treal1 = repmat(treal1(:).', dim, 1); 
               treal1 = treal1(:);
               ddrr1 = reshape(dscal1r1*( (treal1(:).*dtmp(:)).'), dim, dim, k*ncht1);
               ddrr2 = reshape(dscal1r2*( (treal1(:).*dtmp(:)).'), dim, dim, k*ncht1);
               ddrr3 = reshape(dscal1r3*( (treal1(:).*dtmp(:)).'), dim, dim, k*ncht1);

               drr1all = drr1 + ddrr1;
               drr2all = drr2 + ddrr2;
               drr3all = drr3 + ddrr3;
            
               i1grad = igrad(1,i1):igrad(2,i1); 
               i2grad = igrad(1,i2):igrad(2,i2);
               i3grad = igrad(1,i3):igrad(2,i3);
               chnkr.data(i1grad,:,iind1) = reshape(drr1all, dim^2, k, ncht1);
               chnkr.data(i2grad,:,iind1) = reshape(drr2all, dim^2, k, ncht1);
               chnkr.data(i3grad,:,iind1) = reshape(drr3all, dim^2, k, ncht1);
            
               i1sgrad = isgrad(1,i1):isgrad(2,i1);
               i2sgrad = isgrad(1,i2):isgrad(2,i2);
               i3sgrad = isgrad(1,i3):isgrad(2,i3);
            
               ds = sqrt(sum(dtmp.^2,1));
            
               dd2 = sum(dtmp.*d2tmp,1)./ds.*(tused1(:).')/scal1;
               dsr1all = dsr1*scal1 + dscal1r1*ds + dscal1r1*dd2*scal1;
               dsr2all = dsr2*scal1 + dscal1r2*ds + dscal1r2*dd2*scal1;
               dsr3all = dsr3*scal1 + dscal1r3*ds + dscal1r3*dd2*scal1;
               
               sc1 = 1./scal1./ds;
               chnkr.data(i1sgrad,:,iind1) = reshape(dsr1all.*sc1, dim, k, ncht1);
               chnkr.data(i2sgrad,:,iind1) = reshape(dsr2all.*sc1, dim, k, ncht1);
               chnkr.data(i3sgrad,:,iind1) = reshape(dsr3all.*sc1, dim, k, ncht1);

               % copy in rounded corner chunks
               scal2 = l2/w2;
               dscal2r1 = - dw2dr1*l2/(w2^2);
               dscal2r2 = (r2-r3)/(l2*w2) - dw2dr2*l2/(w2^2);
               dscal2r3 = (r3-r2)/(l2*w2) - dw2dr3*l2/(w2^2);
            
               ncht2 = chnkrt2.nch;
               chnkr = chnkr.addchunk(ncht2);
               iind2 = (nch+ncht1+1):(nch+ncht1+ncht2);
               chnkr.r(:,:,iind2) = reshape(chnkrt2.r(:,:), dim, k, ncht2);
               chnkr.d(:,:,iind2) = reshape(chnkrt2.d(:,:), dim, k, ncht2);
               chnkr.d2(:,:,iind2) = reshape(chnkrt2.d2(:,:), dim, k, ncht2);
               chnkr.adj(:,iind2) = chnkrt2.adj + nch + ncht1;
               chnkr.adj(2,nch+ncht1) = nch + ncht1 + 1;
               chnkr.adj(1,nch+ncht1+1) = nch + ncht1;

            % derivative with respect to left vertex/right vertex
            
               tused2 = tab2(1,:) + (t(:)+1)/2.*(tab2(2,:) - tab2(1,:)); 
               tused2 = tused2(:);

               [~, dtmp, d2tmp] = froundnew(tused2, r1, r2, r3, w2, sig);
            
               [drr1, drr2, drr3, dsr1, dsr2, dsr3] = ...
                    froundnewgrad(tused2, r1, r2, r3, w2, ...
                                  dw2dr1, dw2dr2, dw2dr3, sig);            
           
               treal2 = tused2/scal2;
               treal2 = repmat(treal2(:).', dim, 1); 
               treal2 = treal2(:);

               ddrr1 = reshape(dscal2r1*( (treal2(:).*dtmp(:)).'), dim, dim, k*ncht2);
               ddrr2 = reshape(dscal2r2*( (treal2(:).*dtmp(:)).'), dim, dim, k*ncht2);
               ddrr3 = reshape(dscal2r3*( (treal2(:).*dtmp(:)).'), dim, dim, k*ncht2);
            
               drr1all = drr1 + ddrr1;
               drr2all = drr2 + ddrr2;
               drr3all = drr3 + ddrr3;
                        
               chnkr.data(i1grad,:,iind2) = reshape(drr1all, dim^2, k, ncht2);
               chnkr.data(i2grad,:,iind2) = reshape(drr2all, dim^2, k, ncht2);
               chnkr.data(i3grad,:,iind2) = reshape(drr3all, dim^2, k, ncht2);
            
            
               ds = sqrt(sum(dtmp.^2,1));
            
               dd2 = sum(dtmp.*d2tmp,1)./ds.*(tused2(:).')/scal2;
               dsr1all = dsr1*scal2 + dscal2r1*ds + dscal2r1*dd2*scal2;
               dsr2all = dsr2*scal2 + dscal2r2*ds + dscal2r2*dd2*scal2;
               dsr3all = dsr3*scal2 + dscal2r3*ds + dscal2r3*dd2*scal2;

               sc2 = 1./scal2./ds;
               chnkr.data(i1sgrad,:,iind2) = reshape(dsr1all.*sc2, dim, k, ncht2);
               chnkr.data(i2sgrad,:,iind2) = reshape(dsr2all.*sc2, dim, k, ncht2);
               chnkr.data(i3sgrad,:,iind2) = reshape(dsr3all.*sc2, dim, k, ncht2);

               ncht = ncht1 + ncht2;

            end

            if nvals > 0
                ichind = (nch+1):(nch+ncht);
                ds = reshape(sum(chnkr.d(:,:,ichind).^2, 1), k, ncht);
                dsw = bsxfun(@times, w(:), ds);
                dssums = sum(dsw, 1);
                dssums2 = cumsum([0, dssums(1:end-1)]);
                dsint = aint*ds;
                dsint = bsxfun(@plus, dsint, dssums2);
                lencorner = sum(dsw(:));
                ss = -lencorner/2.0 + dsint;
                ss = ss/lencorner*16;
                erfss = erf(ss);
                datass = reshape((val2(:)-val1(:))/2*((erfss(:)).'+1) ...
                    + val1(:), nvals, k, ncht);
                chnkr.data(:,:,ichind) = datass;
            end
        else
            % chunk up corner made by three verts with true corners 
            % and adaptive refinement
            l2 = sqrt(sum((r2-r3).^2));
            v = -v;
            v2 = (r3-r2)/l2;
            
            % left piece
            
            ncht = depth + 1;
            chnkrt = chunker(pref);
            chnkrt = chnkrt.addchunk(ncht);
            hall = hadap*w2/2; hall = repmat(hall(:).',k,1);
            hall=hall(:).';
            chnkrt.r = bsxfun(@times,v,tadap)*w2 + r2;
            chnkrt.d(:,:) = repmat(v,1,k*ncht).*hall;
            chnkrt.d2(:,:) = zeros(dim,k*ncht).*(hall.^2);
            chnkrt.adj = [0, 1:(ncht-1); 2:ncht, 0];
            chnkrt = reverse(chnkrt);
            chnkrt = sort(chnkrt);
            
            % copy in new chunks
            chnkr = chnkr.addchunk(ncht);
            chnkr.r(:,:,nch+1:nch+ncht) = chnkrt.r;
            chnkr.d(:,:,nch+1:nch+ncht) = chnkrt.d;
            chnkr.d2(:,:,nch+1:nch+ncht) = chnkrt.d2;
            chnkr.adj(:,nch+1:nch+ncht) = chnkrt.adj+nch;
            chnkr.adj(2,nch) = nch+1;
            chnkr.adj(1,nch+1) = nch;
            chnkr.adj(2,nch+ncht) = 0;

            if nvals > 0
                chnkr.data(:,:,nch+1:nch+ncht) = repmat(val1,1,k,ncht);
            end
            
            % right piece
            
            nch = chnkr.nch;
            
            ncht = depth + 1;
            chnkrt = chunker(pref);
            chnkrt = chnkrt.addchunk(ncht);
            chnkrt.r = bsxfun(@times,v2,tadap)*w2 + r2;
            chnkrt.d(:,:) = repmat(v2,1,k*ncht).*hall;
            chnkrt.d2(:,:) = zeros(dim,k*ncht).*(hall.^2);
            
            chnkrt.adj = [0, 1:(ncht-1); 2:ncht, 0];

            chnkrt = sort(chnkrt);
            
            % copy in new chunks
            chnkr = chnkr.addchunk(ncht);
            chnkr.r(:,:,nch+1:nch+ncht) = chnkrt.r;
            chnkr.d(:,:,nch+1:nch+ncht) = chnkrt.d;
            chnkr.d2(:,:,nch+1:nch+ncht) = chnkrt.d2;
            chnkr.adj(:,nch+1:nch+ncht) = chnkrt.adj+nch;
            chnkr.adj(2,nch) = nch+1;
            chnkr.adj(1,nch+1) = nch;
            
            chnkr = chnkr.addvert([nch,nch+1]);
            
            if nvals > 0
                chnkr.data(:,:,nch+1:nch+ncht) = repmat(val2,1,k,ncht);
            end
        end            
    end
    
end

if ifclosed
    nch = chnkr.nch;
    chnkr.adj(1,1) = nch;
    chnkr.adj(2,nch) = 1;
end

% update normals
chnkr.nstor(:,:,1:nch) = normals(chnkr);

% update weights
chnkr.wtsstor(:,1:nch) = weights(chnkr);

%chnkr = chnkr.refine();
%chnkr = chnkr.sort();


if smoothwidths && nargout > 1
    varargout{1} = igraddatarows;
end

end



function [r,d,d2] = fround(t,m,h,dim)

[y,dy,d2y] = chnk.spcl.absconvgauss(t,m,0.0,h);

r = zeros(dim,length(t));
d = zeros(dim,length(t));
d2 = zeros(dim,length(t));

r(1,:) = t;
r(2,:) = y;
d(1,:) = 1.0;
d(2,:) = dy;
d2(1,:) = 0.0;
d2(2,:) = d2y;

end



function [r,d,d2] = froundnew(t,r1,r2,r3,w,h)
%
% Make the rounded corner of width w corresponding to
% the vertices r1, r2, r3 (with r2 in the middle)
%
% Input:
%
% t - parameter in [-1,1]
% r1 - first vertex
% r2 - middle vertex
% r3 - last vertex
% w - width of corner (first point should be 
%            r2 + w*(r1-r2)/|r1-r2|)
% h - std dev of Gaussian used in rounding
%      recommend h = 1/8 for ~16,15,14 digits 
%      in position,first,second derivative
%
% Output:
%
% r,d,d2 - position, derivative, second derivative
%           of corner rounded according to specified
%           Gaussian std dev and width of corner
%

dim = size(r1,1);
m1=1; b0=0;
[y,dy,d2y] = chnk.spcl.absconvgauss(t,m1,b0,h);

r12 = r1-r2; u12 = r12/norm(r12);
r32 = r3-r2; u32 = r32/norm(r32);

rhs = w*[u12, u32];
xx = zeros(dim,2);
xx(1:2,1:2) = [-1 1; 1 1];
amat = rhs/xx;

r = zeros(dim,length(t));
d = zeros(dim,length(t));
d2 = zeros(dim,length(t));

r(1,:) = t;
r(2,:) = y;
d(1,:) = 1.0;
d(2,:) = dy;
d2(1,:) = 0.0;
d2(2,:) = d2y;

r = r2+amat*r;
d = amat*d;
d2 = amat*d2;

end

function [drr1,drr2,drr3,dsr1,dsr2,dsr3] = froundnewgrad(t,r1,r2,r3,w,...
                                                     dwdr1,dwdr2,dwdr3,h)
%
% Get the derivatives of the position and arclength density 
% with respect to the vertices for a rounded corner.
%
% Input:
%
% t - parameter in [-1,1]
% r1 - first vertex
% r2 - middle vertex
% r3 - last vertex
% w - width of corner (first point should be 
%            r2 + w*(r1-r2)/|r1-r2|)
% dwdr1 - gradient of w w.r.t. r1
% dwdr2 - gradient of w w.r.t. r2
% dwdr3 - gradient of w w.r.t. r3
% h - std dev of Gaussian used in rounding
%      recommend h = 1/8 for ~16,15,14 digits 
%      in position,first,second derivative
%
% Output:
%
% drr1 - gradient of position with respect to r1
% drr2 - gradient of position with respect to r2
% drr3 - gradient of position with respect to r3
% dsr1 - gradient of arclength density with respect to r1
% dsr2 - gradient of arclength density with respect to r2
% dsr3 - gradient of arclength density with respect to r3
%
%

dim = size(r1,1);
m1=1; b0=0;

[y,dy,~] = chnk.spcl.absconvgauss(t,m1,b0,h);

r12 = r1-r2; nr12 = norm(r12); u12 = r12/nr12;
r32 = r3-r2; nr32 = norm(r32); u32 = r32/nr32;

xx = zeros(dim,2);
xx(1:2,1:2) = [-1 1; 1 1];
rhs = [u12, u32];
amat = w*(rhs/xx);

drr1 = zeros(dim,dim,length(t));
drr2 = zeros(dim,dim,length(t));
drr3 = zeros(dim,dim,length(t));
dsr1 = zeros(dim,length(t));
dsr2 = zeros(dim,length(t));
dsr3 = zeros(dim,length(t));

r = zeros(dim,length(t));
d = zeros(dim,length(t));

for i = 1:dim
    
    onei = zeros(dim,1); onei(i)=1;
    
    du12dr1i = onei/nr12 - r12(i)*r12/(nr12^3);
    du32dr1i = zeros(dim,1);
    dwdr1i = dwdr1(i);
    
    du12dr3i = zeros(dim,1);
    du32dr3i = onei/nr32 - r32(i)*r32/(nr32^3);
    dwdr3i = dwdr3(i);
    
    du12dr2i = -du12dr1i;
    du32dr2i = -du32dr3i;
    dwdr2i = dwdr2(i);
    
    rhsdr1i = w*[du12dr1i, du32dr1i] + dwdr1i*[u12, u32];
    rhsdr2i = w*[du12dr2i, du32dr2i] + dwdr2i*[u12, u32];
    rhsdr3i = w*[du12dr3i, du32dr3i] + dwdr3i*[u12, u32];    
    
    amatdr1i = rhsdr1i/xx;
    amatdr2i = rhsdr2i/xx;
    amatdr3i = rhsdr3i/xx;
    
    r(1,:) = t;
    r(2,:) = y;
    
    d(1,:) = 1.0;
    d(2,:) = dy;
    
    ad = amat*d;
    s = sqrt(sum(ad.^2,1));
        
    
    drr1(i,:,:) = reshape(amatdr1i*r,size(drr1(i,:,:)));
    drr2(i,:,:) = reshape(onei + amatdr2i*r,size(drr2(i,:,:)));
    drr3(i,:,:) = reshape(amatdr3i*r,size(drr2(i,:,:)));
    
    dsr1(i,:) = (sum(ad.*(amatdr1i*d),1))./s;
    dsr2(i,:) = (sum(ad.*(amatdr2i*d),1))./s;
    dsr3(i,:) = (sum(ad.*(amatdr3i*d),1))./s;
end


end

function [sm,dsmda,dsmdb] = smoothmin(a,b)
%SMOOTHMIN smooth minimum of two positive numbers
%
% smoothmin(a,b) = 1/log(exp(1/a)+exp(1/b)-1)
%
% has the properties:
%
% min(a,b)/2 < smoothmin(a,b) < min(a,b)
% smoothmin(a,b) < smoothmin(c,b) if a<c
% smoothmin(a,b) = smoothmin(b,a)

abmin = min(a,b);
e1 = exp(-1./abmin);
ea = exp(1./a-1./abmin);
eb = exp(1./b-1./abmin);

logsafe = 1./abmin+log(ea+eb-e1);

sm = 1./logsafe;

% the derivative with respect to a is 
% 
% e^(1/a)/( log^2(e^(1/a)+e^(1/b)-1) * a^2 * (e^(1/a)+e^(1/b)-1) )

dsmda = 1.0./(logsafe.^2.*a.^2.*(1+exp(1./b-1./a)-exp(-1./a)));
dsmdb = 1.0./(logsafe.^2.*b.^2.*(1+exp(1./a-1./b)-exp(-1./b)));

end

function [w,dwdr1,dwdr2,dwdr3] = smoothminwidth(r1,r2,r3,a)
%SMOOTHMINWIDTH uses the three vertices that  make a corner
% and the fraction of the minimum (smooth min) of the two sidlengths
% to determine (a) the width and (b) the gradients of the width
% with respect to the vertex locations
%

ll = sqrt(sum( (r1-r2).^2));
lr = sqrt(sum( (r3-r2).^2));

[sm,dsmdl,dsmdr] = smoothmin(ll,lr);

w = a*sm;
dwdr1 = a*dsmdl*(r1-r2)/ll;
dwdr3 = a*dsmdr*(r3-r2)/lr;
dwdr2 = a*(dsmdl*(r2-r1)/ll + dsmdr*(r2-r3)/lr);

end
