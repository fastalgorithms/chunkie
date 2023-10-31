function chnkr = roundedpolygon(verts,cparams,pref,edgevals)
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
%    	cparams.ifclosed = true, if it's a closed polygon
%                          false, if it's an open segment (true)
%       cparams may also specify any of the values to send to the
%            refine routine, otherwise it is called with defaults
%
%           ~ Rounding parameters ~
% 	    cparams.eps - resolve curve to tolerance eps
%                    resolve coordinates, arclength,
%          	     and first and second derivs of coordinates
%		         to this tolerance (1.0e-6)
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
%   chnkr = roundedpolygon(barbell_verts); % rounded "barbell" domain
%                                     % with standard options
%   cparams = []; cparams.rounded = false;
%   pref = []; pref.k = 30;
%   chnkr = chunkpoly(barbell_verts,cparams,pref); % not rounded
%
% See also CHUNKERFUNC, CHUNKER, CHUNKERPREF, REFINE

[dimv,nv] = size(verts);

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

autowidths = false;
autowidthsfac = 0.1;
ifclosed = true;
eps = 1e-6;
nover = 0;

if isfield(cparams,'ifclosed')
   ifclosed = cparams.ifclosed;
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

if (autowidths || widths_not_set)
   widths = autowidthsfac*...
   	  min(edges2(1:end-1),edges2(2:end));
end

if ifclosed
    widths = [widths(:); widths(1)];
    verts = [verts, verts(:,1)];
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
chnkr = chnkr.makedatarows(nvals);
k = chnkr.k; dim = chnkr.dim;
[t,w] = lege.exps(k);
onek = ones(k,1);


if nvals > 0
    aint = lege.intmat(chnkr.k);
end

for i = 1:nv-1
    % grab vertices
    r1 = verts(:,i); r2 = verts(:,i+1);
    w1 = widths(i); w2 = widths(i+1);
    l = sqrt(sum((r1-r2).^2));
    assert(l > w1+w2+2*eps(1)*l,'widths too large for side');
    
    % make chunk in middle
    v = (r2-r1)/l;
    ts = w1+(l-w2-w1)*(t+1)/2.0;
    chnkr = chnkr.addchunk();
    nch = chnkr.nch;
    if nvals > 0
        val1 = edgevals(:,i);
        chnkr.data(:,:,nch) = bsxfun(@times,val1,onek(:).');
    end
    chnkr.r(:,:,nch) = r1 + bsxfun(@times,v(:),(ts(:)).');
    chnkr.d(:,:,nch) = repmat(v(:),1,k);
    chnkr.d2(:,:,nch) = zeros(dim,k);
    chnkr.adj(1,nch) = 0; chnkr.adj(2,nch) = 0;
    if nch > 1
        chnkr.adj(1,nch) = nch-1;
        chnkr.adj(2,nch-1) = nch;
    end
    chnkr.h(nch) = (l-w2-w1)/2.0;
    
    if or(i < nv-1,ifclosed)
			% chunk up smoothed corner made by three verts
      if (i==nv-1)
        if nvals > 0
          val2 = edgevals(:,1);
        end
        r3 = verts(:,2);
      else
        if nvals > 0
          val2 = edgevals(:,i+1);
        end
        r3 = verts(:,i+2);
      end
      l2 = sqrt(sum((r2-r3).^2));
      v = -v;
      v2 = (r3-r2)/l2;
      cosphi = dot(v2,v); %cosine of angle between edges
      sinphi = sqrt(1-cosphi*cosphi);
      trange = w2*sqrt((1-cosphi)/2); %parameter space range of gaussian
      hbell = abs(trange)/8.0; % width parameter of gaussian
      m = sinphi/(1-cosphi); % slope of abs approx
      cpt.ta = -trange;
      cpt.tb = trange;
      cpt.eps = eps; cpt.levrestr = 0; cpt.ifclosed = 0;
      chnkrt = sort(chunkerfunc(@(t)fround(t,m,hbell,dim),cpt,pref));
      
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
      chnkr.adj(:,nch+1:nch+ncht) = chnkrt.adj+nch;
      chnkr.adj(2,nch) = nch+1;
      chnkr.adj(1,nch+1) = nch;
      chnkr.h(nch+1:nch+ncht) = chnkrt.h;
      
      if nvals > 0
        ds = bsxfun(@times,reshape(sum((rotmat*chnkrt.d(:,:)).^2,1),k,ncht), ...
                    (chnkrt.h(:)).');
        dsw = bsxfun(@times,w(:),ds);
        dssums = sum(dsw,1);
        dssums2 = cumsum([0,dssums(1:end-1)]);
        dsint = aint*ds;
        dsint = bsxfun(@plus,dsint,dssums2);
        lencorner = sum(dsw(:));
        ss = -lencorner/2.0 + dsint;
        ss = ss/lencorner*16;
        erfss = erf(ss);
        datass = reshape((val2(:)-val1(:))/2*((erfss(:)).'+1) ...
			 +val1(:),nvals,k,ncht);
        chnkr.data(:,:,nch+1:nch+ncht) = datass;
      end
    end
    
end

if ifclosed
    nch = chnkr.nch;
    chnkr.adj(1,1) = nch;
    chnkr.adj(2,nch) = 1;
end

chnkr = chnkr.refine(cparams);

% compute normals

chnkr.n = normals(chnkr);
    
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
